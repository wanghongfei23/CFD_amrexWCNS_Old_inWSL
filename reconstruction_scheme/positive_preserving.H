/**
 * @file
 * @brief 两相高阶重构中的保正限制器接口。
 * @ingroup reconstruction
 */

#pragma once
#include "equation.H"
#include "linear_algebra.H"
#include "physical_bound.H"
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

/**
 * @brief 检查状态数组中是否存在 NaN。
 * @param array 待检查的状态数组。
 * @return 存在 NaN 时返回 true。
 */
AMREX_GPU_HOST_DEVICE
inline bool contains_nan(const VarArray& array)
{
    for (int i = 0; i < NVar; ++i) {
        if (std::isnan(array[i])) {
            return true;
        }
    }
    return false;
}

/**
 * @brief 使用凸组合将高阶重构结果拉回到物理可接受区间。
 * @param W 单元中心参考状态。
 * @param W_interpolated 待修正的高阶重构状态。
 * @param theta 凸组合系数。
 * @param scomp 起始分量编号。
 * @param ncomp 参与修正的分量个数。
 */
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE void apply_convex_combination(
    const VarArray& W, // 原始变量值
    VarArray& W_interpolated, // 高阶插值后的变量值
    amrex::Real theta, // 调整系数
    int scomp, // 起始组件索引
    int ncomp // 组件的数量
)
{
    if (theta == 0.0 || contains_nan(W_interpolated)) {
        for (int i = scomp; i < ncomp; ++i) {
            W_interpolated[i] = W[i];
        }
        return;
    }
    for (int i = scomp; i < ncomp; ++i) {
        W_interpolated[i] = (1.0 - theta) * W[i] + theta * W_interpolated[i];
    }
}

/**
 * @brief 检查状态是否违反保正约束。
 * @param varArray 待检查的状态数组。
 * @param pp 两相物性参数。
 * @param pqf 保正目标函数集合。
 * @return 任一约束失败时返回 true。
 */
AMREX_INLINE AMREX_GPU_HOST_DEVICE bool containsInvalidValues(const VarArray& varArray,
    const ProbParm& pp,
    const PositiveQuantityFunctions& pqf)
{
    auto pq = pqf(varArray, pp); // 计算保正变量值

    for (int ifunc = 0; ifunc < pq.size(); ++ifunc) {
        amrex::Real value = pq[ifunc]; // 获取保正变量值

        if (std::isnan(value) || !std::isfinite(value))
            return true; // 检查是否为NaN或无穷大
        if (value < pqf.eps[ifunc])
            return true; // 检查是否小于阈值
    }

    return false;
}

/**
 * @brief 计算单个状态的保正缩放系数。
 * @param consStar 高阶候选状态。
 * @param consHLLC 一阶或低阶参考状态。
 * @param pp 两相物性参数。
 * @param pqf 保正目标函数集合。
 * @return 限幅系数 `theta`。
 */
AMREX_GPU_HOST_DEVICE
amrex::Real
compute_theta(VarArray const& consStar, VarArray const& consHLLC, ProbParm const& pp, PositiveQuantityFunctions const& pqf);

/**
 * @brief 对重构后的守恒变量应用保正限制器。
 * @param W 单元中心守恒变量。
 * @param W_interpolated 重构后的守恒变量。
 * @param pqf 保正目标函数集合。
 * @param pp 两相物性参数。
 * @param ncomp 处理的状态分量数。
 */
AMREX_GPU_HOST_DEVICE void apply_positive_limiter(
    const VarArray& W, // 守恒变量值
    VarArray& W_interpolated, // 高阶插值后的守恒变量值
    const PositiveQuantityFunctions& pqf, // 保正函数结构体
    const ProbParm& pp, // 问题参数
    int ncomp // 组件的数量
);

/**
 * @brief 为界面通量限制计算局部限幅系数。
 * @param idim 空间方向编号。
 * @param ncomp 状态分量数。
 * @param iPrim 单元原始变量。
 * @param iCons 单元守恒变量。
 * @param iFlux 高阶界面通量。
 * @param fluxFirstOrder 一阶参考通量。
 * @param mult 时间步与网格尺度组合出的权重。
 * @param pp 两相物性参数。
 * @param pqf 保正目标函数集合。
 * @return 限幅系数 `theta`。
 */
AMREX_GPU_HOST_DEVICE
amrex::Real
compute_theta_for_side(int idim, int ncomp,
    const VarArray& iPrim, const VarArray& iCons,
    const VarArray& iFlux, const VarArray& fluxFirstOrder,
    const amrex::Real& mult,
    const ProbParm& pp,
    const PositiveQuantityFunctions& pqf);

/**
 * @brief 在通量差分阶段施加保正通量限制。
 * @param bx 当前计算盒子。
 * @param consin 守恒变量输入。
 * @param primin 原始变量输入。
 * @param olddatain 上一时间层数据。
 * @param ncomp 状态分量数。
 * @param rhs 右端项数组。
 * @param flux 守恒通量数组。
 * @param ncflux 非守恒通量数组。
 * @param dt 当前子步长。
 * @param stage Runge-Kutta 阶段。
 * @param pp 两相物性参数。
 * @param pqf 保正目标函数集合。
 * @param dx 网格间距。
 */
AMREX_GPU_HOST_DEVICE void apply_positivive_preserving_flux_limiter(
    const amrex::Box& bx, // 计算区域
    const amrex::FArrayBox& consin, // 守恒量数组
    const amrex::FArrayBox& primin, // 原始量数组
    const amrex::FArrayBox& olddatain,
    const int ncomp,
    amrex::FArrayBox& rhs, // 右端项数组
    const amrex::GpuArray<amrex::FArrayBox*, BL_SPACEDIM>& flux, // 通量数组
    const amrex::GpuArray<amrex::FArrayBox*, BL_SPACEDIM>& ncflux, // 通量数组
    const amrex::Real dt, // 时间步长
    const int stage, // Runge-Kutta 阶段
    const ProbParm& pp, // 问题参数
    const PositiveQuantityFunctions& pqf, // 正限制函数
    const amrex::GpuArray<amrex::Real, BL_SPACEDIM>& dx // 网格间距
);