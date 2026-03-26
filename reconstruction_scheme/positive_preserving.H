// 保正性头文件
// 包含用于确保物理量保持正值的函数和工具

#pragma once
#include "equation.H"
#include "linear_algebra.H"
#include "physical_bound.H"
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

// 检查数组是否包含NaN值
// 输入：变量数组array
// 输出：是否包含NaN值的布尔值
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

// 实施凸组合调整插值后的变量值
// 输入：原始变量值W，高阶插值后的变量值W_interpolated，调整系数theta，起始组件索引scomp，组件的数量ncomp
// 功能：根据theta值调整插值后的变量值
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

// 检查变量数组是否包含无效值
// 输入：变量数组varArray，问题参数pp，保正函数结构体pqf
// 输出：是否包含无效值的布尔值
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

// 计算theta值
// 输入：守恒变量consStar，守恒变量consHLLC，问题参数pp，保正函数结构体pqf
// 输出：theta值
AMREX_GPU_HOST_DEVICE
amrex::Real
compute_theta(VarArray const& consStar, VarArray const& consHLLC, ProbParm const& pp, PositiveQuantityFunctions const& pqf);

// 应用保正限制器
// 输入：守恒变量值W，高阶插值后的守恒变量值W_interpolated，保正函数结构体pqf，问题参数pp，组件的数量ncomp
// 功能：确保插值后的变量值保持正值
AMREX_GPU_HOST_DEVICE void apply_positive_limiter(
    const VarArray& W, // 守恒变量值
    VarArray& W_interpolated, // 高阶插值后的守恒变量值
    const PositiveQuantityFunctions& pqf, // 保正函数结构体
    const ProbParm& pp, // 问题参数
    int ncomp // 组件的数量
);

// 计算侧面的theta值
// 输入：维度idim，组件数ncomp，原始变量iPrim，守恒变量iCons，通量iFlux，一阶通量fluxFirstOrder，乘数mult，问题参数pp，保正函数pqf
// 输出：theta值
AMREX_GPU_HOST_DEVICE
amrex::Real
compute_theta_for_side(int idim, int ncomp,
    const VarArray& iPrim, const VarArray& iCons,
    const VarArray& iFlux, const VarArray& fluxFirstOrder,
    const amrex::Real& mult,
    const ProbParm& pp,
    const PositiveQuantityFunctions& pqf);

// 应用保正通量限制器
// 输入：计算区域bx，守恒量数组consin，原始量数组primin，旧数据数组olddatain，组件数ncomp，右端项数组rhs，通量数组flux，非守恒通量数组ncflux，时间步长dt，Runge-Kutta阶段stage，问题参数pp，保正函数pqf，网格间距dx
// 功能：确保通量计算后物理量保持正值
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