/**
 * @file
 * @brief 两相 stiffened-gas 模型的物性参数与状态变量编号定义。
 * @ingroup two_phase_model
 */

#pragma once

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

// using namespace amrex::literals;

/**
 * @enum CaseType
 * @brief 算例类型枚举，用于选择不同的初始化算例。
 */
enum CaseType {
    CASE_1D_BASE = 0,          // 1D 基础两相流配置
    CASE_1D_SHOCK_TUBE,         // 1D 两相激波管
    CASE_1D_SHOCK_INTERFACE,    // 1D 激波与静止界面相互作用
    CASE_2D_BASE,               // 2D 基础水柱配置
    CASE_2D_MACH2_4_SHOCK,      // 2D Mach 2.4 激波撞击水柱
    CASE_2D_MACH10_SHOCK,       // 2D Mach 10 激波撞击水柱
    CASE_2D_RMI_SF6_CASE1,      // 2D RMI_SF6 case1-半边截断
    CASE_2D_RMI_SF6_SQUARE,     // 2D RMI_SF6 方块
    CASE_3D_RMI_SF6             // 3D Richtmyer–Meshkov不稳定性
};

/**
 * @brief 两相 stiffened-gas 状态方程参数。
 */
struct ProbParm {
    // 算例类型选择
    CaseType case_type = CASE_2D_RMI_SF6_CASE1; 

    // 两相的比热比
    // amrex::Real gamma[2] = { 6.12, 1.4 }; // 水，空气
    amrex::Real gamma[2] = { 1.4, 1.093 };  // 空气，SF6
    // amrex::Real gamma[2] = { 1.093, 1.4 };  // SF6，空气

    // 两相的参考压力
    // amrex::Real pInf[2] = { 3.43e8, 0.0 }; // 水，空气
    amrex::Real pInf[2] = { 0.0, 0.0 }; // 空气，SF6/SF6，空气

    // 体积分数的小量阈值
    amrex::Real alphaEps = 1.1e-6;

    /**
     * @brief 根据体积分数计算混合区的等效 $\gamma$ 和 $p^\infty$。
     * @param alpha1 第一相体积分数。
     * @return 返回 `{gamma_bar, pInf_bar}`。
     */
    AMREX_FORCE_INLINE
    AMREX_GPU_HOST_DEVICE
    amrex::GpuArray<amrex::Real, 2> cal_gamma(amrex::Real alpha1) const
    {
        amrex::Real alpha2 = 1.0 - alpha1; // 第二相体积分数

        // 混合区
        amrex::Real gamma1 = gamma[0], gamma2 = gamma[1]; // 两相的比热比
        amrex::Real pInf1 = pInf[0], pInf2 = pInf[1]; // 两相的参考压力

        // 提前计算 gamma1 - 1 和 gamma2 - 1
        amrex::Real gamma1_minus_1 = gamma1 - 1;
        amrex::Real gamma2_minus_1 = gamma2 - 1;

        // 计算平均比热比 \overline{\gamma}
        amrex::Real inv_gamma_minus_1 = alpha1 / gamma1_minus_1 + alpha2 / gamma2_minus_1;
        amrex::Real gamma_bar = 1.0 / inv_gamma_minus_1 + 1;

        // 计算平均参考压力 \overline{p^\infty}
        amrex::Real gamma_pInf_bar = alpha1 * gamma1 * pInf1 / gamma1_minus_1 + alpha2 * gamma2 * pInf2 / gamma2_minus_1;
        amrex::Real pInf_bar = gamma_pInf_bar * (gamma_bar - 1) / gamma_bar;

        return { gamma_bar, pInf_bar };
    }
};

/**
 * @enum Variables
 * @brief 两相求解器状态向量中的分量编号。
 */
enum Variables {
    AlphaRho1 = 0, // 第一相密度
    AlphaRho2, // 第二相密度
    AMREX_D_DECL(XMom, YMom, ZMom), // 动量分量
    Energy, // 能量
    Alpha1, // 第一相体积分数
    NVar // 变量总数
};

/// 状态向量类型别名。
using VarArray = amrex::GpuArray<amrex::Real, NVar>;
/// 特征矩阵类型别名。
using VarMatrix = amrex::GpuArray<amrex::Real, NVar * NVar>;

// 定义两相流标志
#define TWO_PHASE 1