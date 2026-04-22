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
    // 1D
    Case_1D_CenteredWater_Air_Moving = 0,          // 中心水区域，整体运动（水/空气）
    Case_1D_ShockTube_Water_Air,                   // 水-空气激波管
    Case_1D_ShockInterface_Water_Air,              // 激波与水-空气界面作用
    // 2D
    Case_2D_WaterColumn_Air_Stationary,            // 静止水柱配置（水/空气）
    Case_2D_ShockImpact_Water_Air_Mach2_4,         // 马赫2.4激波撞击水柱（水/空气）
    Case_2D_ShockImpact_Water_Air_Mach10,          // 马赫10激波撞击水柱（水/空气）
    Case_2D_RMI_Air_SF6_Perturbed,                 // 扰动界面的空气-SF6 RMI
    Case_2D_RMI_SF6_Air_Square,                    // 方形SF6-空气 RMI
    // 3D
    Case_3D_RMI_Air_SF6_Perturbed                  // 三维扰动界面的空气-SF6 RMI
};

/**
 * @brief 两相 stiffened-gas 状态方程参数。
 */
struct ProbParm {
    // 算例类型选择
    CaseType case_type = Case_2D_RMI_Air_SF6_Perturbed; 

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