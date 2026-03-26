// 问题参数头文件（两相流）
// 包含两相流问题相关的参数和类型定义

#pragma once

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

// using namespace amrex::literals;

// 定义包含问题特定变量的结构体
struct ProbParm {
    // 两相的比热比
    amrex::Real gamma[2] = { 6.12, 1.4 };

    // 两相的参考压力
    amrex::Real pInf[2] = { 3.43e8, 0.0 };

    // 体积分数的小量阈值
    amrex::Real alphaEps = 1.1e-6;

    // 计算平均比热比和参考压力
    // 输入：第一相体积分数alpha1
    // 输出：包含平均比热比和参考压力的数组
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

// 变量索引枚举
enum Variables {
    AlphaRho1 = 0, // 第一相密度
    AlphaRho2, // 第二相密度
    AMREX_D_DECL(XMom, YMom, ZMom), // 动量分量
    Energy, // 能量
    Alpha1, // 第一相体积分数
    NVar // 变量总数
};

// 变量数组类型定义
using VarArray = amrex::GpuArray<amrex::Real, NVar>;
// 变量矩阵类型定义
using VarMatrix = amrex::GpuArray<amrex::Real, NVar * NVar>;

// 定义两相流标志
#define TWO_PHASE 1