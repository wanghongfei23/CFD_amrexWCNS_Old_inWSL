/**
 * @file
 * @brief 二维单相 Euler 参考实现中的保正判据与边界辅助声明。
 * @ingroup legacy_euler
 */

#pragma once
#include "Prob_Parm.H"
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

/**
 * @brief 二维 Euler 状态的保正目标函数集合。
 */
struct PositiveQuantityFunctions {
    amrex::Real eps[2]; // 使用 Vector 存储容差值
    amrex::Real epsHardSwitch[2];

    // 构造函数，初始化容差值
    AMREX_GPU_HOST_DEVICE PositiveQuantityFunctions()
    {
        eps[0] = 1e-10; // 密度的容差
        eps[1] = 1e-8; // rho c^2的容差
        epsHardSwitch[0] = 1e-11; // 密度的硬开关容差
        epsHardSwitch[1] = 1e-9; // rho c^2的硬开关容差
    }

    // 成员函数 operator()，计算需要保正的变量
    AMREX_GPU_HOST_DEVICE amrex::GpuArray<amrex::Real, 2> operator()(const VarArray& U, const ProbParm& pp) const
    {
        amrex::GpuArray<amrex::Real, 2> result; // 输出两个保正变量

        // 第一个保正变量：密度
        result[0] = U[0]; // 假设 U[0] 是密度

        // 第二个保正变量：压力
        amrex::Real rho = U[0]; // 密度
        amrex::Real u = U[1] / rho; // x 方向速度
        amrex::Real v = U[2] / rho; // y 方向速度
        amrex::Real E = U[3]; // 总能量
        amrex::Real p = (pp.gamma - 1.0) * (E - 0.5 * rho * (u * u + v * v)); // 压力公式
        result[1] = pp.gamma * p;

        return result;
    }
};

AMREX_GPU_HOST_DEVICE
amrex::GpuArray<amrex::Real, BL_SPACEDIM> computeAlpha(const amrex::GpuArray<amrex::Real,
                                                           BL_SPACEDIM>& dx,
    const VarArray& primitiveVars,
    const ProbParm& pp);