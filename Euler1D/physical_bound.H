// 物理边界条件头文件
// 包含保正函数和边界条件计算相关的声明
#pragma once
#include "Prob_Parm.H"
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

// 结构体，用于存储固定数量的保正函数
// 保正函数确保计算过程中物理量（如密度、能量）保持正值
struct PositiveQuantityFunctions {
    amrex::Real eps[2]; // 容差值数组，用于保正计算
    amrex::Real epsHardSwitch[2]; // 硬开关容差值，用于更严格的保正检查

    // 构造函数，初始化容差值
    AMREX_GPU_HOST_DEVICE PositiveQuantityFunctions()
    {
        eps[0] = 1e-10; // 密度的容差
        eps[1] = 1e-8; // rho c^2的容差
        epsHardSwitch[0] = 1e-11; // 密度的硬开关容差
        epsHardSwitch[1] = 1e-9; // rho c^2的硬开关容差
    }

    // 成员函数 operator()，计算需要保正的变量
    // 输入：守恒变量U和问题参数pp
    // 输出：保正变量数组，包含密度和rho c^2
    AMREX_GPU_HOST_DEVICE amrex::GpuArray<amrex::Real, 2> operator()(const VarArray& U, const ProbParm& pp) const
    {
        amrex::GpuArray<amrex::Real, 2> result; // 输出两个保正变量

        // 第一个保正变量：密度
        result[0] = U[0]; // 假设 U[0] 是密度

        // 第二个保正变量：密度乘声速的平方
        amrex::Real rho = U[0]; // 密度
        amrex::Real u = U[1] / rho; // x 方向速度
        amrex::Real E = U[2]; // 总能量
        amrex::Real p = (pp.gamma - 1.0) * (E - 0.5 * rho * u * u); // 压力公式
        result[1] = pp.gamma * p; // 密度乘声速的平方

        return result;
    }
};

// 计算Alpha值的函数，用于边界条件处理
// 输入：网格步长dx，原始变量primitiveVars，问题参数pp
// 输出：Alpha值数组，用于边界条件计算
AMREX_GPU_HOST_DEVICE
amrex::GpuArray<amrex::Real, BL_SPACEDIM> computeAlpha(const amrex::GpuArray<amrex::Real,
                                                           BL_SPACEDIM>& dx,
    const VarArray& primitiveVars,
    const ProbParm& pp);