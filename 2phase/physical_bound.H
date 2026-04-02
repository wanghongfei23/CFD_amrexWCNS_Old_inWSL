/**
 * @file
 * @brief 两相状态的物理可接受性判据与辅助函数声明。
 * @ingroup two_phase_model
 */

#pragma once
#include "Prob_Parm.H"
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

// 结构体，用于存储固定数量的保正函数
struct PositiveQuantityFunctions {
    amrex::Real eps[5]; // 使用 Vector 存储容差值
    amrex::Real epsHardSwitch[5];

    // 构造函数，初始化容差值
    AMREX_GPU_HOST_DEVICE PositiveQuantityFunctions()
        : eps { 1e-10, 1e-10, 1e-10, 1e-10, 1e-8 }
        // : eps { 0.0, 0.0, 0.0, 0.0, 0.0 }
        , epsHardSwitch { 1e-11, 1e-11, 1e-11, 1e-11, 1e-9 }
    {
    }

    // 成员函数 operator()，计算需要保正的变量
    AMREX_GPU_HOST_DEVICE amrex::GpuArray<amrex::Real, 5> operator()(const VarArray& U, const ProbParm& pp) const;
};

AMREX_GPU_HOST_DEVICE
amrex::GpuArray<amrex::Real, BL_SPACEDIM> computeAlpha(const amrex::GpuArray<amrex::Real,
                                                           BL_SPACEDIM>& dx,
    const VarArray& primitiveVars,
    const ProbParm& pp);