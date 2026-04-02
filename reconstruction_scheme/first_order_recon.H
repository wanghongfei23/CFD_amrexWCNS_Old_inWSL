/**
 * @file
 * @brief 一阶 Godunov 型重构实现。
 * @ingroup reconstruction
 */

#pragma once
#include "equation.H"
#include "linear_algebra.H"

// 一阶重构函数
// 输入：网格坐标i,j,k，原始变量数组primin_array，变量数量NUM_STATE，问题参数pp，法向量norm，网格间距dx
// 输出：左右侧状态值
// 模板参数dir：重构方向

template <int dir = 0>
AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        amrex::GpuArray<VarArray, 2>
        reconstruction(int i, int j, int k,
            amrex::Array4<amrex::Real const> const& primin_array,
            int NUM_STATE, ProbParm const& pp,
            amrex::GpuArray<amrex::Real, 3> norm,
            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
{
    using namespace amrex;
    // 初始化左右侧状态
    amrex::GpuArray<VarArray, 2> res;
    Dim3 offset = amrex::IntVect::TheDimensionVector(dir).dim3(); // 方向偏移

    // 对于一阶重构，左右侧状态就是单元格中心的值
    for (int ivar = 0; ivar < NUM_STATE; ++ivar) {
        res[0][ivar] = primin_array(i - offset.x, j - offset.y, k - offset.z, ivar); // 左侧状态
        res[1][ivar] = primin_array(i, j, k, ivar); // 右侧状态
    }

    return res;
}