/**
 * @file
 * @brief MUSCL 限幅重构实现。
 * @ingroup reconstruction
 */

#pragma once
#include "linear_algebra.H"
#include <AMReX_Array4.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <vector>

// MUSCL插值函数
// 输入：三个相邻点的值q1, q2, q3
// 输出：插值后的结果
AMREX_FORCE_INLINE
amrex::Real musclInterpolation(amrex::Real q1, amrex::Real q2, amrex::Real q3)
{
    using namespace amrex;
    amrex::Real delta; // 插值增量
    amrex::Real deltam, deltap; // 左右差分
    deltam = q2 - q1; // 左侧差分
    deltap = q3 - q2; // 右侧差分

    // minmod限制器
    amrex::Real beta = 1.0; // 限制器参数
    if (deltap > 0) {
        delta = max(0.0, max(min(beta * deltam, deltap), min(deltam, beta * deltap)));
    } else {
        delta = min(0.0, min(max(beta * deltam, deltap), max(deltam, beta * deltap)));
    }
    return q2 + delta * 0.5; // 返回插值结果
}

// 重构函数
// 输入：网格坐标i,j,k，原始变量数组primin_array，变量数量NUM_STATE，问题参数pp，法向量norm
// 输出：左右侧重构后的变量值
// 模板参数dir：重构方向

template <int dir = 0>
AMREX_GPU_HOST_DEVICE
    // AMREX_FORCE_INLINE
    amrex::GpuArray<VarArray, 2>
    reconstruction(int i, int j, int k,
        amrex::Array4<amrex::Real const> const& primin_array,
        int NUM_STATE, ProbParm const& pp, amrex::GpuArray<amrex::Real, 3> norm = { 1, 0, 0 })
{
    amrex::Vector<amrex::Real> res(NUM_STATE * 2); // 结果向量，存储左右侧重构值
    // auto iter = res.begin();
    for (int ivar = 0; ivar < NUM_STATE; ivar++) { // 遍历所有变量
        amrex::GpuArray<amrex::Real, 4> state_temp; // 临时存储四个点的值
        // 沿指定方向遍历，获取四个点的值
        traverse_in_direction_with_ii<dir>(i, j, k, -2, 1, [&](int li, int lj, int lk, int ii) {
            state_temp[ii] = primin_array(li, lj, lk, ivar);
        });
        // 计算左侧重构值
        res[ivar] = musclInterpolation(state_temp[0], state_temp[1], state_temp[2]);
        // 计算右侧重构值
        res[ivar + NUM_STATE] = musclInterpolation(state_temp[3], state_temp[2], state_temp[1]);
    }

    return res;
}
