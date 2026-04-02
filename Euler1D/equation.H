/**
 * @file
 * @brief 一维单相 Euler 参考实现的状态方程与初值定义。
 * @ingroup legacy_euler
 */

#pragma once

#include "Prob_Parm.H"
#include "roe_flux.H"
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>

// 初始化数据函数
// 输入：级别level，时间time，边界box，状态数组state，网格间距dx，问题下界problo，问题参数pp
// 功能：设置初始条件
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void initdata(int level, amrex::Real time,
    const amrex::Box& box, amrex::FArrayBox& state,
    const amrex::Real* dx, const amrex::Real* problo, const ProbParm& pp)
{
    using namespace amrex;
    constexpr Real gamma = 1.4; // 气体绝热指数
    constexpr Real split_pos = 0.5; // 分割位置
    AMREX_ASSERT(BL_SPACEDIM == 1); // 确保是一维问题
    const int* lo = box.loVect(); // 网格下界
    const int* hi = box.hiVect(); // 网格上界
    Array4<Real> const& state_array = state.array(); // 状态数组

// 遍历所有网格点
#ifdef AMREX_USE_OMP
#pragma omp parallel for if (Gpu::notInLaunchRegion())
#endif
    for (int i = lo[0]; i <= hi[0]; ++i) {
        // 计算当前网格中心的物理坐标
        Real x = problo[0] + (i + 0.5) * dx[0];
        // 判断当前网格点是否在分割点左侧或右侧
        Real rho, u, p;
        // sod shock tube
        //  if (x <= split_pos) {
        //      // 左侧初始条件
        //      rho = 1.0;
        //      u = 0.0;
        //      p = 1.0;
        //  } else {
        //      // 右侧初始条件
        //      rho = 0.125;
        //      u = 0.0;
        //      p = 0.1;
        //  }
        // lax shock tube
        // if (x < split_pos) {
        //     // 左侧初始条件 (0 ≤ x < 0.5)
        //     rho = 0.445;
        //     u = 0.698;
        //     p = 3.528;
        // } else {
        //     // 右侧初始条件 (0.5 ≤ x ≤ 1)
        //     rho = 0.5;
        //     u = 0.0;
        //     p = 0.5710;
        // }

        // linear avection
        if (x >= 0.3 && x <= 0.4) {
            // 区间 [0.3, 0.4] 内的初始条件
            rho = 10.0;
            u = 1.0;
            p = 1.0;
        } else {
            // 其他区间的初始条件
            rho = 1.0;
            u = 1.0;
            p = 1.0;
        }

        // if (x >= 0.3 && x <= 0.4) {
        //     // 区间 [0.3, 0.4] 内的初始条件
        //     rho = 20.0;
        //     u = 1.0;
        //     p = 1.0;
        // } else {
        //     // 其他区间的初始条件
        //     rho = 0.5;
        //     u = 1.0;
        //     p = 1.0;
        // }

        // Shu-Osher
        // constexpr Real split_pos_shuOsher = 1.0;
        // if (x < split_pos_shuOsher) {
        //     // 左侧初始条件 (0 ≤ x < 1)
        //     rho = 3.857;
        //     u = 2.629;
        //     p = 10.333;
        // } else {
        //     // 右侧初始条件 (1 ≤ x ≤ 10)
        //     rho = 1.0 + 0.2 * sin(5.0 * x); // ρ = 1 + 0.2 * sin(5x)
        //     u = 0.0;
        //     p = 1.0;
        // }

        // sedov
        // if (i == 200) {
        //     rho = 1;
        //     u = 0;
        //     p = 3200000.0 / dx[0] * (gamma - 1);
        // } else {
        //     rho = 1;
        //     u = 0;
        //     p = 1e-6 / (gamma - 1);
        // }

        // 双稀疏波
        // if (x < 0) {
        //     rho = 1;
        //     u = -2;
        //     p = 0.4;
        // } else {
        //     rho = 1;
        //     u = 2;
        //     p = 0.4;
        // }

        // 根据欧拉方程计算守恒变量
        state_array(i, 0, 0, 0) = rho; // 密度
        state_array(i, 0, 0, 1) = rho * u; // 动量
        state_array(i, 0, 0, 2) = p / (gamma - 1) + 0.5 * rho * u * u; // 总能量
    }
}

// 守恒变量转换为原始变量
// 输入：守恒变量cons，问题参数pp
// 输出：原始变量
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
VarArray cons_to_prim(VarArray cons, ProbParm const& pp)
{
    amrex::Real r, ru, rE;
    r = cons[0]; // 密度
    ru = cons[1]; // 动量
    rE = cons[2]; // 总能量
    amrex::Real u = ru / r; // 速度
    amrex::Real E = rE / r; // 单位质量总能量
    amrex::Real e = -u * u / 2 + E; // 单位质量内能

    amrex::Real RT = (pp.gamma - 1) * e; // 气体常数乘温度
    amrex::Real p = r * RT; // 压力
    return { r, u, p };
}

// 为Array4格式的守恒变量转换为原始变量
// 输入：边界ghost_box，原始变量数组prim_array，状态数组statein_array，问题参数pp
// 功能：将Array4格式的守恒变量转换为原始变量
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void cons_to_prim_for_array4(const amrex::Box& ghost_box,
    amrex::Array4<amrex::Real> const& prim_array,
    amrex::Array4<amrex::Real const> const& statein_array,
    ProbParm const& pp = ProbParm())
{
    ParallelFor(ghost_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        int ncomp = prim_array.nComp();
        VarArray input;
        for (int ivar = 0; ivar < ncomp; ivar++) {
            input[ivar] = statein_array(i, j, k, ivar);
        }

        VarArray primi = cons_to_prim(input, pp);

        for (int ivar = 0; ivar < ncomp; ivar++) {
            prim_array(i, j, k, ivar) = primi[ivar];
        }
    });
}

// 原始变量转换为守恒变量
// 输入：原始变量prim，问题参数pp
// 输出：守恒变量
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
VarArray prim_to_cons(VarArray prim, ProbParm const& pp)
{
    amrex::Real gamma = pp.gamma; // 气体绝热指数
    amrex::Real rt = prim[0], ut = prim[1], pt = prim[2]; // 密度、速度、压力
    amrex::Real ekt = (ut * ut) / 2; // 单位质量动能
    amrex::Real rut = rt * ut; // 动量
    amrex::Real ret = pt / (gamma - 1) + rt * ekt; // 总能量
    return { rt, rut, ret };
}

// Riemann求解器
// 输入：网格坐标i,j,k，原始变量数组prim_array，左侧状态stateL，右侧状态stateR，法向量norm，问题参数pp
// 输出：通量
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
VarArray
Riemann_solver_1(int /*i*/, int /*j*/, int /*k*/, const amrex::Array4<const amrex::Real>& /*prim_array*/, VarArray const& stateL, VarArray const& stateR, amrex::GpuArray<amrex::Real, 3> norm, ProbParm const& pp)
{
    return hllc_flux_1D(stateL[0], stateR[0], stateL[1], stateR[1], stateL[2], stateR[2], pp.gamma);
}

// 原始变量转换为通量
// 输入：原始变量prim，问题参数pp，法向量norm
// 输出：通量
 template <int idim = 0>
AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        VarArray
        prim_to_flux(VarArray prim, ProbParm pp, amrex::GpuArray<amrex::Real, 3> norm = { 1, 0, 0 })
{
    amrex::Real r = prim[0], u = prim[1], p = prim[2]; // 密度、速度、压力
    amrex::Real gamma = pp.gamma; // 气体绝热指数

    // 计算通量
    VarArray flux; // 通量：1. 密度通量，2. 动量通量，3. 能量通量
    flux[0] = r * u;
    flux[1] = r * u * u + p;
    flux[2] = (gamma / (gamma - 1) * p + r * (u * u) / 2) * u;
    // std::cout << flux[0] << ' ' << flux[1] << ' ' << flux[2] << '\n';

    return flux;
}

// 计算最大特征速度
// 输入：网格坐标i,j,k，状态数组state，最大速度数组umax，问题参数pp，网格间距dx
// 功能：计算每个网格点的最大特征速度
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void cal_u_max(
    int i, int j, int k,
    amrex::Array4<amrex::Real const> const& state,
    amrex::Array4<amrex::Real> const& umax,
    ProbParm const& pp, // 可选参数，支持不同气体模型
    amrex::GpuArray<amrex::Real, BL_SPACEDIM> dx) noexcept
{
    amrex::Real gamma = pp.gamma; // 气体绝热指数
    amrex::Real rho = state(i, j, k, 0); // 密度
    // if (rho <= 1.e-10) {
    //     umax(i, j, k) = 0.0;
    //     return;
    // }
    amrex::Real rhou = state(i, j, k, 1); // 动量
    amrex::Real rhoE = state(i, j, k, 2); // 总能量

    amrex::Real u = rhou / rho; // 速度
    amrex::Real E = rhoE / rho; // 单位质量总能量
    amrex::Real e = E - 0.5 * u * u; // 单位质量内能

    amrex::Real RT = (gamma - 1) * e; // 气体常数乘温度
    amrex::Real p = rho * RT; // 压力

    umax(i, j, k) = (sqrt(gamma * RT) + std::abs(u)) / dx[0]; // 最大特征速度
}