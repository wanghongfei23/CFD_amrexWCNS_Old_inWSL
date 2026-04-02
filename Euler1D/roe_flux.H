/**
 * @file
 * @brief 一维单相 Euler 参考实现中的 Roe 与 HLLC 通量函数。
 * @ingroup legacy_euler
 */

#pragma once
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

/**
 * @brief 计算一维 Roe 数值通量。
 * @param rl 左状态密度。
 * @param rr 右状态密度。
 * @param ul 左状态速度。
 * @param ur 右状态速度。
 * @param pl 左状态压力。
 * @param pr 右状态压力。
 * @param gamma 比热比。
 * @return 三分量守恒通量。
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real, 3> roe_flux_1d(amrex::Real rl, amrex::Real rr, amrex::Real ul, amrex::Real ur, amrex::Real pl, amrex::Real pr, amrex::Real gamma)
{
    // reference: https://blog.csdn.net/Tankrun1997/article/details/132743487
    amrex::GpuArray<amrex::Real, 3> res; // 结果通量
    enum { L, R }; // 左侧和右侧

    // 计算左右侧的总焓
    amrex::GpuArray<amrex::Real, 2> H = { pl / rl * gamma / (gamma - 1) + (ul * ul) / 2, pr / rr * gamma / (gamma - 1) + (ur * ur) / 2 };

    // 计算Roe平均
    amrex::Real rBar = sqrt(rl * rr); // 密度平均
    amrex::Real uBar = (ul * sqrt(rl) + ur * sqrt(rr)) / (sqrt(rl) + sqrt(rr)); // 速度平均
    amrex::Real HBar = (H[L] * sqrt(rl) + H[R] * sqrt(rr)) / (sqrt(rl) + sqrt(rr)); // 焓平均
    amrex::Real cBar = sqrt((gamma - 1) * (HBar - uBar * uBar / 2)); // 声速平均
    amrex::Real cBar2 = (gamma - 1) * (HBar - uBar * uBar / 2); // 声速平方

    // 计算左右侧的差值
    amrex::Real dr = rr - rl; // 密度差
    amrex::Real du = ur - ul; // 速度差
    amrex::Real dp = pr - pl; // 压力差

    // 计算特征向量
    amrex::Real K1[3] = { 1, uBar, uBar * uBar / 2 };
    amrex::Real K2[3] = { 1, uBar - cBar, HBar - uBar * cBar };
    amrex::Real K3[3] = { 1, uBar + cBar, HBar + uBar * cBar };

    // 计算特征值
    amrex::Real lambda[3] = {
        std::abs(uBar), std::abs(uBar + cBar), std::abs(uBar - cBar)
    };

    // 特征值限制
    amrex::Real eps = 0.1 * (std::abs(uBar) + cBar);
    for (int i = 0; i < 3; i++) {
        if (lambda[i] < eps) {
            lambda[i] = (lambda[i] * lambda[i] + eps * eps) / (2.0 * eps);
        }
    }

    // 计算波强度
    amrex::Real alpha[5] = {
        lambda[0] * (dr - dp / cBar2),
        lambda[1] / (2.0 * cBar2) * (dp + rBar * cBar * du),
        lambda[2] / (2.0 * cBar2) * (dp - rBar * cBar * du)
    };
    alpha[3] = alpha[0] + alpha[1] + alpha[2];
    alpha[4] = cBar * (alpha[1] - alpha[2]);

    // 计算左右侧通量
    amrex::Real FL[3] = {
        rl * ul, // 左侧密度通量
        rl * ul * ul + pl, // 左侧动量通量
        rl * H[L] * ul // 左侧能量通量
    };
    amrex::Real FR[3] = {
        rr * ur, // 右侧密度通量
        rr * ur * ur + pr, // 右侧动量通量
        rr * H[R] * ur // 右侧能量通量
    };

    // 计算Roe通量
    for (int i = 0; i < 3; i++) {
        res[i] = (FL[i] + FR[i]) / 2.0;
    }
    res[0] -= 0.5 * alpha[3];
    res[1] -= 0.5 * (alpha[3] * uBar + alpha[4]);
    res[2] -= 0.5 * (HBar * alpha[3] + uBar * alpha[4] - cBar2 * alpha[0] / (gamma - 1));

    return res;
}

/**
 * @brief 计算一维 HLLC 数值通量。
 * @param rl 左状态密度。
 * @param rr 右状态密度。
 * @param ul 左状态速度。
 * @param ur 右状态速度。
 * @param pl 左状态压力。
 * @param pr 右状态压力。
 * @param gamma 比热比。
 * @return 三分量守恒通量。
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real, 3> hllc_flux_1D(amrex::Real rl, amrex::Real rr, amrex::Real ul, amrex::Real ur, amrex::Real pl, amrex::Real pr, amrex::Real gamma)
{
    // Reference: https://zhuanlan.zhihu.com/p/583555029
    amrex::GpuArray<amrex::Real, 3> res; // 结果通量
    enum { L, R }; // 左侧和右侧

    // 计算左右侧的总焓
    amrex::GpuArray<amrex::Real, 2> H = {
        pl / rl * gamma / (gamma - 1) + (ul * ul) / 2,
        pr / rr * gamma / (gamma - 1) + (ur * ur) / 2
    };

    // 计算左右侧的RT值
    amrex::GpuArray<amrex::Real, 2> RT = { pl / rl, pr / rr };

    // 计算左右侧的声速
    amrex::Real cl = std::sqrt(gamma * RT[L]);
    amrex::Real cr = std::sqrt(gamma * RT[R]);

    // 计算Roe平均
    amrex::Real rBar = std::sqrt(rl * rr); // 密度平均
    amrex::Real uBar = (ul * std::sqrt(rl) + ur * std::sqrt(rr)) / (std::sqrt(rl) + std::sqrt(rr)); // 速度平均
    amrex::Real HBar = (H[L] * std::sqrt(rl) + H[R] * std::sqrt(rr)) / (std::sqrt(rl) + std::sqrt(rr)); // 焓平均
    amrex::Real cBar = std::sqrt((gamma - 1) * (HBar - uBar * uBar / 2)); // 声速平均

    // 计算波速
    amrex::Real SL = std::min(uBar - cBar, ul - cl); // 左侧波速
    amrex::Real SR = std::max(uBar + cBar, ur + cr); // 右侧波速

    // 计算中间波速
    amrex::Real Sstar = ((pr - pl) + (rl * ul * (SL - ul) - rr * ur * (SR - ur))) / (rl * (SL - ul) - rr * (SR - ur));

    // 根据波速位置计算通量
    if (SL >= 0) {
        // 全部在左侧
        res[0] = rl * ul; // 密度通量
        res[1] = rl * ul * ul + pl; // 动量通量
        res[2] = rl * H[L] * ul; // 能量通量
    } else if (Sstar >= 0) {
        // 在左侧和中间之间
        amrex::Real pStar = pl + rl * (SL - ul) * (Sstar - ul); // 中间压力
        amrex::GpuArray<amrex::Real, 3> U = { rl, rl * ul, pl / (gamma - 1) + rl * ul * ul / 2 }; // 左侧守恒变量
        amrex::GpuArray<amrex::Real, 3> F = { rl * ul, rl * ul * ul + pl, rl * H[L] * ul }; // 左侧通量
        amrex::GpuArray<amrex::Real, 3> D = { 0, 1, Sstar }; // 中间状态向量

        for (int i = 0; i < 3; i++) {
            res[i] = (Sstar * (SL * U[i] - F[i]) + SL * pStar * D[i]) / (SL - Sstar);
        }
    } else if (SR >= 0) {
        // 在中间和右侧之间
        amrex::Real pStar = pr + rr * (SR - ur) * (Sstar - ur); // 中间压力
        amrex::GpuArray<amrex::Real, 3> U = { rr, rr * ur, pr / (gamma - 1) + rr * ur * ur / 2 }; // 右侧守恒变量
        amrex::GpuArray<amrex::Real, 3> F = { rr * ur, rr * ur * ur + pr, rr * H[R] * ur }; // 右侧通量
        amrex::GpuArray<amrex::Real, 3> D = { 0, 1, Sstar }; // 中间状态向量

        for (int i = 0; i < 3; i++) {
            res[i] = (Sstar * (SR * U[i] - F[i]) + SR * pStar * D[i]) / (SR - Sstar);
        }
    } else {
        // 全部在右侧
        res[0] = rr * ur; // 密度通量
        res[1] = rr * ur * ur + pr; // 动量通量
        res[2] = rr * H[R] * ur; // 能量通量
    }

    return res;
}