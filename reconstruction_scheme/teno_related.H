/**
 * @file
 * @brief TENO 平滑度判据与相关辅助函数。
 * @ingroup reconstruction
 */

#pragma once
#include <AMReX_Array.H>
#include <AMReX_Math.H>
#include <algorithm>
#include <cmath>

// Teno5_CongZ_Flag函数
// 输入：五个点的数组q
// 输出：是否需要使用TENO格式的标志
AMREX_GPU_HOST_DEVICE
bool Teno5_CongZ_Flag(const amrex::GpuArray<amrex::Real, 5>& q)
{
    using namespace amrex;
    Real eps = 1e-8; // 小量，避免除以零
    GpuArray<amrex::Real, 3> beta = {
        1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),
        1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),
        1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)
    };

    unsigned short minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin(); // 找到最小beta的索引
    unsigned short maxBeta = std::max_element(beta.begin(), beta.end()) - beta.begin(); // 找到最大beta的索引

    constexpr amrex::Real power = 1.0 / 6.0;
    constexpr amrex::Real CT = std::pow(1.5 * 1e-3, power); // 阈值参数
    // constexpr amrex::Real CT = 0.157042; // O1 test
    constexpr amrex::Real CT_1 = 1 - CT;

    amrex::Real tau = std::abs(beta[2] - beta[0]); // 计算tau
    amrex::Real rr = CT * tau - CT_1 * beta[minBeta];
    amrex::Real ll = tau * beta[minBeta];

    // 判断条件
    bool cond1 = (ll + eps >= rr * beta[0]); // 如果原条件为 false，则 cond1 为 true
    bool cond2 = (ll + eps >= rr * beta[1]); // 如果原条件为 false，则 cond2 为 true
    bool cond3 = (ll + eps >= rr * beta[2]); // 如果原条件为 false，则 cond3 为 true
    bool res = !(cond1 && cond2 && cond3); // 如果任何一个条件为 false，则返回 true

    // 如果任何一个条件为 false，则返回 true
    return res;
}

// template <int SpaceDim>
