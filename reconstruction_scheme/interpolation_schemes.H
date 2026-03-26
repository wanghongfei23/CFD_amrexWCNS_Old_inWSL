// 插值方案头文件
// 包含各种数值插值方法的实现，如WENO、TENO等

#pragma once
#include <AMReX_Array.H>
#include <AMReX_Math.H>
// AMREX_GPU_HOST_DEVICE
// AMREX_FORCE_INLINE
// amrex::Real weno5_Z(amrex::GpuArray<amrex::Real, 5> q)
// {
//     amrex::Real eps = 1e-40;
//     amrex::GpuArray<amrex::Real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
//     amrex::GpuArray<amrex::Real, 3> beta;
//     beta[0] = 1.0 / 1.0 * std::pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * std::pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);
//     beta[1] = 1.0 / 1.0 * std::pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * std::pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);
//     beta[2] = 1.0 / 1.0 * std::pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * std::pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

//     amrex::GpuArray<amrex::Real, 3> u;
//     u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
//     u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
//     u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

//     // amrex::GpuArray<amrex::Real(*)(amrex::Real,amrex::Real,amrex::Real),3> u={&u1,&u2,&u3};

//     amrex::Real sumbeta = 0, result = 0;
//     amrex::Real C = 1, qq = 2, tau = std::abs(beta[2] - beta[0]);
//     for (int i = 0; i < 3; i++) {
//         beta[i] = gamma[i] * (C + std::pow(tau / (beta[i] + eps), qq));
//         sumbeta += beta[i];
//     }
//     for (int i = 0; i < 3; i++)
//         result += beta[i] * u[i];
//     return result / sumbeta;
// }

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real weno5_Z(amrex::GpuArray<amrex::Real, 5> q)
{
    // 定义 sqr lambda 表达式，并增加浮点溢出检查
    auto sqr = [&](amrex::Real x) -> amrex::Real {
        return x * x;
    };

    amrex::Real eps = 1e-40;
    amrex::GpuArray<amrex::Real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
    amrex::GpuArray<amrex::Real, 3> beta;
    beta[0] = 1.0 / 1.0 * sqr(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2]) + 1.0 / 4.0 * sqr(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2]);
    beta[1] = 1.0 / 1.0 * sqr(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3]) + 1.0 / 4.0 * sqr(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3]);
    beta[2] = 1.0 / 1.0 * sqr(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4]) + 1.0 / 4.0 * sqr(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4]);

    amrex::GpuArray<amrex::Real, 3> u;
    u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

    amrex::Real sumbeta = 0, result = 0;
    amrex::Real C = 1, qq = 2, tau = std::abs(beta[2] - beta[0]);
    for (int i = 0; i < 3; i++) {
        beta[i] = gamma[i] * (C + sqr(tau / (beta[i] + eps)));
        sumbeta += beta[i];
    }
    for (int i = 0; i < 3; i++)
        result += beta[i] * u[i];
    return result / sumbeta;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real weno5_Z_norm(amrex::GpuArray<amrex::Real, 5> q)
{
    // 定义 sqr lambda 表达式，并增加浮点溢出检查
    auto sqr = [&](amrex::Real x) -> amrex::Real {
        return x * x;
    };

    amrex::Real eps = 1e-40;
    amrex::GpuArray<amrex::Real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
    amrex::GpuArray<amrex::Real, 3> beta;

    // 对 q[0], q[1], q[3], q[4] 进行代数操作，除以 q[2]
    amrex::Real q0_norm = q[0] / q[2];
    amrex::Real q1_norm = q[1] / q[2];
    amrex::Real q3_norm = q[3] / q[2];
    amrex::Real q4_norm = q[4] / q[2];

    // 计算 beta 时使用归一化的值
    beta[0] = 1.0 / 1.0 * sqr(1.0 * q0_norm - 2.0 * q1_norm + 1.0) + 1.0 / 4.0 * sqr(1.0 * q0_norm - 4.0 * q1_norm + 3.0);
    beta[1] = 1.0 / 1.0 * sqr(1.0 * q1_norm - 2.0 + 1.0 * q3_norm) + 1.0 / 4.0 * sqr(1.0 * q1_norm + 0.0 - 1.0 * q3_norm);
    beta[2] = 1.0 / 1.0 * sqr(1.0 - 2.0 * q3_norm + 1.0 * q4_norm) + 1.0 / 4.0 * sqr(3.0 - 4.0 * q3_norm + 1.0 * q4_norm);

    amrex::GpuArray<amrex::Real, 3> u;
    u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

    amrex::Real sumbeta = 0, result = 0;
    amrex::Real C = 1, qq = 2, tau = std::abs(beta[2] - beta[0]);
    for (int i = 0; i < 3; i++) {
        beta[i] = gamma[i] * (C + sqr(tau / (beta[i] + eps)));
        sumbeta += beta[i];
    }
    for (int i = 0; i < 3; i++)
        result += beta[i] * u[i];
    return result / sumbeta;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real weno5_JSchen(amrex::GpuArray<amrex::Real, 5> q)
{
    amrex::Real eps = 1e-6;
    amrex::GpuArray<amrex::Real, 3> gamma = { 1.0 / 16.0, 5.0 / 8.0, 5.0 / 16.0 };
    amrex::GpuArray<amrex::Real, 3> beta;
    beta[0] = 1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);
    beta[1] = 1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2);
    beta[2] = 1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2);

    amrex::GpuArray<amrex::Real, 3> u;
    u[0] = 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    u[1] = -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
    u[2] = 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];

    amrex::Real sumbeta = 0, result = 0;
    for (int i = 0; i < 3; i++) {
        beta[i] = gamma[i] / pow(eps + beta[i], 2.0);
        sumbeta += beta[i];
    }
    for (int i = 0; i < 3; i++)
        result += beta[i] * u[i];
    return result / sumbeta;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real weno_is(amrex::GpuArray<amrex::Real, 5> q)
{
    // Constants
    constexpr amrex::Real epsilon = 1.0e-15;
    constexpr amrex::Real d0 = 15.0 / 32.0;
    constexpr amrex::Real d1 = 5.0 / 32.0;
    constexpr amrex::Real d2 = 5.0 / 16.0;
    constexpr amrex::Real d3 = 1.0 / 16.0;

    // Interpolated values from sub-stencils
    amrex::Real u0 = 0.5 * (q[2] + q[3]); // S0
    amrex::Real u1 = 0.5 * (-q[1] + 3.0 * q[2]); // S1
    amrex::Real u2 = (3.0 * q[2] + 6.0 * q[3] - q[4]) / 8.0; // S2
    amrex::Real u3 = (3.0 * q[0] - 10.0 * q[1] + 15.0 * q[2]) / 8.0; // S3

    // Smoothness indicators
    amrex::Real beta0 = std::pow(q[2] - q[3], 2);
    amrex::Real beta1 = std::pow(q[1] - q[2], 2);
    amrex::Real beta01 = (13.0 / 12.0) * std::pow(q[1] - 2.0 * q[2] + q[3], 2) + (1.0 / 4.0) * std::pow(q[1] - q[3], 2);
    amrex::Real beta2 = (13.0 / 12.0) * std::pow(q[2] - 2.0 * q[3] + q[4], 2) + (1.0 / 4.0) * std::pow(3.0 * q[2] - 4.0 * q[3] + q[4], 2);
    amrex::Real beta3 = (13.0 / 12.0) * std::pow(q[0] - 2.0 * q[1] + q[2], 2) + (1.0 / 4.0) * std::pow(3.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2);

    // Reference smoothness indicator
    amrex::Real tau5 = (13.0 / 12.0) * std::pow(q[4] - 4.0 * q[3] + 6.0 * q[2] - 4.0 * q[1] + q[0], 2) + (1.0 / 4.0) * std::pow(q[4] - 2.0 * q[3] + 2.0 * q[1] - q[0], 2);

    // Nonlinear weights
    amrex::Real eta0 = d0 * (1.0 + tau5 / (beta0 + epsilon) + tau5 / (beta01 + epsilon));
    amrex::Real eta1 = d1 * (1.0 + tau5 / (beta1 + epsilon) + tau5 / (beta01 + epsilon));
    amrex::Real eta2 = d2 * (1.0 + tau5 / (beta2 + epsilon));
    amrex::Real eta3 = d3 * (1.0 + tau5 / (beta3 + epsilon));

    // Normalize weights
    amrex::Real sum_eta = eta0 + eta1 + eta2 + eta3;
    amrex::Real omega0 = eta0 / sum_eta;
    amrex::Real omega1 = eta1 / sum_eta;
    amrex::Real omega2 = eta2 / sum_eta;
    amrex::Real omega3 = eta3 / sum_eta;

    // WENO interpolation
    amrex::Real u_weno = omega0 * u0 + omega1 * u1 + omega2 * u2 + omega3 * u3;

    return u_weno;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real WENO3(amrex::GpuArray<amrex::Real, 3> q)
{
    // q[0] = u_{i-1}
    // q[1] = u_i
    // q[2] = u_{i+1}

    const amrex::Real u_im1 = q[0];
    const amrex::Real u_i = q[1];
    const amrex::Real u_ip1 = q[2];

    // Candidate interpolations for u_{i+1/2}
    const amrex::Real u0 = 0.5 * u_im1 + 1.5 * u_i; // left-biased
    const amrex::Real u1 = -0.5 * u_i + 1.5 * u_ip1; // right-biased

    // Smoothness indicators (simplified for 3-point stencil)
    const amrex::Real beta0 = (u_i - u_im1) * (u_i - u_im1);
    const amrex::Real beta1 = (u_ip1 - u_i) * (u_ip1 - u_i);

    // Ideal weights for third-order accuracy
    constexpr amrex::Real d0 = 1.0 / 3.0;
    constexpr amrex::Real d1 = 2.0 / 3.0;

    // Small epsilon to avoid division by zero
    constexpr amrex::Real eps = 1e-6;

    // Nonlinear weights
    const amrex::Real alpha0 = d0 / ((eps + beta0) * (eps + beta0));
    const amrex::Real alpha1 = d1 / ((eps + beta1) * (eps + beta1));

    const amrex::Real w0 = alpha0 / (alpha0 + alpha1);
    const amrex::Real w1 = alpha1 / (alpha0 + alpha1);

    // Final interpolated value at x_{i+1/2}
    return w0 * u0 + w1 * u1;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real Teno5_CongZ(amrex::GpuArray<amrex::Real, 5> q)
{
    amrex::Real eps = 1e-10;
    amrex::GpuArray<amrex::Real, 3> beta = {
        1.0 / 1.0 * pow(1.0 * q[0] - 2.0 * q[1] + 1.0 * q[2], 2) + 1.0 / 4.0 * pow(1.0 * q[0] - 4.0 * q[1] + 3.0 * q[2], 2),
        1.0 / 1.0 * pow(1.0 * q[1] - 2.0 * q[2] + 1.0 * q[3], 2) + 1.0 / 4.0 * pow(1.0 * q[1] + 0.0 * q[2] - 1.0 * q[3], 2),
        1.0 / 1.0 * pow(1.0 * q[2] - 2.0 * q[3] + 1.0 * q[4], 2) + 1.0 / 4.0 * pow(3.0 * q[2] - 4.0 * q[3] + 1.0 * q[4], 2)
    };

    unsigned short minBeta = std::min_element(beta.begin(), beta.end()) - beta.begin();
    // constexpr amrex::Real CT = std::pow(1.5 * 1e-5, 1.0 / 6.0);
    constexpr amrex::Real CT = 0.2;
    constexpr amrex::Real CT_1 = 1 - CT;
    amrex::Real tau = std::abs(beta[2] - beta[0]);
    amrex::Real rr = CT * tau - CT_1 * beta[minBeta];
    amrex::Real ll = tau * beta[minBeta] + eps;
    unsigned short flag = 0;

    if (ll < rr * beta[0])
        flag += 1;
    if (ll < rr * beta[1])
        flag += 2;
    if (ll < rr * beta[2])
        flag += 4;

    switch (flag) {
    case 0:
        /* 1,1,1 */
        return 3.0 / 128.0 * q[0] - 5.0 / 32.0 * q[1] + 45.0 / 64.0 * q[2] + 15.0 / 32.0 * q[3] - 5.0 / 128.0 * q[4];
    case 1:
        /* 0,1,1 */
        return -1.0 / 16.0 * q[1] + 9.0 / 16.0 * q[2] + 9.0 / 16.0 * q[3] - 1.0 / 16.0 * q[4];
    case 2:
        /* 1,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    case 3:
        /* 0,0,1 */
        return 3.0 / 8.0 * q[2] + 3.0 / 4.0 * q[3] - 1.0 / 8.0 * q[4];
    case 4:
        /* 1,1,0 */
        return 1.0 / 16.0 * q[0] - 5.0 / 16.0 * q[1] + 15.0 / 16.0 * q[2] + 5.0 / 16.0 * q[3];
    case 5:
        /* 0,1,0 */
        return -1.0 / 8.0 * q[1] + 3.0 / 4.0 * q[2] + 3.0 / 8.0 * q[3];
        // return WENO3({ q[1], q[2], q[3] });
    case 6:
        /* 1,0,0 */
        return 3.0 / 8.0 * q[0] - 5.0 / 4.0 * q[1] + 15.0 / 8.0 * q[2];
    default:
        /* 0,0,0 */
        return q[2];
    }
}