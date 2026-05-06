/**
 * @file
 * @brief 五阶 THINC/WENO 混合重构与界面捕捉实现。
 * @ingroup reconstruction
 */

#pragma once
#include "eigen_system.H"
#include "interpolation_schemes.H"
#include "linear_algebra.H"
#include "positive_preserving.H"
#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
// #include <vector>

// phi_THINC函数
// 输入：三个相邻点的phi值phi_im1, phi_i, phi_ip1，位置x，参数beta
// 输出：THINC重构后的phi值
AMREX_GPU_HOST_DEVICE
amrex::Real phi_THINC(amrex::Real phi_im1, amrex::Real phi_i, amrex::Real phi_ip1, amrex::Real x, amrex::Real beta = 1.8)
{
    auto sgn = [](amrex::Real val) -> amrex::Real {
        return val > 0 ? 1.0 : (val < 0 ? -1.0 : 0.0);
    };

    amrex::Real min_phi = std::min(phi_im1, phi_ip1); // 最小值
    amrex::Real max_phi = std::max(phi_im1, phi_ip1); // 最大值
    amrex::Real phi_bar = (max_phi + min_phi) / 2; // 平均值
    amrex::Real phi_delta = (max_phi - min_phi) / 2; // 半差值
    amrex::Real theta = sgn(phi_ip1 - phi_im1); // 符号

    if ((phi_im1 - phi_i) * (phi_i - phi_ip1) > 0) { //(phi_i > min_phi && phi_i < max_phi) {
        amrex::Real di = atanh(1.0 - 1e-8 - (2.0 - 2e-8) * phi_i) / beta; // 计算界面位置
        // amrex::Real di = -(1.0 / (2.0 * beta)) * std::log((std::exp(2.0 * beta * phi_i) - 1.0) / (std::exp(2.0 * beta) - std::exp(2.0 * beta * phi_i)));
        return 0.5 + (0.5 - 1e-8) * theta * tanh(beta * (x - theta * di)); // THINC重构

        // amrex::Real di = atanh((phi_bar - phi_i) / phi_delta) / beta;
        // return phi_bar + phi_delta * theta * tanh(beta * (x - theta * di));
    } else {
        return phi_i; // 不满足条件时直接返回中心值
    }
}

AMREX_GPU_HOST_DEVICE
amrex::Real phi_THINC3(amrex::Real phi_im1, amrex::Real phi_i, amrex::Real phi_ip1, amrex::Real x, amrex::Real beta = 1.8)
{
    auto sgn = [](amrex::Real val) -> amrex::Real {
        return val > 0 ? 1.0 : (val < 0 ? -1.0 : 0.0);
    };

    amrex::Real min_phi = std::min(phi_im1, phi_ip1);
    amrex::Real max_phi = std::max(phi_im1, phi_ip1);
    amrex::Real phi_bar = (max_phi + min_phi) / 2;
    amrex::Real phi_delta = (max_phi - min_phi) / 2;
    amrex::Real theta = sgn(phi_ip1 - phi_im1);

    if ((phi_im1 - phi_i) * (phi_i - phi_ip1) > 0) { // 单调性条件
        // 使用新的重构公式
        amrex::Real kappa = beta; // 将 beta 视为 kappa 参数
        amrex::Real tanh_k4 = std::tanh(kappa / 4.0);
        amrex::Real tanh_k2 = std::tanh(kappa / 2.0);

        amrex::Real numerator = 2.0 * tanh_k4 * (phi_ip1 - phi_i) * (phi_i - phi_im1);
        amrex::Real denominator = (phi_ip1 - phi_i) * (tanh_k2 - tanh_k4) + (phi_i - phi_im1) * (tanh_k2 + tanh_k4);

        amrex::Real delta_u = numerator / denominator;

        return phi_i + delta_u;
    } else {
        return phi_i; // 如果不满足单调性条件，直接返回 phi_i
    }
}

AMREX_GPU_HOST_DEVICE
amrex::Real phi_THINC_bounded(amrex::Real phi_im1, amrex::Real phi_i, amrex::Real phi_ip1, amrex::Real x, amrex::Real beta = 1.8)
{
    auto sgn = [](amrex::Real val) -> amrex::Real {
        return val > 0 ? 1.0 : (val < 0 ? -1.0 : 0.0);
    };

    amrex::Real min_phi = std::min(phi_im1, phi_ip1);
    amrex::Real max_phi = std::max(phi_im1, phi_ip1);
    amrex::Real phi_bar = (max_phi + min_phi) / 2;
    amrex::Real phi_delta = (max_phi - min_phi) / 2;
    amrex::Real theta = sgn(phi_ip1 - phi_im1);

    if ((phi_im1 - phi_i) * (phi_i - phi_ip1) > 0) { //(phi_i > min_phi && phi_i < max_phi) {
        // amrex::Real di = atanh(1.0 - 2.0 * phi_i) / beta;
        // // amrex::Real di = -(1.0 / (2.0 * beta)) * std::log((std::exp(2.0 * beta * phi_i) - 1.0) / (std::exp(2.0 * beta) - std::exp(2.0 * beta * phi_i)));
        // return 0.5 + 0.5 * theta * tanh(beta * (x - theta * di));
        amrex::Real di = atanh((phi_bar - phi_i) / phi_delta) / beta;
        return phi_bar + phi_delta * theta * tanh(beta * (x - theta * di));
    } else {
        return phi_i;
    }
}

// 原始变量重构函数
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::GpuArray<VarArray, 2>
primitive_reconstruction(int i, int j, int k,
    amrex::Array4<amrex::Real const> const& primin_array,
    int NUM_STATE, ProbParm const& pp,
    amrex::GpuArray<amrex::Real, 3> norm,
    amrex::Dim3 offset,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
{
    // 初始化 qL 和 qR 用于存储左/右状态的物理量
    amrex::GpuArray<amrex::GpuArray<amrex::Real, 5>, NVar> qL, qR;

    // 遍历指定方向，计算每个位置的物理量
    for (int ii = -3; ii <= 2; ii++) {
        VarArray primj;
        for (int ivar = 0; ivar < NUM_STATE; ivar++) {
            primj[ivar] = primin_array(i + ii * offset.x, j + ii * offset.y, k + ii * offset.z, ivar);
        }

        int left_index = ii + 3;
        if (left_index < 5) {
            for (int ivar = 0; ivar < NUM_STATE; ivar++) {
                qL[ivar][left_index] = primj[ivar];
            }
        }

        int right_index = 5 - left_index;
        if (right_index < 5) {
            for (int ivar = 0; ivar < NUM_STATE; ivar++) {
                qR[ivar][right_index] = primj[ivar];
            }
        }
    }

    // 对左、右状态的数据进行 WENO5 重构
    // 【临时定位】thinc+格式选择
    VarArray primReconL, primReconR;
    for (int ivar = 0; ivar < NUM_STATE; ivar++) {
        primReconL[ivar] = TCNS_AS(qL[ivar]);
        primReconR[ivar] = TCNS_AS(qR[ivar]);
    }
    // TCNS_AS
    // weno_is

    VarArray primCellL, primCellR;

    for (int ivar = 0; ivar < NUM_STATE; ivar++) {
        primCellL[ivar] = primin_array(i - 1 * offset.x, j - 1 * offset.y, k - 1 * offset.z, ivar);
        primCellR[ivar] = primin_array(i - 0 * offset.x, j - 0 * offset.y, k - 0 * offset.z, ivar);
    }

    // 将原始变量转换为守恒变量
    VarArray consCellL = prim_to_cons(primCellL, pp);
    VarArray consCellR = prim_to_cons(primCellR, pp);

    VarArray consReconL = prim_to_cons(primReconL, pp);
    VarArray consReconR = prim_to_cons(primReconR, pp);

    // 应用保正限制器
    const PositiveQuantityFunctions pqf;
    apply_positive_limiter(consCellL, consReconL, pqf, pp, NUM_STATE);
    apply_positive_limiter(consCellR, consReconR, pqf, pp, NUM_STATE);

    // 将保正后的守恒变量转换回原始变量
    VarArray primL = cons_to_prim(consReconL, pp);
    VarArray primR = cons_to_prim(consReconR, pp);

    // 返回拼接后的结果
    return { primL, primR };
}

// 特征变量重构函数
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::GpuArray<VarArray, 2>
characteristic_reconstruction(int i, int j, int k,
    amrex::Array4<amrex::Real const> const& primin_array,
    int NUM_STATE, ProbParm const& pp,
    amrex::GpuArray<amrex::Real, 3> norm,
    amrex::Dim3 offset,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
{
    // 初始化 primCellL 和 primCellR，用于存储从 primin_array 中提取的物理量
    VarArray primCellL, primCellR;

    for (int ivar = 0; ivar < NUM_STATE; ivar++) {
        primCellL[ivar] = primin_array(i - 1 * offset.x, j - 1 * offset.y, k - 1 * offset.z, ivar);
        primCellR[ivar] = primin_array(i - 0 * offset.x, j - 0 * offset.y, k - 0 * offset.z, ivar);
    }

    // 使用 Roe 平均法对 primLR 进行处理，得到 primRef
    EigenSystem eig(primCellL, primCellR, norm, pp);

    // 初始化变量 it，qL 和 qR 用于存储左/右特征矩阵的变换结果
    amrex::GpuArray<amrex::GpuArray<amrex::Real, 5>, NVar> qL, qR;

    // 遍历指定方向，计算每个位置的物理量，并进行特征矩阵变换
    for (int ii = -3; ii <= 2; ii++) {
        VarArray primj;
        for (int ivar = 0; ivar < NUM_STATE; ivar++) {
            primj[ivar] = primin_array(i + ii * offset.x, j + ii * offset.y, k + ii * offset.z, ivar);
        }

        int left_index = ii + 3;
        auto charTemp = eig.prim_to_char(primj, pp);

        if (left_index < 5) {
            for (int ivar = 0; ivar < NUM_STATE; ivar++) {
                qL[ivar][left_index] = charTemp[ivar];
            }
        }

        int right_index = 5 - left_index;
        if (right_index < 5) {
            for (int ivar = 0; ivar < NUM_STATE; ivar++) {
                qR[ivar][right_index] = charTemp[ivar];
            }
        }
    }

    // 对左、右特征矩阵的数据进行 WENO5 重构
    VarArray charReconL, charReconR;
    for (int ivar = 0; ivar < NUM_STATE; ivar++) {
        charReconL[ivar] = weno_is(qL[ivar]);
        charReconR[ivar] = weno_is(qR[ivar]);
    }
    // charReconL[0] = phi_THINC(qL[0][1], qL[0][2], qL[0][3], 1.0, 1.5);
    // charReconR[0] = phi_THINC(qR[0][1], qR[0][2], qR[0][3], 1.0, 1.5);
    // 使用 LA_M_X 进行矩阵-向量乘法，将计算结果存入 res
    VarArray consL = eig.char_to_cons(charReconL, pp);
    VarArray consR = eig.char_to_cons(charReconR, pp);

    // 将原始变量 primCellL 和 primCellR 转换为守恒变量
    VarArray consCellL = prim_to_cons(primCellL, pp);
    VarArray consCellR = prim_to_cons(primCellR, pp);

    const PositiveQuantityFunctions pqf;

    // 对 consL 和 consR 应用保正限制器
    apply_positive_limiter(consCellL, consL, pqf, pp, NUM_STATE);
    apply_positive_limiter(consCellR, consR, pqf, pp, NUM_STATE);

    // 将保正后的守恒变量转换回原始变量
    VarArray primL = cons_to_prim(consL, pp);
    VarArray primR = cons_to_prim(consR, pp);

    // 返回拼接后的结果
    return { primL, primR };
}
constexpr amrex::Real switchd_thredhold = 1e-3;
// // 自定义重构函数，其中体积分数和相密度采用其他方法重构
// AMREX_GPU_HOST_DEVICE
// AMREX_FORCE_INLINE
// amrex::GpuArray<VarArray, 2>
// mixture_reconstruction(int i, int j, int k,
//     amrex::Array4<amrex::Real const> const& primin_array,
//     int NUM_STATE, ProbParm const& pp,
//     amrex::GpuArray<amrex::Real, 3> norm,
//     amrex::Dim3 offset,
//     amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
// {
//     // 初始化 qL 和 qR 用于存储左/右状态的物理量
//     amrex::GpuArray<amrex::GpuArray<amrex::Real, 5>, NVar> qL, qR;

//     VarArray primCellL, primCellR;
//     for (int ivar = 0; ivar < NUM_STATE; ivar++) {
//         primCellL[ivar] = primin_array(i - 1 * offset.x, j - 1 * offset.y, k - 1 * offset.z, ivar);
//         primCellR[ivar] = primin_array(i - 0 * offset.x, j - 0 * offset.y, k - 0 * offset.z, ivar);
//     }

//     // 获取左右单元的体积分数
//     amrex::Real alpha1_L = primCellL[Alpha1];
//     amrex::Real alpha1_R = primCellR[Alpha1];

//     // 计算对应的另一相体积分数
//     amrex::Real alpha2_L = 1.0 - alpha1_L;
//     amrex::Real alpha2_R = 1.0 - alpha1_R;

//     // 检测是否需要使用原始变量重构
//     // bool use_primitive_reconstruction = ((alpha1_L > 1e-3 && alpha2_L > 1e-3) || (alpha1_R > 1e-3 && alpha2_R > 1e-3));

//     // 遍历指定方向，计算每个位置的物理量
//     for (int ii = -3; ii <= 2; ii++) {
//         VarArray primj;
//         for (int ivar = 0; ivar < NUM_STATE; ivar++) {
//             primj[ivar] = primin_array(i + ii * offset.x, j + ii * offset.y, k + ii * offset.z, ivar);
//         }

//         int left_index = ii + 3;
//         if (left_index < 5) {
//             for (int ivar = 0; ivar < NUM_STATE; ivar++) {
//                 qL[ivar][left_index] = primj[ivar];
//             }
//         }

//         int right_index = 5 - left_index;
//         if (right_index < 5) {
//             for (int ivar = 0; ivar < NUM_STATE; ivar++) {
//                 qR[ivar][right_index] = primj[ivar];
//             }
//         }
//     }

//     // 对左、右状态的数据进行 WENO5 重构，但跳过体积分数和相密度
//     VarArray primReconL, primReconR;
//     constexpr amrex::Real beta = 1.2;

//     amrex::GpuArray<amrex::Real, 3> interface_normal_L = compute_normal_vector(i - offset.x, j - offset.y, k - offset.z, primin_array, dx);
//     amrex::Real normL_dir = interface_normal_L[0] * norm[0] + interface_normal_L[1] * norm[1] + interface_normal_L[2] * norm[2];
//     amrex::Real betaL = std::abs(normL_dir) * beta + 0.01;
//     primReconL[Alpha1] = phi_THINC(qL[Alpha1][1], qL[Alpha1][2], qL[Alpha1][3], 1.0, betaL);
//     // 2. 重构相密度 AlphaRho1 和 AlphaRho2，依赖于体积分数的重构结果
//     if ((alpha1_L > switchd_thredhold && alpha2_L > switchd_thredhold)) { // && ((qL[Alpha1][1] - qL[Alpha1][2]) * (qL[Alpha1][2] - qL[Alpha1][3]) > 0)) {
//         // 1. 重构体积分数 Alpha1
//         primReconL[AlphaRho1] = primCellL[AlphaRho1] + (primCellL[AlphaRho1] / primCellL[Alpha1]) * (primReconL[Alpha1] - primCellL[Alpha1]);
//         primReconL[AlphaRho2] = primCellL[AlphaRho2] + (primCellL[AlphaRho2] / (1.0 - primCellL[Alpha1])) * (primCellL[Alpha1] - primReconL[Alpha1]);

//         // 3. 对其他变量使用 WENO5 重构
//         for (int ivar = 0; ivar < NUM_STATE; ivar++) {
//             if (ivar == Alpha1 || ivar == AlphaRho1 || ivar == AlphaRho2) {
//                 // 跳过已经处理的变量
//                 continue;
//             }
//             primReconL[ivar] = weno_is(qL[ivar]);
//         }
//     } else {
//         for (int ivar = 0; ivar < NUM_STATE; ivar++) {
//             if (ivar == Alpha1) {
//                 // 跳过已经处理的变量
//                 continue;
//             }
//             primReconL[ivar] = weno_is(qL[ivar]);
//         }
//         // primReconL[AlphaRho1] = primReconL[AlphaRho1] < 0 ? 0 : primReconL[AlphaRho1];
//         // primReconL[AlphaRho2] = primReconL[AlphaRho2] < 0 ? 0 : primReconL[AlphaRho2];
//         // primReconL[Alpha1] = primReconL[Alpha1] < 0 ? 0 : primReconL[Alpha1];
//     }

//     // 1. 重构体积分数 Alpha1
//     amrex::GpuArray<amrex::Real, 3> normal_vector_R = compute_normal_vector(i, j, k, primin_array, dx);
//     amrex::Real normR_dir = normal_vector_R[0] * norm[0] + normal_vector_R[1] * norm[1] + normal_vector_R[2] * norm[2];
//     amrex::Real betaR = std::abs(normR_dir) * beta + 0.01;
//     primReconR[Alpha1] = phi_THINC(qR[Alpha1][1], qR[Alpha1][2], qR[Alpha1][3], 1.0, betaR);
//     if ((alpha1_R > switchd_thredhold && alpha2_R > switchd_thredhold)) { // && ((qR[Alpha1][1] - qR[Alpha1][2]) * (qR[Alpha1][2] - qR[Alpha1][3]) > 0)) {
//         // 2. 重构相密度 AlphaRho1 和 AlphaRho2，依赖于体积分数的重构结果
//         primReconR[AlphaRho1] = primCellR[AlphaRho1] + (primCellR[AlphaRho1] / primCellR[Alpha1]) * (primReconR[Alpha1] - primCellR[Alpha1]);
//         primReconR[AlphaRho2] = primCellR[AlphaRho2] + (primCellR[AlphaRho2] / (1.0 - primCellR[Alpha1])) * (primCellR[Alpha1] - primReconR[Alpha1]);
//         // 3. 对其他变量使用 WENO5 重构
//         for (int ivar = 0; ivar < NUM_STATE; ivar++) {
//             if (ivar == Alpha1 || ivar == AlphaRho1 || ivar == AlphaRho2) {
//                 // 跳过已经处理的变量
//                 continue;
//             }
//             primReconR[ivar] = weno_is(qR[ivar]);
//         }
//     } else {
//         for (int ivar = 0; ivar < NUM_STATE; ivar++) {
//             if (ivar == Alpha1) {
//                 // 跳过已经处理的变量
//                 continue;
//             }
//             primReconR[ivar] = weno_is(qR[ivar]);
//         }
//         // primReconR[AlphaRho1] = primReconR[AlphaRho1] < 0 ? 0 : primReconR[AlphaRho1];
//         // primReconR[AlphaRho2] = primReconR[AlphaRho2] < 0 ? 0 : primReconR[AlphaRho2];
//         // primReconR[Alpha1] = primReconR[Alpha1] < 0 ? 0 : primReconR[Alpha1];
//     }

//     // 将原始变量转换为守恒变量
//     VarArray consCellL = prim_to_cons(primCellL, pp);
//     VarArray consCellR = prim_to_cons(primCellR, pp);

//     VarArray consReconL = prim_to_cons(primReconL, pp);
//     VarArray consReconR = prim_to_cons(primReconR, pp);

//     // 应用保正限制器
//     const PositiveQuantityFunctions pqf;
//     apply_positive_limiter(consCellL, consReconL, pqf, pp, NUM_STATE);
//     apply_positive_limiter(consCellR, consReconR, pqf, pp, NUM_STATE);

//     // 将保正后的守恒变量转换回原始变量
//     VarArray primL = cons_to_prim(consReconL, pp);
//     VarArray primR = cons_to_prim(consReconR, pp);

//     // 返回拼接后的结果
//     return { primL, primR };
// }

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::GpuArray<VarArray, 2>
mixture_reconstruction(int i, int j, int k,
    amrex::Array4<amrex::Real const> const& primin_array,
    int NUM_STATE, ProbParm const& pp,
    amrex::GpuArray<amrex::Real, 3> norm,
    amrex::Dim3 offset,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
{
    // #ifndef AlphaRho1
    // #error "Compilation terminated due to missing AlphaRho1 definition"
    // #endif

    // ... rest of code ...
    // 初始化 qL 和 qR 用于存储左/右状态的物理量
    amrex::GpuArray<amrex::GpuArray<amrex::Real, 5>, NVar> qL, qR;

    VarArray primCellL, primCellR;
    for (int ivar = 0; ivar < NUM_STATE; ivar++) {
        primCellL[ivar] = primin_array(i - 1 * offset.x, j - 1 * offset.y, k - 1 * offset.z, ivar);
        primCellR[ivar] = primin_array(i - 0 * offset.x, j - 0 * offset.y, k - 0 * offset.z, ivar);
    }

    // 获取左右单元的体积分数
    amrex::Real alpha1_L = primCellL[Alpha1];
    amrex::Real alpha1_R = primCellR[Alpha1];

    // 计算对应的另一相体积分数
    amrex::Real alpha2_L = 1.0 - alpha1_L;
    amrex::Real alpha2_R = 1.0 - alpha1_R;

    // 检测是否需要使用原始变量重构
    // bool use_primitive_reconstruction = ((alpha1_L > 1e-3 && alpha2_L > 1e-3) || (alpha1_R > 1e-3 && alpha2_R > 1e-3));

    // 遍历指定方向，计算每个位置的物理量
    for (int ii = -3; ii <= 2; ii++) {
        VarArray primj;
        for (int ivar = 0; ivar < NUM_STATE; ivar++) {
            primj[ivar] = primin_array(i + ii * offset.x, j + ii * offset.y, k + ii * offset.z, ivar);
        }

        int left_index = ii + 3;
        if (left_index < 5) {
            for (int ivar = 0; ivar < NUM_STATE; ivar++) {
                qL[ivar][left_index] = primj[ivar];
            }
        }

        int right_index = 5 - left_index;
        if (right_index < 5) {
            for (int ivar = 0; ivar < NUM_STATE; ivar++) {
                qR[ivar][right_index] = primj[ivar];
            }
        }
    }

    // 对左、右状态的数据进行 WENO5 重构，但跳过体积分数和相密度
    VarArray primReconL, primReconR;
    constexpr amrex::Real beta = 1.0;
    // constexpr amrex::Real betaRho = 1.5;

    // if ((alpha1_L > switchd_thredhold && alpha2_L > switchd_thredhold)) { // && ((qL[Alpha1][1] - qL[Alpha1][2]) * (qL[Alpha1][2] - qL[Alpha1][3]) > 0)) {
    //     // if (true) {
    //     // 1. 重构体积分数 Alpha1
    //     amrex::GpuArray<amrex::Real, 3> interface_normal_L = compute_normal_vector(i - offset.x, j - offset.y, k - offset.z, primin_array, dx);
    //     amrex::Real normL_dir = interface_normal_L[0] * norm[0] + interface_normal_L[1] * norm[1] + interface_normal_L[2] * norm[2];
    //     amrex::Real betaL = std::abs(normL_dir) * beta + 0.01;
    //     primReconL[Alpha1] = phi_THINC3(qL[Alpha1][1], qL[Alpha1][2], qL[Alpha1][3], 1.0, betaL);
    //     // 2. 重构相密度 AlphaRho1 和 AlphaRho2，依赖于体积分数的重构结果
    //     primReconL[AlphaRho1] = primCellL[AlphaRho1] + (primCellL[AlphaRho1] / primCellL[Alpha1]) * (primReconL[Alpha1] - primCellL[Alpha1]);
    //     primReconL[AlphaRho2] = primCellL[AlphaRho2] + (primCellL[AlphaRho2] / (1.0 - primCellL[Alpha1])) * (primCellL[Alpha1] - primReconL[Alpha1]);

    //     // 3. 对其他变量使用 WENO5 重构
    //     for (int ivar = 0; ivar < NUM_STATE; ivar++) {
    //         if (ivar == Alpha1 || ivar == AlphaRho1 || ivar == AlphaRho2) { //|| ivar == AlphaRho2
    //             // 跳过已经处理的变量
    //             continue;
    //         }
    //         primReconL[ivar] = qL[ivar][2]; // weno_is(qL[ivar])
    //         // primReconL[ivar] = weno_is(qL[ivar]);
    //     }
    // } else {
    //     for (int ivar = 0; ivar < NUM_STATE; ivar++) {
    //         primReconL[ivar] = weno_is(qL[ivar]);
    //     }
    // }
    if ((alpha1_L > switchd_thredhold && alpha2_L > switchd_thredhold) && ((qL[Alpha1][1] - qL[Alpha1][2]) * (qL[Alpha1][2] - qL[Alpha1][3]) > 1e-15)) {
        // 1. 计算混合权重函数 f(x) = 4x(1-x)
        amrex::Real avg_alpha1 = 0.5 * (qL[Alpha1][2] + qL[Alpha1][3]); // 计算 qL[Alpha1][2] 和 qL[Alpha1][3] 的平均值
        // amrex::Real mixing_weight = 2.0 * sqrt(avg_alpha1 * (1.0 - avg_alpha1)); // 带入 4x(1-x)
        // amrex::Real mixing_weight = 4.0 * avg_alpha1 * (1.0 - avg_alpha1); // 带入 4x(1-x)
        amrex::Real mixing_weight = 1.0; // 带入 4x(1-x)

        // 2. 重构体积分数 Alpha1
        amrex::GpuArray<amrex::Real, 3> interface_normal_L = compute_normal_vector(i - offset.x, j - offset.y, k - offset.z, primin_array, dx);
        amrex::Real normL_dir = interface_normal_L[0] * norm[0] + interface_normal_L[1] * norm[1] + interface_normal_L[2] * norm[2];
        amrex::Real betaL = std::abs(normL_dir) * beta + 0.01;

        // 当前 THINC 重构结果
        amrex::Real thinc_recon_alpha1 = phi_THINC_bounded(qL[Alpha1][1], qL[Alpha1][2], qL[Alpha1][3], 1.0, betaL);
        // amrex::Real thinc_recon_alpha1 = qL[Alpha1][3];

        // WENO5 重构结果
        amrex::Real weno_recon_alpha1
            = weno_is(qL[Alpha1]);

        // 凸组合：mixing_weight 控制 THINC 和 WENO5 的权重
        primReconL[Alpha1] = mixing_weight * thinc_recon_alpha1 + (1.0 - mixing_weight) * weno_recon_alpha1;

        // 3. 分项重构 AlphaRho1 和 AlphaRho2
        // 3.1 THINC 重构结果
        // amrex::Real thinc_recon_AlphaRho1 = primCellL[AlphaRho1] + (primCellL[AlphaRho1] / primCellL[Alpha1]) * (thinc_recon_alpha1 - primCellL[Alpha1]);
        // amrex::Real thinc_recon_AlphaRho2 = primCellL[AlphaRho2] + (primCellL[AlphaRho2] / (1.0 - primCellL[Alpha1])) * (primCellL[Alpha1] - thinc_recon_alpha1);
        amrex::Real thinc_recon_AlphaRho1 = phi_THINC_bounded(qL[AlphaRho1][1], qL[AlphaRho1][2], qL[AlphaRho1][3], 1.0, betaL);
        amrex::Real thinc_recon_AlphaRho2 = phi_THINC_bounded(qL[AlphaRho2][1], qL[AlphaRho2][2], qL[AlphaRho2][3], 1.0, betaL);

        // 3.2 WENO5 重构结果
        amrex::Real weno_recon_AlphaRho1 = weno_is(qL[AlphaRho1]);
        amrex::Real weno_recon_AlphaRho2 = weno_is(qL[AlphaRho2]);

        // 3.3 凸组合
        primReconL[AlphaRho1] = mixing_weight * thinc_recon_AlphaRho1 + (1.0 - mixing_weight) * weno_recon_AlphaRho1;
        primReconL[AlphaRho2] = mixing_weight * thinc_recon_AlphaRho2 + (1.0 - mixing_weight) * weno_recon_AlphaRho2;

        // 4. 对其他变量使用凸组合
        for (int ivar = 0; ivar < NUM_STATE; ivar++) {
            if (ivar == Alpha1 || ivar == AlphaRho1 || ivar == AlphaRho2) {
                // 跳过已经处理的变量
                continue;
            }

            // 当前 THINC 重构结果（假设为 qL[ivar][2]）
            amrex::Real thinc_recon_var = qL[ivar][2];
            // amrex::Real thinc_recon_var = phi_THINC3(qL[ivar][1], qL[ivar][2], qL[ivar][3], 1.0, betaL);

            // WENO5 重构结果
            amrex::Real weno_recon_var = weno_is(qL[ivar]);

            // 凸组合：mixing_weight 控制 THINC 和 WENO5 的权重
            primReconL[ivar] = mixing_weight * thinc_recon_var + (1.0 - mixing_weight) * weno_recon_var;
        }
    } else {
        for (int ivar = 0; ivar < NUM_STATE; ivar++) {
            primReconL[ivar] = weno_is(qL[ivar]);
        }
    }

    // if ((alpha1_R > switchd_thredhold && alpha2_R > switchd_thredhold)) { // && ((qR[Alpha1][1] - qR[Alpha1][2]) * (qR[Alpha1][2] - qR[Alpha1][3]) > 0)) {
    //     // if (true) {
    //     // 1. 重构体积分数 Alpha1
    //     amrex::GpuArray<amrex::Real, 3> normal_vector_R = compute_normal_vector(i, j, k, primin_array, dx);
    //     amrex::Real normR_dir = normal_vector_R[0] * norm[0] + normal_vector_R[1] * norm[1] + normal_vector_R[2] * norm[2];
    //     amrex::Real betaR = std::abs(normR_dir) * beta + 0.01;
    //     primReconR[Alpha1] = phi_THINC3(qR[Alpha1][1], qR[Alpha1][2], qR[Alpha1][3], 1.0, betaR);
    //     // 2. 重构相密度 AlphaRho1 和 AlphaRho2，依赖于体积分数的重构结果
    //     primReconR[AlphaRho1] = primCellR[AlphaRho1] + (primCellR[AlphaRho1] / primCellR[Alpha1]) * (primReconR[Alpha1] - primCellR[Alpha1]);
    //     primReconR[AlphaRho2] = primCellR[AlphaRho2] + (primCellR[AlphaRho2] / (1.0 - primCellR[Alpha1])) * (primCellR[Alpha1] - primReconR[Alpha1]);
    //     // 3. 对其他变量使用 WENO5 重构
    //     for (int ivar = 0; ivar < NUM_STATE; ivar++) {
    //         if (ivar == Alpha1 || ivar == AlphaRho1 || ivar == AlphaRho2) { //|| ivar == AlphaRho2
    //             // 跳过已经处理的变量
    //             continue;
    //         }
    //         primReconR[ivar] = qR[ivar][2]; // weno_is(qR[ivar]);
    //         // primReconL[ivar] = weno_is(qL[ivar]);
    //     }
    // } else {
    //     for (int ivar = 0; ivar < NUM_STATE; ivar++) {
    //         primReconR[ivar] = weno_is(qR[ivar]);
    //     }
    // }

    if ((alpha1_R > switchd_thredhold && alpha2_R > switchd_thredhold) && ((qR[Alpha1][1] - qR[Alpha1][2]) * (qR[Alpha1][2] - qR[Alpha1][3]) > 1e-15)) {
        // 1. 计算混合权重函数 f(x) = 4x(1-x)
        amrex::Real avg_alpha1 = 0.5 * (qR[Alpha1][2] + qR[Alpha1][3]); // 计算 qR[Alpha1][2] 和 qR[Alpha1][3] 的平均值
        // amrex::Real mixing_weight = 2.0 * sqrt(avg_alpha1 * (1.0 - avg_alpha1)); // 带入 4x(1-x)
        // amrex::Real mixing_weight = 4.0 * avg_alpha1 * (1.0 - avg_alpha1); // 带入 4x(1-x)
        amrex::Real mixing_weight = 1.0; // 带入 4x(1-x)

        // 2. 重构体积分数 Alpha1
        amrex::GpuArray<amrex::Real, 3> normal_vector_R = compute_normal_vector(i, j, k, primin_array, dx);
        amrex::Real normR_dir = normal_vector_R[0] * norm[0] + normal_vector_R[1] * norm[1] + normal_vector_R[2] * norm[2];
        amrex::Real betaR = std::abs(normR_dir) * beta + 0.01;

        // 当前 THINC 重构结果
        amrex::Real thinc_recon_alpha1 = phi_THINC_bounded(qR[Alpha1][1], qR[Alpha1][2], qR[Alpha1][3], 1.0, betaR);
        // amrex::Real thinc_recon_alpha1 = qR[Alpha1][3];

        // WENO5 重构结果
        amrex::Real weno_recon_alpha1 = weno_is(qR[Alpha1]);

        // 凸组合：mixing_weight 控制 THINC 和 WENO5 的权重
        primReconR[Alpha1] = mixing_weight * thinc_recon_alpha1 + (1.0 - mixing_weight) * weno_recon_alpha1;

        // 3. 分项重构 AlphaRho1 和 AlphaRho2
        // 3.1 THINC 重构结果
        // amrex::Real thinc_recon_AlphaRho1 = primCellR[AlphaRho1] + (primCellR[AlphaRho1] / primCellR[Alpha1]) * (thinc_recon_alpha1 - primCellR[Alpha1]);
        // amrex::Real thinc_recon_AlphaRho2 = primCellR[AlphaRho2] + (primCellR[AlphaRho2] / (1.0 - primCellR[Alpha1])) * (primCellR[Alpha1] - thinc_recon_alpha1);
        amrex::Real thinc_recon_AlphaRho1 = phi_THINC_bounded(qR[AlphaRho1][1], qR[AlphaRho1][2], qR[AlphaRho1][3], 1.0, betaR);
        amrex::Real thinc_recon_AlphaRho2 = phi_THINC_bounded(qR[AlphaRho2][1], qR[AlphaRho2][2], qR[AlphaRho2][3], 1.0, betaR);

        // 3.2 WENO5 重构结果
        amrex::Real weno_recon_AlphaRho1 = weno_is(qR[AlphaRho1]);
        amrex::Real weno_recon_AlphaRho2 = weno_is(qR[AlphaRho2]);

        // 3.3 凸组合
        primReconR[AlphaRho1] = mixing_weight * thinc_recon_AlphaRho1 + (1.0 - mixing_weight) * weno_recon_AlphaRho1;
        primReconR[AlphaRho2] = mixing_weight * thinc_recon_AlphaRho2 + (1.0 - mixing_weight) * weno_recon_AlphaRho2;

        // 4. 对其他变量使用凸组合
        for (int ivar = 0; ivar < NUM_STATE; ivar++) {
            if (ivar == Alpha1 || ivar == AlphaRho1 || ivar == AlphaRho2) {
                // 跳过已经处理的变量
                continue;
            }

            // 当前 THINC 重构结果（假设为 qR[ivar][2]）
            amrex::Real thinc_recon_var = qR[ivar][2];
            // amrex::Real thinc_recon_var = phi_THINC3(qR[ivar][1], qR[ivar][2], qR[ivar][3], 1.0, betaR);

            // WENO5 重构结果
            amrex::Real weno_recon_var = weno_is(qR[ivar]);

            // 凸组合：mixing_weight 控制 THINC 和 WENO5 的权重
            primReconR[ivar] = mixing_weight * thinc_recon_var + (1.0 - mixing_weight) * weno_recon_var;
        }
    } else {
        for (int ivar = 0; ivar < NUM_STATE; ivar++) {
            primReconR[ivar] = weno_is(qR[ivar]);
        }
    }

    // 将原始变量转换为守恒变量
    VarArray consCellL = prim_to_cons(primCellL, pp);
    VarArray consCellR = prim_to_cons(primCellR, pp);

    VarArray consReconL = prim_to_cons(primReconL, pp);
    VarArray consReconR = prim_to_cons(primReconR, pp);

    // 应用保正限制器
    const PositiveQuantityFunctions pqf;
    apply_positive_limiter(consCellL, consReconL, pqf, pp, NUM_STATE);
    apply_positive_limiter(consCellR, consReconR, pqf, pp, NUM_STATE);

    // 将保正后的守恒变量转换回原始变量
    VarArray primL = cons_to_prim(consReconL, pp);
    VarArray primR = cons_to_prim(consReconR, pp);

    // 返回拼接后的结果
    return { primL, primR };
}

// 模板函数 reconstruction：进行数值重构，直接对原始变量进行 WENO5 重构
// 参数：
//   i, j, k: 网格的索引位置
//   primin_array: 存储物理量（例如，压力、速度等）的 4D 数组
//   NUM_STATE: 物理状态量的数量（例如，速度、温度等）
//   norm: 用于计算方向的归一化向量，默认为 {1, 0, 0}，表示沿 i 方向的单位向量
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
    // 定义 offset 来控制某一维度的操作
    amrex::Dim3 offset = amrex::IntVect::TheDimensionVector(dir).dim3();

    // 初始化 primCellL 和 primCellR，用于存储从 primin_array 中提取的物理量
    VarArray primCellL, primCellR;

    for (int ivar = 0; ivar < NUM_STATE; ivar++) {
        primCellL[ivar] = primin_array(i - 1 * offset.x, j - 1 * offset.y, k - 1 * offset.z, ivar);
        primCellR[ivar] = primin_array(i - 0 * offset.x, j - 0 * offset.y, k - 0 * offset.z, ivar);
    }

    // 获取左右单元的体积分数
    amrex::Real alpha1_L = primCellL[Alpha1];
    amrex::Real alpha1_R = primCellR[Alpha1];

    // 计算对应的另一相体积分数
    amrex::Real alpha2_L = 1.0 - alpha1_L;
    amrex::Real alpha2_R = 1.0 - alpha1_R;

    // 检测是否需要使用原始变量重构
    bool use_primitive_reconstruction = ((alpha1_L > switchd_thredhold && alpha2_L > switchd_thredhold) && (alpha1_R > switchd_thredhold && alpha2_R > switchd_thredhold));
    // bool use_primitive_reconstruction = std::abs(alpha1_L - alpha1_R) > 1e-3;
    // 根据条件选择重构方法
    if (use_primitive_reconstruction) {
        // 使用原始变量重构
        // return primitive_reconstruction(i, j, k, primin_array, NUM_STATE, pp, norm, offset, dx);
        return mixture_reconstruction(i, j, k, primin_array, NUM_STATE, pp, norm, offset, dx);
    } else {
        // 使用特征变量重构
        return characteristic_reconstruction(i, j, k, primin_array, NUM_STATE, pp, norm, offset, dx);
    }
}
