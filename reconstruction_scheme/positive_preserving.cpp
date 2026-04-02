/**
 * @file
 * @brief 两相高阶重构保正限制器的实现。
 * @ingroup reconstruction
 */
#include "positive_preserving.H"
#include "deriv.H"

using namespace amrex;

// AMREX_GPU_HOST_DEVICE void apply_positivive_preserving_flux_limiter(
//     const amrex::Box& bx, // 计算区域
//     const amrex::FArrayBox& consin, // 守恒量数组
//     const amrex::FArrayBox& primin, // 原始量数组
//     const amrex::FArrayBox& olddatain, // 原始量数组
//     const int ncomp,
//     amrex::FArrayBox& rhs, // 右端项数组
//     const amrex::GpuArray<amrex::FArrayBox*, BL_SPACEDIM>& flux, // 通量数组
//     const amrex::GpuArray<amrex::FArrayBox*, BL_SPACEDIM>& ncflux, // 通量数组
//     const amrex::Real dt, // 时间步长
//     const int stage, // Runge-Kutta 阶段
//     const ProbParm& pp, // 问题参数
//     const PositiveQuantityFunctions& pqf, // 正限制函数
//     const amrex::GpuArray<amrex::Real, BL_SPACEDIM>& dx // 网格间距
// )
// {

//     // 该函数已废弃，不再使用，直接在通量计算后进行处理。
//     //  return;
//     //  获取数组视图
//     Array4<Real const> const& consin_array = consin.const_array();
//     Array4<Real const> const& oldatain_array = olddatain.const_array();
//     Array4<Real const> const& primin_array = primin.const_array();
//     Array4<Real> const& rhs_array = rhs.array();

//     // 模块 1: 标记违反正性条件的单元
//     constexpr int nRange = 1;
//     amrex::Box bx_with_ghost = amrex::grow(bx, nRange >= 1 ? nRange : 1); // 增加一层 ghost cells
//     amrex::BaseFab<bool> tags(bx_with_ghost, 1); // 标记数组

//     tags.setVal(true); // 初始化所有值为 true
//     Array4<bool> const& tags_array = tags.array();

//     ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
//         VarArray iCons;
//         // 这里注意，Runge-Kutta法为多个前向欧拉法的凸组合
//         // 因此为了保持正性，需要保持的是每个凸组合的正性，即u^n+k_m，其中k_m为rk中某个阶段计算的增量
//         // 而不是使用计算k_m时的u^m计算保正
//         // 这里的dt也要经过rk系数的对应处理
//         for (int ivar = 0; ivar < ncomp; ivar++) {
//             iCons[ivar] = consin_array(i, j, k, ivar) - dt * rhs_array(i, j, k, ivar);
//             // iCons[ivar] = oldatain_array(i, j, k, ivar) - dt * rhs_array(i, j, k, ivar);
//         }
//         auto pqValues = pqf(iCons, pp);
//         for (int idx = 0; idx < pqValues.size(); idx++) {
//             Real pqValue = pqValues[idx];
//             if (pqValue < pqf.eps[idx] || isnan(pqValue)) {
//                 // 统一的三维循环逻辑（兼容1D/2D/3D）
//                 for (int di = -nRange; di <= nRange; ++di) {
//                     for (int dj = -nRange; dj <= nRange; ++dj) {
//                         for (int dk = -nRange; dk <= nRange; ++dk) {
//                             // 根据维度跳过超出范围的迭代
//                             if (BL_SPACEDIM < 3 && dk != 0)
//                                 continue;
//                             if (BL_SPACEDIM < 2 && dj != 0)
//                                 continue;

//                             // 计算邻居单元的索引
//                             const int ni = i + di;
//                             const int nj = j + dj;
//                             const int nk = k + dk;

//                             // 标记邻居单元（无需边界检查，tags_array 带有 ghost cells）
//                             tags_array(ni, nj, nk) = false;
//                         }
//                     }
//                 }
//                 // tags_array(i, j, k) = false;
//                 break; // 如果检测到问题，直接退出循环
//             }
//         }
//     });

//     // 新增模块: 声明并初始化记录 theta 值的 FArrayBox 数组
//     // 每个方向的 theta 值对应一个面（xlo, xhi, ylo, yhi）
//     amrex::FArrayBox theta_fab[BL_SPACEDIM]; // 每个方向一个 FArrayBox

//     for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
//         // 使用 flux[dir] 的 box 来定义 theta 的 box
//         const amrex::Box& flux_box = flux[dir]->box();
//         theta_fab[dir].resize(flux_box, 1); // 初始化大小，1 表示单分量
//         theta_fab[dir].setVal(1.0); // 初始化所有值为 1.0
//     }

//     // 模块 2: 处理通量并计算 theta 值
//     for (int idim = 0; idim < BL_SPACEDIM; ++idim) {

//         Array4<Real> const& flux_array = flux[idim]->array();
//         Array4<Real> const& ncflux_array = ncflux[idim]->array();

//         Array4<Real> theta_array = theta_fab[idim].array();

//         Dim3 offset = amrex::IntVect::TheDimensionVector(idim).dim3();
//         amrex::Box const& flux_box = flux[idim]->box();
//         amrex::GpuArray<amrex::Real, 3U> norm = get_direction(idim);

//         ParallelFor(flux_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
//             bool is_left_positive = tags_array(i - offset.x, j - offset.y, k - offset.z);
//             bool is_right_positive = tags_array(i, j, k);
//             if (is_left_positive && is_right_positive)
//                 return;

//             // 取出一阶原始变量以计算保正使用的一阶通量
//             // VarArray 为 amrex::GpuArray<amrex::Real, NVar>
//             // NVar在定义方程时已定义
//             VarArray iPrimL, iPrimR;
//             for (int ivar = 0; ivar < ncomp; ivar++) {
//                 iPrimL[ivar] = primin_array(i - offset.x, j - offset.y, k - offset.z, ivar);
//                 iPrimR[ivar] = primin_array(i, j, k, ivar);
//             }

//             // 计算一阶守恒通量
//             VarArray fluxFirstOrder = Riemann_solver_1(iPrimL, iPrimR, norm, pp);

//             // 计算一阶非守恒通量，以alpha du/dx为例
//             // 导数部分通量左右单元相同为u_i+1/2
//             // Coef随左右单元不同，例如为alpha_i和alpha_i+1
//             // 注意,AdvectionTerm给出的是由用户定义的对流项的数组
//             // 需要把它压缩成一个与方程组件数相同的数组，
//             // 如alpha du/dx在一维两相流方程中加在第五项上，没有其他的对流项
//             // 那么这个非守恒通量数组只有第五项不为零
//             AdvectionTerm advTerm;
//             auto ncFluxFirstOrder = advTerm.computeFlux(iPrimL, iPrimR, pp, norm);
//             auto ncIndex = advTerm.getTargetTermIndices();

//             VarArray incFluxFirstOrderL = { 0.0 }, incFluxFirstOrderR = { 0.0 };
//             auto ncCoefL = advTerm.computeCoefficient(iPrimL, pp, norm);
//             auto ncCoefR = advTerm.computeCoefficient(iPrimR, pp, norm);

//             // 这里把不同的alpha和u组合起来组成非守恒的左右通量
//             for (int idx = 0; idx < advTerm.num_terms; idx++) {
//                 incFluxFirstOrderL[ncIndex[idx]] += ncFluxFirstOrder[idx] * ncCoefL[idx];
//                 incFluxFirstOrderR[ncIndex[idx]] += ncFluxFirstOrder[idx] * ncCoefR[idx];
//             }

//             // 取出前面计算的高阶守恒和非守恒通量
//             // 注意非守恒通量的数据结构为每个面上存储了2*组件数的变量，前一半是左通量，后一半是右通量
//             VarArray iFlux,
//                 incFluxL, incFluxR;
//             for (int ivar = 0; ivar < ncomp; ivar++) {
//                 iFlux[ivar] = flux_array(i, j, k, ivar);
//                 incFluxL[ivar] = ncflux_array(i, j, k, ivar);
//                 incFluxR[ivar] = ncflux_array(i, j, k, ncomp + ivar);
//             }

//             Real thetaL = 1.0, thetaR = 1.0;

//             // if (!is_left_positive) {
//             // VarArray iConsL = prim_to_cons(iPrimL, pp);
//             VarArray fluxFirstOrderL, iFluxL, iConsL;
//             for (int ivar = 0; ivar < ncomp; ivar++) {
//                 fluxFirstOrderL[ivar] = fluxFirstOrder[ivar] + incFluxFirstOrderL[ivar];
//                 iFluxL[ivar] = iFlux[ivar] + incFluxL[ivar];
//                 iConsL[ivar] = consin_array(i - offset.x, j - offset.y, k - offset.z, ivar);
//                 // iConsL[ivar] = oldatain_array(i - offset.x, j - offset.y, k - offset.z, ivar);
//             }

//             Real alphaL = computeAlpha(dx, iPrimL, pp)[idim];
//             thetaL = compute_theta_for_side(
//                 idim, ncomp, iPrimL, iConsL, iFluxL, fluxFirstOrderL, 2. * dt / dx[idim] / alphaL, pp, pqf);
//             // volatile Real thetaL_temp = compute_theta_for_side(
//             //     idim, ncomp, iPrimL, iConsL, iFluxL, fluxFirstOrderL, 2. * dt / dx[idim] / alphaL, pp, pqf);
//             // }

//             // if (!is_right_positive) {
//             // 初始化 iConsR
//             // VarArray iConsR = prim_to_cons(iPrimR, pp);
//             VarArray fluxFirstOrderR, iFluxR, iConsR;
//             for (int ivar = 0; ivar < ncomp; ivar++) {
//                 fluxFirstOrderR[ivar] = fluxFirstOrder[ivar] + incFluxFirstOrderR[ivar];
//                 iFluxR[ivar] = iFlux[ivar] + incFluxR[ivar];
//                 iConsR[ivar] = consin_array(i, j, k, ivar);
//                 // iConsR[ivar] = oldatain_array(i, j, k, ivar);
//                 // }

//                 Real alphaR = computeAlpha(dx, iPrimR, pp)[idim];
//                 thetaR = compute_theta_for_side(
//                     idim, ncomp, iPrimR, iConsR, iFluxR, fluxFirstOrderR, -2. * dt / dx[idim] / alphaR, pp, pqf);

//                 // volatile Real thetaR_temp = compute_theta_for_side(
//                 //     idim, ncomp, iPrimR, iConsR, iFluxR, fluxFirstOrderR, -2. * dt / dx[idim] / alphaR, pp, pqf);
//             }

//             // 应用凸组合调整插值后的变量值
//             Real theta = min(thetaL, thetaR);
//             theta_array(i, j, k) = theta;
//             apply_convex_combination(fluxFirstOrder, iFlux, theta, 0, ncomp);
//             apply_convex_combination(incFluxFirstOrderL, incFluxL, theta, 0, ncomp);
//             apply_convex_combination(incFluxFirstOrderR, incFluxR, theta, 0, ncomp);
//             for (int ivar = 0; ivar < ncomp; ivar++) {
//                 flux_array(i, j, k, ivar) = iFlux[ivar];
//                 ncflux_array(i, j, k, ivar) = incFluxL[ivar];
//                 ncflux_array(i, j, k, ivar + ncomp) = incFluxR[ivar];
//             }
//         });
//     }

//     // 模块 3: 更新右端项
//     ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
//         bool is_marked = !tags_array(i, j, k);
//         if (!is_marked) {
//             for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
//                 Dim3 offset = amrex::IntVect::TheDimensionVector(idim).dim3();
//                 if (!tags_array(i - offset.x, j - offset.y, k - offset.z) || !tags_array(i + offset.x, j + offset.y, k + offset.z)) {
//                     is_marked = true;
//                     break;
//                 }
//             }
//         }

//         if (is_marked) {
//             VarArray rhs_temp = { 0.0 };
//             // 计算重新差分后的 rhs 值
//             for (int idim = 0; idim < BL_SPACEDIM; ++idim) {
//                 Array4<Real> const& flux_array = flux[idim]->array(); // 守恒通量
//                 Array4<Real> const& ncflux_array = ncflux[idim]->array(); // 非守恒通量

//                 Dim3 offset = amrex::IntVect::TheDimensionVector(idim).dim3();

//                 for (int ivar = 0; ivar < ncomp; ++ivar) {
//                     Real flux_diff = (flux_array(i + offset.x, j + offset.y, k + offset.z, ivar) - flux_array(i, j, k, ivar)) / dx[idim];
//                     // 注意要取左面的右值和右面的左值
//                     Real ncflux_diff = (ncflux_array(i + offset.x, j + offset.y, k + offset.z, ivar) - ncflux_array(i, j, k, ivar + ncomp)) / dx[idim];
//                     rhs_temp[ivar] += flux_diff + ncflux_diff;
//                 }
//             }

//             // 更新 rhs_array
//             for (int ivar = 0; ivar < ncomp; ++ivar) {
//                 rhs_array(i, j, k, ivar) = rhs_temp[ivar];
//             }
//         }
//     });

//     // 创建一个标志变量数组，用于记录每个线程的错误状态
//     // bool any_error = false;

//     // ParallelFor(bx, [&] AMREX_GPU_DEVICE(int i, int j, int k) {
//     //     VarArray iCons;
//     //     for (int ivar = 0; ivar < ncomp; ivar++) {
//     //         iCons[ivar] = consin_array(i, j, k, ivar) - dt * rhs_array(i, j, k, ivar);
//     //         // iCons[ivar] = oldatain_array(i, j, k, ivar) - dt * rhs_array(i, j, k, ivar);
//     //     }
//     //     auto pqValues = pqf(iCons, pp);

//     //     bool hasNegativeOrNan = false; // 标志变量，用于判断是否需要记录错误
//     //     for (int idx = 0; idx < pqValues.size(); idx++) {
//     //         Real pqValue = pqValues[idx];
//     //         if (pqValue < 0.0 || isnan(pqValue)) {
//     //             hasNegativeOrNan = true; // 检测到负值或 NaN
//     //             break;
//     //         }
//     //     }

//     //     if (hasNegativeOrNan) {
//     //         // 获取当前线程的全局索引
//     //         int global_idx = (i - bx.smallEnd(0)) + (j - bx.smallEnd(1)) * bx.length(0) + (k - bx.smallEnd(2)) * bx.length(0) * bx.length(1);

//     //         any_error = true;

//     //         // 打印当前单元坐标
//     //         printf("Error detected at (%d, %d, %d)\n", i, j, k);

//     //         // 打印距离为1的单元的 tags_array 情况
//     //         printf("Tags array around (%d, %d, %d):\n", i, j, k);

//     //         for (int di = -1; di <= 1; ++di) {
//     //             for (int dj = -1; dj <= 1; ++dj) {
//     //                 for (int dk = -1; dk <= 1; ++dk) {

//     //                     if (BL_SPACEDIM < 3 && dk != 0)
//     //                         continue;
//     //                     if (BL_SPACEDIM < 2 && dj != 0)
//     //                         continue;

//     //                     int ni = i + di;
//     //                     int nj = j + dj;
//     //                     int nk = k + dk;

//     //                         // 打印每个邻居单元的 tag 值
//     //                         printf("%d ", tags_array(ni, nj, nk));
//     //                 }
//     //             }
//     //             // 在不同层之间增加空行以便区分
//     //             printf("\n");
//     //         }

//     //                     // 打印相邻面的 theta_fab 数据
//     //         printf("Theta_fab values around (%d, %d, %d):\n", i, j, k);
//     //         for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
//     //             printf("Direction %d:\n", dir);

//     //             // 定义当前单元和相邻单元的索引偏移
//     //             int offsets[2] = {0, 1}; // 每个方向有两个面：当前面和相邻面

//     //             for (int offset : offsets) {
//     //                 int ni = i, nj = j, nk = k;

//     //                 // 根据方向调整索引
//     //                 if (dir == 0) ni += offset; // x 方向
//     //                 if (dir == 1) nj += offset; // y 方向
//     //                 if (dir == 2 && BL_SPACEDIM >= 3) nk += offset; // z 方向（仅当维度为3时有效）

//     //                 // 确保索引在有效范围内
//     //                 printf("%.6e ", theta_fab[dir].array()(ni, nj, nk));
//     //             }

//     //             // 换行以便区分不同方向
//     //             printf("\n");
//     //         }

//     //         for (int idx = 0; idx < pqValues.size(); idx++) {
//     //             Real pqValue = pqValues[idx];
//     //             std::cout<<pqValue<<' ';
//     //         }

//     //         std::cout << std::endl;
//     //     }
//     // });

//     // 如果检测到任何错误，则调用 amrex::Abort
//     // if (any_error) {
//     //     amrex::Abort("Detected negative or NaN value in positive-preserving check. Stopping execution.");
//     // }
// }

// 计算侧面的theta值
// 输入：维度idim，组件数ncomp，原始变量iPrim，守恒变量iCons，通量iFlux，一阶通量fluxFirstOrder，乘数mult，问题参数pp，保正函数pqf
// 输出：theta值
AMREX_GPU_HOST_DEVICE
amrex::Real
compute_theta_for_side(int idim, int ncomp,
    const VarArray& iPrim, const VarArray& iCons,
    const VarArray& iFlux,
    const VarArray& fluxFirstOrder,
    const amrex::Real& mult,
    const ProbParm& pp,
    const PositiveQuantityFunctions& pqf)
{
    VarArray consStar, consHLLC;
    // 计算原始变量对应的通量
    VarArray flux = prim_to_flux(iPrim, pp, get_direction(idim));
    
    // 计算守恒变量的中间值
    for (int ivar = 0; ivar < ncomp; ivar++) {
        consStar[ivar] = iCons[ivar] - mult * (iFlux[ivar] - flux[ivar]);
        consHLLC[ivar] = iCons[ivar] - mult * (fluxFirstOrder[ivar] - flux[ivar]);
    }
    
    // 计算theta值
    Real theta = compute_theta(consStar, consHLLC, pp, pqf);
    return theta;
}

// 应用保正限制器
// 输入：守恒变量值W，高阶插值后的守恒变量值W_interpolated，保正函数结构体pqf，问题参数pp，组件的数量ncomp
// 输出：无，直接修改W_interpolated
AMREX_GPU_HOST_DEVICE void apply_positive_limiter(
    const VarArray& W, // 守恒变量值
    VarArray& W_interpolated, // 高阶插值后的守恒变量值
    const PositiveQuantityFunctions& pqf, // 保正函数结构体
    const ProbParm& pp, // 问题参数
    int ncomp // 组件的数量
)
{
    // 计算theta值
    amrex::Real theta = compute_theta(W_interpolated, W, pp, pqf);
    // 应用凸组合
    apply_convex_combination(W, W_interpolated, theta, 0, ncomp);
}

// 计算theta值
// 输入：守恒变量consStar，守恒变量consHLLC，问题参数pp，保正函数pqf
// 输出：theta值
AMREX_GPU_HOST_DEVICE
amrex::Real
compute_theta(VarArray const& consStar, VarArray const& consHLLC, ProbParm const& pp, PositiveQuantityFunctions const& pqf)
{
    // 计算theta_0的lambda函数
    auto compute_theta_0 = [] AMREX_GPU_HOST_DEVICE(
                               amrex::Real positiveQuantity_W, // 原始变量中的正值
                               amrex::Real positiveQuantity_W_interpolated, // 插值后变量中的正值
                               amrex::Real EPSILON // 正值的阈值
                               ) -> amrex::Real {
        return (EPSILON - positiveQuantity_W) / (positiveQuantity_W_interpolated - positiveQuantity_W);
    };
    
    // 计算保正变量值
    auto pqStar = pqf(consStar, pp), pqHLLC = pqf(consHLLC, pp);
    amrex::Real theta = 1.;
    
    // 遍历所有保正变量
    for (int ifunc = 0; ifunc < pqStar.size(); ifunc++) {
        // 检查保正变量是否小于阈值
        if (pqStar[ifunc] < pqf.eps[ifunc]) {
            // 0.99的系数可以防止截断误差出问题
            amrex::Real temptheta = compute_theta_0(pqHLLC[ifunc], pqStar[ifunc], pqf.eps[ifunc]);
            theta = amrex::min(theta, temptheta * 0.99);
        }
    }

    return theta;
}
