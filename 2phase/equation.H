/**
 * 两相五方程模型的状态方程、变量变换、Riemann 求解与初值设置
 * 
 * 本文件实现了两相流数值模拟中的核心功能，包括：
 * - 状态方程相关物理量的计算
 * - 守恒变量与原始变量之间的变换
 * - Riemann 求解器的实现
 * - 不同维度下的算例初始化
 */

#pragma once

#include "Prob_Parm.H"
// #include "initialization.H"
// #include "linear_algebra.H"

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
// #include <AMReX_SPACE.H>

/**
 * 指定从状态方程中提取的目标物理量。
 */
enum EOSType {
    Pressure = 0,       // 压力
    SoundSpeed,         // 声速
    SoundSpeedSqrDensity, // 声速平方乘以密度
    TotalEnergy         // 总能量
};
/**
 * 从守恒变量中恢复状态方程相关物理量
 * @tparam eType 目标物理量类型
 * @param cons 守恒变量
 * @param pp 两相物性参数
 * @return 根据 eType 返回压力、声速或辅助量
 */
template <EOSType eType>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
cons_to_eos(VarArray cons, ProbParm const& pp)
{
    using namespace amrex;

    Real rho1 = cons[AlphaRho1];
    Real rho2 = cons[AlphaRho2];
    Real AMREX_D_DECL(rhou = cons[XMom], rhov = cons[YMom], rhow = cons[ZMom]);
    Real rhoE = cons[Energy];
    Real alpha1 = cons[Alpha1];

    Real halfRhoU_sqr = (AMREX_D_TERM(rhou * rhou, +rhov * rhov, +rhow * rhow)) / 2;

    auto computeState = [&](Real gamma, Real pInf, Real rho) -> Real {
        Real internalEnergy = rhoE - halfRhoU_sqr / rho;
        Real pressure = (gamma - 1) * internalEnergy - gamma * pInf;
        Real soundSpeedSqrDensity = gamma * (gamma - 1) * (internalEnergy - pInf);
        Real soundSpeed = sqrt(soundSpeedSqrDensity / rho);

        if constexpr (eType == Pressure) {
            return pressure;
        } else if constexpr (eType == SoundSpeed) {
            return soundSpeed;
        } else if constexpr (eType == SoundSpeedSqrDensity) {
            return soundSpeedSqrDensity;
        }
    };

    amrex::GpuArray<amrex::Real, 2> gamma_pInf = pp.cal_gamma(alpha1);
    Real gamma_bar = gamma_pInf[0];
    Real pInf_bar = gamma_pInf[1];

    Real rho = rho1 + rho2;
    return computeState(gamma_bar, pInf_bar, rho);
}

/**
 * 从原始变量中恢复状态方程相关物理量
 * @tparam eType 目标物理量类型
 * @param prim 原始变量
 * @param pp 两相物性参数
 * @return 根据 eType 返回总能量、声速或辅助量
 */
template <EOSType eType>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
prims_to_eos(const VarArray& prim, ProbParm const& pp)
{
    using namespace amrex;

    Real rho1 = prim[AlphaRho1];
    Real rho2 = prim[AlphaRho2];
    Real AMREX_D_DECL(u = prim[XMom], v = prim[YMom], w = prim[ZMom]);
    Real pressure = prim[Energy];
    Real alpha1 = prim[Alpha1];

    Real rho = rho1 + rho2;
    Real kineticEnergy = rho * (AMREX_D_TERM(u * u, +v * v, +w * w)) / 2;

    auto computeOutput = [&](Real gamma, Real pInf, Real rho, Real pressure) -> Real {
        if constexpr (eType == EOSType::TotalEnergy) {
            Real internalEnergy = (pressure + gamma * pInf) / (gamma - 1);
            return internalEnergy + kineticEnergy;
        } else if constexpr (eType == EOSType::SoundSpeed) {
            Real soundSpeedSqr = gamma * (pressure + pInf) / rho;
            return sqrt(soundSpeedSqr);
        } else if constexpr (eType == EOSType::SoundSpeedSqrDensity) {
            return gamma * (pressure + pInf);
        }
    };

    amrex::GpuArray<amrex::Real, 2> gamma_pInf = pp.cal_gamma(alpha1);
    Real gamma_bar = gamma_pInf[0];
    Real pInf_bar = gamma_pInf[1];

    return computeOutput(gamma_bar, pInf_bar, rho, pressure);
}

/**
 * 将两相守恒变量转换为原始变量
 * @param cons 守恒变量
 * @param pp 两相物性参数
 * @return 原始变量数组
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
VarArray cons_to_prim(VarArray cons, ProbParm const& pp)
{
    using namespace amrex;

    Real rho1 = cons[AlphaRho1];
    Real rho2 = cons[AlphaRho2];
    Real AMREX_D_DECL(rhou = cons[XMom], rhov = cons[YMom], rhow = cons[ZMom]);
    Real rhoE = cons[Energy];
    Real alpha1 = cons[Alpha1];

    Real rho = rho1 + rho2;
    Real AMREX_D_DECL(u = rhou / rho, v = rhov / rho, w = rhow / rho);

    Real pressure = cons_to_eos<Pressure>(cons, pp);

    return { rho1, rho2, AMREX_D_DECL(u, v, w), pressure, alpha1 };
}

/**
 * 将两相原始变量转换为守恒变量
 * @param prim 原始变量
 * @param pp 两相物性参数
 * @return 守恒变量数组
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
VarArray prim_to_cons(VarArray prim, ProbParm const& pp)
{
    using namespace amrex;

    Real rho1 = prim[AlphaRho1];
    Real rho2 = prim[AlphaRho2];
    Real AMREX_D_DECL(u = prim[XMom], v = prim[YMom], w = prim[ZMom]);
    Real pressure = prim[Pressure];
    Real alpha1 = prim[Alpha1];

    Real rho = rho1 + rho2;
    Real rhoE = prims_to_eos<TotalEnergy>(prim, pp);
    Real AMREX_D_DECL(rhou = rho * u, rhov = rho * v, rhow = rho * w);

    return { rho1, rho2, AMREX_D_DECL(rhou, rhov, rhow), rhoE, alpha1 };
}

/**
 * 在给定盒子范围内批量执行守恒量到原始量的转换
 * @param ghost_box 需要填充的区域，通常包含 ghost cell
 * @param prim_array 原始变量输出数组
 * @param statein_array 守恒变量输入数组
 * @param pp 两相物性参数
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void cons_to_prim_for_array4(const amrex::Box& ghost_box,
    amrex::Array4<amrex::Real> const& prim_array,
    amrex::Array4<amrex::Real const> const& statein_array,
    ProbParm const& pp = ProbParm())
{
    ParallelFor(ghost_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        int ncomp = NVar;
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

/**
 * 在给定法向上计算两相守恒通量
 * @param prim 原始变量
 * @param pp 两相物性参数
 * @param norm 面法向量
 * @return 与状态变量顺序一致的通量数组
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
VarArray
prim_to_flux(const VarArray& prim, ProbParm const& pp, amrex::GpuArray<amrex::Real, 3> norm)
{
    using namespace amrex;

    Real rho1 = prim[AlphaRho1];
    Real rho2 = prim[AlphaRho2];
    Real AMREX_D_DECL(u = prim[XMom], v = prim[YMom], w = prim[ZMom]);
    Real pressure = prim[Energy];
    Real alpha1 = prim[Alpha1];

    Real rho = rho1 + rho2;
    Real Vn = AMREX_D_TERM(u * norm[0], +v * norm[1], +w * norm[2]);
    Real rhoE = prims_to_eos<TotalEnergy>(prim, pp);

    VarArray flux = {
        rho1 * Vn,
        rho2 * Vn,
        AMREX_D_DECL(
            rho * u * Vn + pressure * norm[0],
            rho * v * Vn + pressure * norm[1],
            rho * w * Vn + pressure * norm[2]
        ),
        (rhoE + pressure) * Vn,
        0.
    };

    return flux;
}

/**
 * 根据体积分数梯度估计界面法向量
 * 
 * 参数约束：
 * - dx 网格间距必须为正数，且不为零
 * - primin_array 必须包含有效的体积分数数据
 * 
 * @param i x 方向索引
 * @param j y 方向索引
 * @param k z 方向索引
 * @param primin_array 原始变量数组
 * @param dx 网格间距
 * @return 单位法向量；若梯度退化则返回零向量
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real, 3>
compute_normal_vector1(int i, int j, int k,
    amrex::Array4<amrex::Real const> const& primin_array,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
{
    if constexpr (BL_SPACEDIM == 1) {
        return { 1.0, 0.0, 0.0 };
    }

    amrex::Real dpsi[3] = { 0.0 };

    AMREX_D_TERM(
        dpsi[0] = (primin_array(i + 1, j, k, Alpha1) - primin_array(i - 1, j, k, Alpha1)) / (2.0 * dx[0]);,
        dpsi[1] = (primin_array(i, j + 1, k, Alpha1) - primin_array(i, j - 1, k, Alpha1)) / (2.0 * dx[1]);,
        dpsi[2] = (primin_array(i, j, k + 1, Alpha1) - primin_array(i, j, k - 1, Alpha1)) / (2.0 * dx[2]););

    amrex::Real norm_factor = 0.0;
    for (int dim = 0; dim < BL_SPACEDIM; ++dim) {
        norm_factor += dpsi[dim] * dpsi[dim];
    }
    norm_factor = std::sqrt(norm_factor);

    amrex::GpuArray<amrex::Real, 3> n_l = { 0.0 };
    if (norm_factor > 1e-8) {
        for (int dim = 0; dim < BL_SPACEDIM; ++dim) {
            n_l[dim] = -dpsi[dim] / norm_factor;
        }
    }

    return n_l;
}

/**
 * 在给定方向上求解两相界面的近似 Riemann 问题
 * 
 * 算法说明：
 * - 使用 HLLC (Harten-Lax-van Leer Contact) 通量计算方法
 * - 计算左右状态的波速估计
 * - 根据波速确定中间状态
 * - 计算最终的界面通量
 * - 包含激波稳定方案，提高计算稳定性
 * 
 * @tparam idim 空间方向编号
 * @param i x 方向索引
 * @param j y 方向索引
 * @param k z 方向索引
 * @param prim_array 原始变量数组，用于局部平滑度与辅助判据计算
 * @param stateL 界面左状态
 * @param stateR 界面右状态
 * @param norm 面法向量
 * @param pp 两相物性参数
 * @param dx 网格间距
 * @return 界面通量
 */
template <int idim>
AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
VarArray
Riemann_solver_1(int i, int j, int k, const amrex::Array4<const amrex::Real>& prim_array, VarArray const& stateL, VarArray const& stateR, amrex::GpuArray<amrex::Real, 3> norm, ProbParm const& pp, amrex::GpuArray<amrex::Real, BL_SPACEDIM> dx)
{
    using namespace amrex;
    enum VelComponents {
        Ar1 = 0,
        Ar2,
        AMREX_D_DECL(U, V, W),
        Pressure,
        Alpha1
    };
    Long ncomp = NVar;
    enum { L, R };
    
    Real ar1l = stateL[Ar1], ar2l = stateL[Ar2];
    Real ar1r = stateR[Ar1], ar2r = stateR[Ar2];
    Real rl = ar1l + ar2l;
    Real rr = ar1r + ar2r;
    Real pl = stateL[Pressure], pr = stateR[Pressure];
    Real alpha1l = stateL[Alpha1], alpha1r = stateR[Alpha1];

    Real AMREX_D_DECL(ul = stateL[U], vl = stateL[V], wl = stateL[W]);
    Real AMREX_D_DECL(uR = stateR[U], vR = stateR[V], wR = stateR[W]);

    Real kappa = 0.0;
    Real alpha_avg = 0.5 * (alpha1l + alpha1r);
    Real S_switch = alpha_avg * (1.0 - alpha_avg);

    Real u_comp = kappa;
    Real n_dir = (alpha1r > alpha1l) ? 1.0 : -1.0;

    GpuArray<Real, 2> Vn = { AMREX_D_TERM(ul * norm[0] + u_comp * n_dir, +vl * norm[1], +wl * norm[2]),
        AMREX_D_TERM(uR * norm[0] + u_comp * n_dir, +vR * norm[1], +wR * norm[2]) };
    Real cl = prims_to_eos<SoundSpeed>(stateL, pp);
    Real cr = prims_to_eos<SoundSpeed>(stateR, pp);

#if BL_SPACEDIM == 2
    Real hll_flag = false;
    Long norm_velocity_direction;
    Dim3 offset_along, offset_norm;

    if (norm[0] == 1.0) {
        norm_velocity_direction = YMom;
        offset_along = amrex::IntVect::TheDimensionVector(0).dim3();
        offset_norm = amrex::IntVect::TheDimensionVector(1).dim3();
    } else {
        norm_velocity_direction = XMom;
        offset_along = amrex::IntVect::TheDimensionVector(1).dim3();
        offset_norm = amrex::IntVect::TheDimensionVector(0).dim3();
    }

    auto cal_numerator = [=](int i, int j, int k) -> Real {
        Real delta_v_up = prim_array(i + offset_norm.x, j + offset_norm.y, k + offset_norm.z, norm_velocity_direction)
            - prim_array(i, j, k, norm_velocity_direction);
        Real delta_v_down = prim_array(i, j, k, norm_velocity_direction)
            - prim_array(i - offset_norm.x, j - offset_norm.y, k - offset_norm.z, norm_velocity_direction);
        return std::min(std::abs(delta_v_up), std::abs(delta_v_down));
    };

    Real alpha_ip1 = cal_numerator(i, j, k) / (cr + std::max(uR, vR));
    Real alpha_i = cal_numerator(i - offset_along.x, j - offset_along.y, k - offset_along.z) / (cl + std::max(ul, vl));

    hll_flag = std::min(alpha_ip1, alpha_i) - 0.1 > 0;
#endif

    Real sqrt_rl = std::sqrt(rl);
    Real sqrt_rr = std::sqrt(rr);

    Real VnBar = (Vn[L] * sqrt_rl + Vn[R] * sqrt_rr) / (sqrt_rl + sqrt_rr);
    Real cBar = sqrt((cl * cl * sqrt_rl + cr * cr * sqrt_rr) / (sqrt_rl + sqrt_rr)
        + 0.5 * (sqrt_rl * sqrt_rr) / ((sqrt_rl + sqrt_rr) * (sqrt_rl + sqrt_rr)) * (Vn[R] - Vn[L]) * (Vn[R] - Vn[L]));

    Real Ma_local = std::max(std::abs(Vn[L] / cl), std::abs(Vn[R] / cr)), Ma_limit = 0.1;
    Real phi = sin(std::min(1., Ma_local / Ma_limit) * M_PI_2);
    Real SL = std::min(Vn[L] - cl, VnBar - cBar);
    Real SR = std::max(Vn[R] + cr, VnBar + cBar);


    amrex::Real Sstar = ((pr - pl) + (rl * Vn[L] * (SL - Vn[L]) - rr * Vn[R] * (SR - Vn[R]))) / (rl * (SL - Vn[L]) - rr * (SR - Vn[R]));

    VarArray res;
    // HLLC 通量计算
    if (SL >= 0) {
        res = prim_to_flux(stateL, pp, norm);
        res[Alpha1] = 0;
    } else if (Sstar > 0) {
        VarArray UL = prim_to_cons(stateL, pp);
        VarArray FL = prim_to_flux(stateL, pp, norm);
        amrex::Real Usf = (SL - Vn[L]) / (SL - Sstar);
        VarArray Ustar = { Usf * ar1l,
            Usf * ar2l,
            AMREX_D_DECL(Usf * (ul + norm[0] * (Sstar - Vn[L])) * rl,
                Usf * (vl + norm[1] * (Sstar - Vn[L])) * rl,
                Usf * (wl + norm[2] * (Sstar - Vn[L])) * rl),
            Usf * (UL[Pressure] + (Sstar - Vn[L]) * (rl * Sstar + pl / (SL - Vn[L]))),
            alpha1l };
        for (int ivar = 0; ivar < ncomp - 1; ivar++) {
            res[ivar] = FL[ivar] + SL * (Ustar[ivar] - UL[ivar]);
        }
#if BL_SPACEDIM == 2
        if (hll_flag) {
            const VarArray UR = prim_to_cons(stateR, pp);
            const VarArray FR = prim_to_flux(stateR, pp, norm);

            const auto correct = [&](const std::initializer_list<int>& vars) {
                for (int ivar : vars) {
                    res[ivar] = (SR * FL[ivar] - SL * FR[ivar] + SL * SR * (UR[ivar] - UL[ivar])) / (SR - SL);
                }
            };

            if (norm[0] == 1.0) {
                correct({ Ar1, Ar2, V });
            } else {
                correct({ Ar1, Ar2, U });
            }
        }
#endif
        res[Alpha1] = 0;
    } else if (SR > 0) {
        VarArray UR = prim_to_cons(stateR, pp);
        VarArray FR = prim_to_flux(stateR, pp, norm);
        amrex::Real Usf = (SR - Vn[R]) / (SR - Sstar);
        VarArray Ustar = { Usf * ar1r,
            Usf * ar2r,
            AMREX_D_DECL(Usf * (uR + norm[0] * (Sstar - Vn[R])) * rr,
                Usf * (vR + norm[1] * (Sstar - Vn[R])) * rr,
                Usf * (wR + norm[2] * (Sstar - Vn[R])) * rr),
            Usf * (UR[Pressure] + (Sstar - Vn[R]) * (rr * Sstar + pr / (SR - Vn[R]))),
            alpha1r };
        for (int ivar = 0; ivar < ncomp - 1; ivar++) {
            res[ivar] = FR[ivar] + SR * (Ustar[ivar] - UR[ivar]);
        }
#if BL_SPACEDIM == 2
        if (hll_flag) {
            const VarArray UL = prim_to_cons(stateL, pp);
            const VarArray FL = prim_to_flux(stateL, pp, norm);

            const auto correct = [&](const std::initializer_list<int>& vars) {
                for (int ivar : vars) {
                    res[ivar] = (SR * FL[ivar] - SL * FR[ivar] + SL * SR * (UR[ivar] - UL[ivar])) / (SR - SL);
                }
            };

            if (norm[0] == 1.0) {
                correct({ Ar1, Ar2, V });
            } else {
                correct({ Ar1, Ar2, U });
            }
        }
#endif
        res[Alpha1] = 0;
    } else {
        res = prim_to_flux(stateR, pp, norm);
        res[Alpha1] = 0;
    }

    return res;
}

/**
 * 求解两相界面的近似 Riemann 问题（简化版本）
 * 
 * 算法说明：
 * - Riemann_solver_1 的简化版本
 * - 使用 HLLC (Harten-Lax-van Leer Contact) 通量计算方法
 * - 计算左右状态的波速估计
 * - 根据波速确定中间状态
 * - 计算最终的界面通量
 * - 移除了激波稳定方案，适用于简单场景
 * 
 * @param i x 方向索引
 * @param j y 方向索引
 * @param k z 方向索引
 * @param prim_array 原始变量数组
 * @param stateL 界面左状态
 * @param stateR 界面右状态
 * @param norm 面法向量
 * @param pp 两相物性参数
 * @return 界面通量
 */
AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
VarArray
Riemann_solver_2(int i, int j, int k, const amrex::Array4<const amrex::Real>& prim_array, VarArray const& stateL, VarArray const& stateR, amrex::GpuArray<amrex::Real, 3> norm, ProbParm const& pp)
{
    using namespace amrex;
    enum VelComponents {
        Ar1 = 0,
        Ar2,
        AMREX_D_DECL(U, V, W),
        Pressure,
        Alpha1
    };
    Long ncomp = NVar;
    enum { L, R };
    
    Real ar1l = stateL[Ar1], ar2l = stateL[Ar2];
    Real ar1r = stateR[Ar1], ar2r = stateR[Ar2];
    Real rl = ar1l + ar2l;
    Real rr = ar1r + ar2r;

    Real AMREX_D_DECL(ul = stateL[U], vl = stateL[V], wl = stateL[W]);
    Real AMREX_D_DECL(uR = stateR[U], vR = stateR[V], wR = stateR[W]);

    GpuArray<Real, 2> Vn = { AMREX_D_TERM(ul * norm[0], +vl * norm[1], +wl * norm[2]),
        AMREX_D_TERM(uR * norm[0], +vR * norm[1], +wR * norm[2]) };
    Real cl = prims_to_eos<SoundSpeed>(stateL, pp);
    Real cr = prims_to_eos<SoundSpeed>(stateR, pp);

    Real pl = stateL[Pressure], pr = stateR[Pressure];
    Real alpha1l = stateL[Alpha1], alpha1r = stateR[Alpha1];
    
    Real sqrt_rl = std::sqrt(rl);
    Real sqrt_rr = std::sqrt(rr);

    Real VnBar = (Vn[L] * sqrt_rl + Vn[R] * sqrt_rr) / (sqrt_rl + sqrt_rr);
    Real cBar = sqrt((cl * cl * sqrt_rl + cr * cr * sqrt_rr) / (sqrt_rl + sqrt_rr)
        + 0.5 * (sqrt_rl * sqrt_rr) / ((sqrt_rl + sqrt_rr) * (sqrt_rl + sqrt_rr)) * (Vn[R] - Vn[L]) * (Vn[R] - Vn[L]));

    Real Ma_local = std::max(std::abs(Vn[L] / cl), std::abs(Vn[R] / cr)), Ma_limit = 0.1;
    Real phi = sin(std::min(1., Ma_local / Ma_limit) * M_PI_2);
    Real SL = std::min(Vn[L] - cl, VnBar - cBar);
    Real SR = std::max(Vn[R] + cr, VnBar + cBar);

    amrex::Real Sstar = ((pr - pl) + (rl * Vn[L] * (SL - Vn[L]) - rr * Vn[R] * (SR - Vn[R]))) / (rl * (SL - Vn[L]) - rr * (SR - Vn[R]));
    // HLLC 通量计算
    if (SL >= 0) {
        VarArray res = prim_to_flux(stateL, pp, norm);
        res[Alpha1] = 0;
        return res;
    } else if (Sstar > 0) {
        VarArray UL = prim_to_cons(stateL, pp);
        VarArray FL = prim_to_flux(stateL, pp, norm);
        amrex::Real Usf = (SL - Vn[L]) / (SL - Sstar);
        VarArray Ustar = { Usf * ar1l,
            Usf * ar2l,
            AMREX_D_DECL(Usf * (ul + norm[0] * (Sstar - Vn[L])) * rl,
                Usf * (vl + norm[1] * (Sstar - Vn[L])) * rl,
                Usf * (wl + norm[2] * (Sstar - Vn[L])) * rl),
            Usf * (UL[Pressure] + (Sstar - Vn[L]) * (rl * Sstar + pl / (SL - Vn[L]))),
            alpha1l };
        VarArray res;
        for (int ivar = 0; ivar < ncomp - 1; ivar++) {
            res[ivar] = FL[ivar] + SL * (Ustar[ivar] - UL[ivar]);
        }
        res[Alpha1] = 0;
        return res;
    } else if (SR > 0) {
        VarArray UR = prim_to_cons(stateR, pp);
        VarArray FR = prim_to_flux(stateR, pp, norm);
        amrex::Real Usf = (SR - Vn[R]) / (SR - Sstar);
        VarArray Ustar = { Usf * ar1r,
            Usf * ar2r,
            AMREX_D_DECL(Usf * (uR + norm[0] * (Sstar - Vn[R])) * rr,
                Usf * (vR + norm[1] * (Sstar - Vn[R])) * rr,
                Usf * (wR + norm[2] * (Sstar - Vn[R])) * rr),
            Usf * (UR[Pressure] + (Sstar - Vn[R]) * (rr * Sstar + pr / (SR - Vn[R]))),
            alpha1r };
        VarArray res;
        for (int ivar = 0; ivar < ncomp - 1; ivar++) {
            res[ivar] = FR[ivar] + SR * (Ustar[ivar] - UR[ivar]);
        }
        res[Alpha1] = 0;
        return res;
    } else {
        VarArray res = prim_to_flux(stateR, pp, norm);
        res[Alpha1] = 0;
        return res;
    }
}

/**
 * 计算每个方向的最大特征速度
 * @param i x 方向索引
 * @param j y 方向索引
 * @param k z 方向索引
 * @param state 状态变量数组
 * @param umax 最大特征速度输出数组
 * @param pp 两相物性参数
 * @param dx 网格间距
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void cal_u_max(
    int i, int j, int k,
    amrex::Array4<amrex::Real const> const& state,
    amrex::Array4<amrex::Real> const& umax,
    ProbParm const& pp,
    amrex::GpuArray<amrex::Real, BL_SPACEDIM> dx) noexcept
{
    using namespace amrex;

    VarArray cons;
    for (int ivar = 0; ivar < NVar; ivar++) {
        cons[ivar] = state(i, j, k, ivar);
    }
    VarArray prim = cons_to_prim(cons, pp);

    Real c = cons_to_eos<SoundSpeed>(cons, pp);

    AMREX_D_TERM(
        amrex::Real lambda_x = (c + std::abs(prim[XMom])) / dx[0];,
        amrex::Real lambda_y = (c + std::abs(prim[YMom])) / dx[1];
        ,
        amrex::Real lambda_z = (c + std::abs(prim[ZMom])) / dx[2];)

    AMREX_D_TERM(
        umax(i, j, k, 0) = lambda_x;,
        umax(i, j, k, 1) = lambda_y;
        ,
        umax(i, j, k, 2) = lambda_z;)
}



