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

#if BL_SPACEDIM == 1
/**
 * 初始化一维两相算例的守恒量场
 * @param level 网格层级
 * @param time 当前时间
 * @param box 网格盒子
 * @param state 状态变量数组
 * @param dx 网格间距
 * @param problo 计算域左下角坐标
 * @param pp 两相物性参数
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE void
initdata(int level, amrex::Real time,
    const amrex::Box& box, amrex::FArrayBox& state,
    const amrex::Real* dx, const amrex::Real* problo, const ProbParm& pp)
{
    using namespace amrex;

    // 原始变量数组
    VarArray prim;

    // 创建状态数组的 Array4 视图
    Array4<Real> state_arr = state.array();

    // 遍历网格框中的所有单元格（一维情况下仅 x 方向）
    const int i_lo = box.smallEnd(0);
    const int i_hi = box.bigEnd(0);

    for (int i = i_lo; i <= i_hi; ++i) {
        // 计算单元格中心在 x 方向的物理坐标
        Real x = problo[0] + (i + 0.5) * dx[0];

        // 根据选择的算例类型设置原始变量
        switch (pp.case_type) {
            case CASE_1D_BASE: {
                //【算例】1：基础一维两相流配置
                // 描述：中心区域为水，两侧为空气，整体以相同速度运动
                constexpr Real epss = 1e-8;
                if (x >= 0.25 && x < 0.75) {
                    // 区域：0.25 <= x < 0.75（水）
                    prim[AlphaRho1] = 1000.0; // alpha1 * rho1 (kg/m^3)
                    prim[AlphaRho2] = epss; // alpha2 * rho2 (kg/m^3)
                    prim[XMom] = 100.0;       // u (m/s)
                    prim[Energy] = 101325.0;   // p (Pa)
                    prim[Alpha1] = 1.0 - epss; // alpha1
                } else {
                    // 区域：其他（空气）
                    prim[AlphaRho1] = epss;  // alpha1 * rho1 (kg/m^3)
                    prim[AlphaRho2] = 1.204;   // alpha2 * rho2 (kg/m^3)
                    prim[XMom] = 100.0;        // u (m/s)
                    prim[Energy] = 101325.0;    // p (Pa)
                    prim[Alpha1] = epss;     // alpha1
                }
                break;
            }
            case CASE_1D_SHOCK_TUBE: {
                //【算例】2：二维两相激波管（2phase sod tube）
                // 描述：左侧为高压水，右侧为低压空气，模拟激波管问题
                constexpr Real epss = 1e-8;
                if (x < 0.8) {
                    // 区域：x < 0.8（高压水）
                    prim[AlphaRho1] = 1000.0; // alpha1 * rho1 (kg/m^3)
                    prim[AlphaRho2] = epss;   // alpha2 * rho2 (kg/m^3)
                    prim[XMom] = 0.0;         // u (m/s)
                    prim[Energy] = 1.0e9;     // p (Pa) - 高压
                    prim[Alpha1] = 1.0 - epss; // alpha1
                } else {
                    // 区域：x >= 0.8（低压空气）
                    prim[AlphaRho1] = epss;    // alpha1 * rho1 (kg/m^3)
                    prim[AlphaRho2] = 20.0;    // alpha2 * rho2 (kg/m^3)
                    prim[XMom] = 0.0;          // u (m/s)
                    prim[Energy] = 1.0e5;      // p (Pa) - 低压
                    prim[Alpha1] = epss;       // alpha1
                }
                break;
            }
            case CASE_1D_SHOCK_INTERFACE: {
                //【算例】3：激波与静止界面相互作用（shock static interface）
                // 描述：三段区域配置，包含高速水流、静止水和空气
                constexpr Real epss = 1e-8;
                if (x > 0.0 && x < 1.31) {
                    // 区域：0.0 < x < 1.31（高速水流）
                    prim[AlphaRho1] = 1230.0; // rho1 * alpha1 (kg/m^3)
                    prim[AlphaRho2] = epss;   // rho2 * alpha2 (kg/m^3)
                    prim[XMom] = 432.7;       // u (m/s) - 高速
                    prim[Energy] = 1.0e9;     // p (Pa) - 高压
                    prim[Alpha1] = 1.0 - epss; // alpha1
                } else if (x > 1.31 && x < 2.0) {
                    // 区域：1.31 < x < 2.0（静止水）
                    prim[AlphaRho1] = 1000.0; // rho1 * alpha1 (kg/m^3)
                    prim[AlphaRho2] = epss; // rho2 * alpha2 (kg/m^3)
                    prim[XMom] = 0.0;         // u (m/s) - 静止
                    prim[Energy] = 1.0e5;     // p (Pa) - 常压
                    prim[Alpha1] = 1.0 - epss; // alpha1
                } else {
                    // 区域：其他（空气）
                    prim[AlphaRho1] = epss;    // rho1 * alpha1 (kg/m^3)
                    prim[AlphaRho2] = 1.0;     // rho2 * alpha2 (kg/m^3)
                    prim[XMom] = 0.0;          // u (m/s) - 静止
                    prim[Energy] = 1.0e5;      // p (Pa) - 常压
                    prim[Alpha1] = epss;       // alpha1
                }
                break;
            }
            default: {
                // 默认算例：二维两相激波管
                constexpr Real epss = 1e-8;
                if (x < 0.8) {
                    // 区域：x < 0.8（高压水）
                    prim[AlphaRho1] = 1000.0; // alpha1 * rho1 (kg/m^3)
                    prim[AlphaRho2] = epss;   // alpha2 * rho2 (kg/m^3)
                    prim[XMom] = 0.0;         // u (m/s)
                    prim[Energy] = 1.0e9;     // p (Pa) - 高压
                    prim[Alpha1] = 1.0 - epss; // alpha1
                } else {
                    // 区域：x >= 0.8（低压空气）
                    prim[AlphaRho1] = epss;    // alpha1 * rho1 (kg/m^3)
                    prim[AlphaRho2] = 20.0;    // alpha2 * rho2 (kg/m^3)
                    prim[XMom] = 0.0;          // u (m/s)
                    prim[Energy] = 1.0e5;      // p (Pa) - 低压
                    prim[Alpha1] = epss;       // alpha1
                }
                break;
            }
        }

        // 将原始变量转换为守恒变量
        VarArray cons = prim_to_cons(prim, pp);
        // 使用 Array4 将守恒变量存储到状态数组中
        state_arr(i, 0, 0, AlphaRho1) = cons[AlphaRho1];
        state_arr(i, 0, 0, AlphaRho2) = cons[AlphaRho2];
        state_arr(i, 0, 0, XMom) = cons[XMom];
        state_arr(i, 0, 0, Energy) = cons[Energy];
        state_arr(i, 0, 0, Alpha1) = cons[Alpha1];
    }
}
#endif

#if BL_SPACEDIM == 2
/**
 * 初始化二维两相算例的守恒量场
 * @param level 网格层级
 * @param time 当前时间
 * @param box 网格盒子
 * @param state 状态变量数组
 * @param dx 网格间距
 * @param problo 计算域左下角坐标
 * @param pp 两相物性参数
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE void
initdata(int level, amrex::Real time,
    const amrex::Box& box, amrex::FArrayBox& state,
    const amrex::Real* dx, const amrex::Real* problo, const ProbParm& pp)
{
    using namespace amrex;

    // 原始变量数组
    VarArray prim;

    // 创建状态数组的 Array4 视图
    Array4<Real> state_arr = state.array();

    // 遍历网格框中的所有单元格（二维情况）
    const int i_lo = box.smallEnd(0);
    const int i_hi = box.bigEnd(0);
    const int j_lo = box.smallEnd(1);
    const int j_hi = box.bigEnd(1);

    for (int j = j_lo; j <= j_hi; ++j) {
        for (int i = i_lo; i <= i_hi; ++i) {
            // 计算单元格中心在 x 和 y 方向的物理坐标
            Real x = problo[0] + (i + 0.5) * dx[0];
            Real y = problo[1] + (j + 0.5) * dx[1];

            // 根据选择的算例类型设置原始变量
            switch (pp.case_type) {
                case CASE_2D_BASE: {
                    //【算例】1：基础二维水柱配置
                    // 描述：简单的水柱和空气区域划分，无激波
                    constexpr Real water_cylinder_x_center = 0.04; // 水柱中心 x 坐标 [m]
                    constexpr Real water_cylinder_y_center = 0.00; // 水柱中心 y 坐标 [m]
                    constexpr Real water_cylinder_radius = 0.011; // 水柱半径 [m]
                    constexpr Real epss = 1e-8;
                    
                    Real distance = std::sqrt((x - water_cylinder_x_center) * (x - water_cylinder_x_center) + 
                                           (y - water_cylinder_y_center) * (y - water_cylinder_y_center));
                                           
                    if (distance <= water_cylinder_radius) {
                        // 区域：水柱内部
                        prim[AlphaRho1] = 1000.0; // α1ρ1 (kg/m^3)
                        prim[AlphaRho2] = epss; // α2ρ2 (kg/m^3)
                        prim[XMom] = 0.0;         // u (m/s)
                        prim[YMom] = 0.0;         // v (m/s)
                        prim[Energy] = 1.01e5;     // p (Pa)
                        prim[Alpha1] = 1.0 - epss; // α1
                    } else { 
                        // 区域：水柱外部（空气）
                        prim[AlphaRho1] = epss;  // α1ρ1 (kg/m^3)
                        prim[AlphaRho2] = 1.17;    // α2ρ2 (kg/m^3)
                        prim[XMom] = 0.0;          // u (m/s)
                        prim[YMom] = 0.0;          // v (m/s)
                        prim[Energy] = 1.01e5;      // p (Pa)
                        prim[Alpha1] = epss;     // α1
                    }
                    break;
                }
                case CASE_2D_MACH2_4_SHOCK: {
                    //【算例】2：Mach 2.4 激波撞击水柱
                    // 描述：左侧为高速高压空气（激波后），中间为水柱，右侧为静止空气
                    // 几何参数：水柱中心 (0.04, 0.00) m，半径 0.011 m
                    // 激波位置：x = 0.029 m
                    constexpr Real water_cylinder_x_center = 0.04; // 水柱中心 x 坐标 [m]
                    constexpr Real water_cylinder_y_center = 0.00; // 水柱中心 y 坐标 [m]
                    constexpr Real water_cylinder_radius = 0.011; // 水柱半径 [m]
                    constexpr Real epss = 1e-8;

                    Real distance = std::sqrt((x - water_cylinder_x_center) * (x - water_cylinder_x_center) + 
                                           (y - water_cylinder_y_center) * (y - water_cylinder_y_center));
                                           
                    if (x <= 0.029) {
                        // 区域：x <= 0.029（激波后空气）
                        prim[AlphaRho1] = epss;    // α1ρ1 (kg/m^3)
                        prim[AlphaRho2] = 3.7579;  // α2ρ2 (kg/m^3)
                        prim[XMom] = 574.57;       // u (m/s) - 高速
                        prim[YMom] = 0.0;          // v (m/s)
                        prim[Energy] = 6.6189e5;    // p (Pa) - 高压
                        prim[Alpha1] = epss;       // α1
                    } else if (distance <= water_cylinder_radius + 1e-8) {
                        // 区域：水柱内部
                        prim[AlphaRho1] = 1000.0; // α1ρ1 (kg/m^3)
                        prim[AlphaRho2] = epss;   // α2ρ2 (kg/m^3)
                        prim[XMom] = 0.0;         // u (m/s) - 静止
                        prim[YMom] = 0.0;         // v (m/s) - 静止
                        prim[Energy] = 1.01e5;     // p (Pa) - 常压
                        prim[Alpha1] = 1.0 - epss; // α1
                    } else {
                        // 区域：水柱外部（静止空气）
                        prim[AlphaRho1] = epss;    // α1ρ1 (kg/m^3)
                        prim[AlphaRho2] = 1.17;    // α2ρ2 (kg/m^3)
                        prim[XMom] = 0.0;          // u (m/s) - 静止
                        prim[YMom] = 0.0;          // v (m/s) - 静止
                        prim[Energy] = 1.01e5;      // p (Pa) - 常压
                        prim[Alpha1] = epss;       // α1
                    }
                    break;
                }
                case CASE_2D_MACH10_SHOCK: {
                    //【算例】3：马赫 10 激波撞击水柱
                    // 描述：左侧为极高速度高压空气（激波后），中间为小水柱，右侧为静止空气
                    // 几何参数：水柱中心 (0.013, 0.00) m，半径 0.004 m
                    // 激波位置：x = 0.009 m
                    constexpr Real water_cylinder_x_center = 0.013; // 水柱中心 x 坐标 [m]
                    constexpr Real water_cylinder_y_center = 0.00; // 水柱中心 y 坐标 [m]
                    constexpr Real water_cylinder_radius = 0.004;  // 水柱半径 [m]
                    constexpr Real epss = 1e-8;

                    Real distance = std::sqrt((x - water_cylinder_x_center) * (x - water_cylinder_x_center) + 
                                           (y - water_cylinder_y_center) * (y - water_cylinder_y_center));
                                           
                    if (distance <= water_cylinder_radius) {
                        // 区域：水柱内部
                        prim[AlphaRho1] = 1000.0; // α1ρ1 (kg/m^3)
                        prim[AlphaRho2] = epss;   // α2ρ2 (kg/m^3)
                        prim[XMom] = 0.0;         // u (m/s) - 静止
                        prim[YMom] = 0.0;         // v (m/s) - 静止
                        prim[Energy] = 1.0e5;      // p (Pa) - 常压
                        prim[Alpha1] = 1.0 - epss; // α1
                    } else if (x > 0.009 - 1e-8) {
                        // 区域：x > 0.009（激波前空气）
                        prim[AlphaRho1] = epss;    // α1ρ1 (kg/m^3)
                        prim[AlphaRho2] = 1.2;     // α2ρ2 (kg/m^3)
                        prim[XMom] = 0.0;          // u (m/s) - 静止
                        prim[YMom] = 0.0;          // v (m/s) - 静止
                        prim[Energy] = 1.0e5;      // p (Pa) - 常压
                        prim[Alpha1] = epss;       // α1
                    } else {
                        // 区域：x <= 0.009（激波后空气）
                        prim[AlphaRho1] = epss;     // α1ρ1 (kg/m^3)
                        prim[AlphaRho2] = 6.8571;   // α2ρ2 (kg/m^3)
                        prim[XMom] = 2.8179e3;      // u (m/s) - 极高速度
                        prim[YMom] = 0.0;           // v (m/s)
                        prim[Energy] = 1.165e7;     // p (Pa) - 极高压力
                        prim[Alpha1] = epss;        // α1
                    }
                    break;
                }
                case CASE_2D_RMI_SF6_CASE1: {
                    // 【算例】4： RMI_SF6 case1-半边截断
                    // 描述：2D Richtmyer–Meshkov 不稳定性（空气–SF6）
                    // 激波位置：x = 3.2（激波从右向左传播）
                    constexpr Real x_shock = 3.2;          // 激波面位置 x_shock
                    constexpr Real x_1 = 2.9;          // 基准位置 x_1
                    constexpr Real varepsilon = 0.2;          // 扰动振幅 varepsilon
                    constexpr Real epss = 1e-8;            // 小量用于避免浮点误差

                    // 计算扰动界面 x 坐标（作为 y 的函数）
                    Real x_interface = x_1 - varepsilon * std::sin(2.0 * M_PI * (y + 0.25));

                    // 判断空间位置 → 三区域划分
                    if (x < x_interface) {
                        // 区域 1: SF6（重流体，左侧）
                        prim[AlphaRho1] = epss;          // α₁ρ₁ ≈ 0 (空气几乎为0)
                        prim[AlphaRho2] = 5.04;          // α₂ρ₂ = ρ_SF6
                        prim[XMom]      = 0.0;           // u = 0
                        prim[YMom]      = 0.0;           // v = 0
                        prim[Energy]    = 1.0;           // p = 1.0
                        prim[Alpha1]    = epss;          // α₁ ≈ 0
                    } else if (x > x_shock + epss) {
                        // 区域 2: 已受冲击空气（激波后，右侧）
                        prim[AlphaRho1] = 1.411;         // α₁ρ₁ = ρ_post_air
                        prim[AlphaRho2] = epss;          // α₂ρ₂ ≈ 0
                        prim[XMom]      = -0.39;         // u = -0.39（向左）
                        prim[YMom]      = 0.0;           // v = 0
                        prim[Energy]    = 1.628;         // p = 1.628
                        prim[Alpha1]    = 1.0 - epss;    // α₁ ≈ 1
                    } else {
                        // 区域 3: 未受冲击空气（激波前，介于界面与激波面之间）
                        prim[AlphaRho1] = 1.0;           // α₁ρ₁ = ρ_air
                        prim[AlphaRho2] = epss;          // α₂ρ₂ ≈ 0
                        prim[XMom]      = 0.0;           // u = 0
                        prim[YMom]      = 0.0;           // v = 0
                        prim[Energy]    = 1.0;           // p = 1.0
                        prim[Alpha1]    = 1.0 - epss;    // α₁ ≈ 1
                    }
                    break;
                }
                case CASE_2D_RMI_SF6_SQUARE: {
                    // 【算例】5：RMI_SF6 方块
                    // 描述：2D Richtmyer–Meshkov 不稳定性（SF6方块–空气）
                    // 激波马赫数：1.2
                    // 激波位置：0.2
                    // 相1，相2：  SF6，空气
                    constexpr Real x_shock = 0.2;     // 激波初始位置（左侧）
                    constexpr Real epss    = 1e-8;    // 小量

                    // -------------------------
                    // 理想方块界面参数
                    // -------------------------
                    constexpr Real x_c = 0.5;         // 方块中心 x 坐标
                    constexpr Real y_c = 0.0;         // 方块中心 y 坐标
                    constexpr Real L   = 0.4;         // 方块边长

                    // 理想方块内部判定
                    bool in_square = (std::abs(x - x_c) <= 0.5 * L) &&
                                    (std::abs(y - y_c) <= 0.5 * L);

                    // 分区：
                    // x < x_shock         -> 激波后空气（左侧）
                    // x_shock ~ square    -> 激波前空气
                    // in_square           -> SF6方块
                    if (in_square) {
                        // 区域 1: SF6方块内部（重流体）
                        prim[AlphaRho1] = 5.04;       // α₁ρ₁ = ρ_SF6 (相1现在是SF6)
                        prim[AlphaRho2] = epss;       // α₂ρ₂ ≈ 0
                        prim[XMom]      = 0.0;
                        prim[YMom]      = 0.0;
                        prim[Energy]    = 1.0;
                        prim[Alpha1]    = 1.0 - epss; // α₁ ≈ 1 (SF6体积分数接近1)
                    }
                    else if (x < x_shock - epss) {
                        // 区域 2: 已受冲击空气（激波后，左侧，轻流体）
                        prim[AlphaRho1] = epss;       // α₁ρ₁ ≈ 0 (SF6几乎为0)
                        prim[AlphaRho2] = 1.411;      // α₂ρ₂ = ρ_post_air (相2是空气)
                        prim[XMom]      = 0.39;       // 向右传播时，速度应为正
                        prim[YMom]      = 0.0;
                        prim[Energy]    = 1.628;
                        prim[Alpha1]    = epss;       // α₁ ≈ 0 (SF6体积分数接近0)
                    }
                    else {
                        // 区域 3: 未受冲击空气（激波前，轻流体）
                        prim[AlphaRho1] = epss;       // α₁ρ₁ ≈ 0
                        prim[AlphaRho2] = 1.0;        // α₂ρ₂ = ρ_air (相2是空气)
                        prim[XMom]      = 0.0;
                        prim[YMom]      = 0.0;
                        prim[Energy]    = 1.0;
                        prim[Alpha1]    = epss;       // α₁ ≈ 0
                    }
                    break;
                }
                default: {
                    // 默认算例：RMI_SF6 方块
                    constexpr Real x_shock = 0.2;     // 激波初始位置（左侧）
                    constexpr Real epss    = 1e-8;    // 小量

                    // -------------------------
                    // 理想方块界面参数
                    // -------------------------
                    constexpr Real x_c = 0.5;         // 方块中心 x 坐标
                    constexpr Real y_c = 0.0;         // 方块中心 y 坐标
                    constexpr Real L   = 0.4;         // 方块边长

                    // 理想方块内部判定
                    bool in_square = (std::abs(x - x_c) <= 0.5 * L) &&
                                    (std::abs(y - y_c) <= 0.5 * L);

                    // 分区：
                    // x < x_shock         -> 激波后空气（左侧）
                    // x_shock ~ square    -> 激波前空气
                    // in_square           -> SF6方块
                    if (in_square) {
                        // 区域 1: SF6方块内部（重流体）
                        prim[AlphaRho1] = 5.04;       // α₁ρ₁ = ρ_SF6 (相1现在是SF6)
                        prim[AlphaRho2] = epss;       // α₂ρ₂ ≈ 0
                        prim[XMom]      = 0.0;
                        prim[YMom]      = 0.0;
                        prim[Energy]    = 1.0;
                        prim[Alpha1]    = 1.0 - epss; // α₁ ≈ 1 (SF6体积分数接近1)
                    }
                    else if (x < x_shock - epss) {
                        // 区域 2: 已受冲击空气（激波后，左侧，轻流体）
                        prim[AlphaRho1] = epss;       // α₁ρ₁ ≈ 0 (SF6几乎为0)
                        prim[AlphaRho2] = 1.411;      // α₂ρ₂ = ρ_post_air (相2是空气)
                        prim[XMom]      = 0.39;       // 向右传播时，速度应为正
                        prim[YMom]      = 0.0;
                        prim[Energy]    = 1.628;
                        prim[Alpha1]    = epss;       // α₁ ≈ 0 (SF6体积分数接近0)
                    }
                    else {
                        // 区域 3: 未受冲击空气（激波前，轻流体）
                        prim[AlphaRho1] = epss;       // α₁ρ₁ ≈ 0
                        prim[AlphaRho2] = 1.0;        // α₂ρ₂ = ρ_air (相2是空气)
                        prim[XMom]      = 0.0;
                        prim[YMom]      = 0.0;
                        prim[Energy]    = 1.0;
                        prim[Alpha1]    = epss;       // α₁ ≈ 0
                    }
                    break;
                }
            }

            // 将原始变量转换为守恒变量
            VarArray cons = prim_to_cons(prim, pp);

            // 使用 Array4 将守恒变量存储到状态数组中
            state_arr(i, j, 0, AlphaRho1) = cons[AlphaRho1];
            state_arr(i, j, 0, AlphaRho2) = cons[AlphaRho2];
            state_arr(i, j, 0, XMom) = cons[XMom];
            state_arr(i, j, 0, YMom) = cons[YMom];
            state_arr(i, j, 0, Energy) = cons[Energy];
            state_arr(i, j, 0, Alpha1) = cons[Alpha1];
        }
    }
}
#endif

#if BL_SPACEDIM == 3
/**
 * 初始化三维两相算例的守恒量场
 * @param level 网格层级
 * @param time 当前时间
 * @param box 网格盒子
 * @param state 状态变量数组
 * @param dx 网格间距
 * @param problo 计算域左下角坐标
 * @param pp 两相物性参数
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE void
initdata(int level, amrex::Real time,
    const amrex::Box& box, amrex::FArrayBox& state,
    const amrex::Real* dx, const amrex::Real* problo, const ProbParm& pp)
{
    using namespace amrex;

    // 原始变量数组
    VarArray prim;

    // 创建状态数组的 Array4 视图
    Array4<Real> state_arr = state.array();

    // 遍历网格框中的所有单元格（三维情况）
    const int i_lo = box.smallEnd(0);
    const int i_hi = box.bigEnd(0);
    const int j_lo = box.smallEnd(1);
    const int j_hi = box.bigEnd(1);
    const int k_lo = box.smallEnd(2);
    const int k_hi = box.bigEnd(2);

    for (int k = k_lo; k <= k_hi; ++k) {
        for (int j = j_lo; j <= j_hi; ++j) {
            for (int i = i_lo; i <= i_hi; ++i) {
                // 计算单元格中心的物理坐标
                Real x = problo[0] + (i + 0.5) * dx[0];
                Real y = problo[1] + (j + 0.5) * dx[1];
                Real z = problo[2] + (k + 0.5) * dx[2];

                // 根据选择的算例类型设置原始变量
                switch (pp.case_type) {
                    case CASE_3D_RMI_SF6: {
                        //【算例】1：三维Richtmyer–Meshkov不稳定性
                        // 描述：3D Richtmyer–Meshkov 不稳定性（空气–SF6）
                        // 激波马赫数：
                        // 激波位置：x = 3.2
                        // 相1，相2：空气，SF6
                        constexpr Real x_shock = 3.2;          // 激波面位置 x_shock
                        constexpr Real x_1 = 2.9;          // 基准位置 x_1
                        constexpr Real varepsilon = 0.2;          // 扰动振幅 varepsilon
                        constexpr Real epss = 1e-8;            // 小量用于避免浮点误差
                        constexpr Real lamba = 1.0;            // 激波波速 lamba

                        // 计算扰动界面 x 坐标（作为 y 和 z 的函数，三维扰动）
                        Real x_interface = x_1 - varepsilon * (std::sin(2.0 * M_PI * (y + 0.25)/lamba) * std::sin(2.0 * M_PI * (z + 0.25)/lamba));

                        // 判断空间位置 → 三区域划分
                        if (x < x_interface) {
                            // 区域 1: SF6（重流体，左侧）
                            prim[AlphaRho1] = epss;          // α₁ρ₁ ≈ 0 (空气几乎为0)
                            prim[AlphaRho2] = 5.04;          // α₂ρ₂ = ρ_SF6
                            prim[XMom]      = 0.0;           // u = 0
                            prim[YMom]      = 0.0;           // v = 0
                            prim[ZMom]      = 0.0;           // w = 0
                            prim[Energy]    = 1.0;           // p = 1.0
                            prim[Alpha1]    = epss;          // α₁ ≈ 0
                        } else if (x > x_shock + epss) {
                            // 区域 2: 已受冲击空气（激波后，右侧）
                            prim[AlphaRho1] = 1.411;         // α₁ρ₁ = ρ_post_air
                            prim[AlphaRho2] = epss;          // α₂ρ₂ ≈ 0
                            prim[XMom]      = -0.39;         // u = -0.39（向左）
                            prim[YMom]      = 0.0;           // v = 0
                            prim[ZMom]      = 0.0;           // w = 0
                            prim[Energy]    = 1.628;         // p = 1.628
                            prim[Alpha1]    = 1.0 - epss;    // α₁ ≈ 1
                        } else {
                            // 区域 3: 未受冲击空气（激波前，介于界面与激波面之间）
                            prim[AlphaRho1] = 1.0;           // α₁ρ₁ = ρ_air
                            prim[AlphaRho2] = epss;          // α₂ρ₂ ≈ 0
                            prim[XMom]      = 0.0;           // u = 0
                            prim[YMom]      = 0.0;           // v = 0
                            prim[ZMom]      = 0.0;           // w = 0
                            prim[Energy]    = 1.0;           // p = 1.0
                            prim[Alpha1]    = 1.0 - epss;    // α₁ ≈ 1
                        }
                        break;
                    }
                    default: {
                        // 默认算例：三维Richtmyer–Meshkov不稳定性
                        constexpr Real x_shock = 3.2;          // 激波面位置 x_shock
                        constexpr Real x_1 = 2.9;          // 基准位置 x_1
                        constexpr Real varepsilon = 0.2;          // 扰动振幅 varepsilon
                        constexpr Real epss = 1e-8;            // 小量用于避免浮点误差
                        constexpr Real lamba = 1.0;            // 激波波速 lamba

                        // 计算扰动界面 x 坐标（作为 y 和 z 的函数，三维扰动）
                        Real x_interface = x_1 - varepsilon * (std::sin(2.0 * M_PI * (y + 0.25)/lamba) * std::sin(2.0 * M_PI * (z + 0.25)/lamba));

                        // 判断空间位置 → 三区域划分
                        if (x < x_interface) {
                            // 区域 1: SF6（重流体，左侧）
                            prim[AlphaRho1] = epss;          // α₁ρ₁ ≈ 0 (空气几乎为0)
                            prim[AlphaRho2] = 5.04;          // α₂ρ₂ = ρ_SF6
                            prim[XMom]      = 0.0;           // u = 0
                            prim[YMom]      = 0.0;           // v = 0
                            prim[ZMom]      = 0.0;           // w = 0
                            prim[Energy]    = 1.0;           // p = 1.0
                            prim[Alpha1]    = epss;          // α₁ ≈ 0
                        } else if (x > x_shock + epss) {
                            // 区域 2: 已受冲击空气（激波后，右侧）
                            prim[AlphaRho1] = 1.411;         // α₁ρ₁ = ρ_post_air
                            prim[AlphaRho2] = epss;          // α₂ρ₂ ≈ 0
                            prim[XMom]      = -0.39;         // u = -0.39（向左）
                            prim[YMom]      = 0.0;           // v = 0
                            prim[ZMom]      = 0.0;           // w = 0
                            prim[Energy]    = 1.628;         // p = 1.628
                            prim[Alpha1]    = 1.0 - epss;    // α₁ ≈ 1
                        } else {
                            // 区域 3: 未受冲击空气（激波前，介于界面与激波面之间）
                            prim[AlphaRho1] = 1.0;           // α₁ρ₁ = ρ_air
                            prim[AlphaRho2] = epss;          // α₂ρ₂ ≈ 0
                            prim[XMom]      = 0.0;           // u = 0
                            prim[YMom]      = 0.0;           // v = 0
                            prim[ZMom]      = 0.0;           // w = 0
                            prim[Energy]    = 1.0;           // p = 1.0
                            prim[Alpha1]    = 1.0 - epss;    // α₁ ≈ 1
                        }
                        break;
                    }
                }

                // 将原始变量转换为守恒变量
                VarArray cons = prim_to_cons(prim, pp);

                // 使用 Array4 将守恒变量存储到状态数组中
                state_arr(i, j, k, AlphaRho1) = cons[AlphaRho1];
                state_arr(i, j, k, AlphaRho2) = cons[AlphaRho2];
                state_arr(i, j, k, XMom) = cons[XMom];
                state_arr(i, j, k, YMom) = cons[YMom];
                state_arr(i, j, k, ZMom) = cons[ZMom];
                state_arr(i, j, k, Energy) = cons[Energy];
                state_arr(i, j, k, Alpha1) = cons[Alpha1];
            }
        }
    }
}
#endif

