// // deriv/AdvectionTerm.H
#pragma once
#include "equation.H"

class FirstOrderDerivativeTermUsr {
public:
    // 编译期确定的元信息
    static constexpr int num_terms = 2;
    static constexpr amrex::GpuArray<int, num_terms> term_indices { Alpha1, Alpha1 }; // 硬编码目标方程索引

    // 无参数构造函数
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE explicit FirstOrderDerivativeTermUsr()
    {
    }

    // 自动获取目标项索引（返回空数组）
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE static amrex::GpuArray<int, num_terms> getTermIndices()
    {
        return term_indices;
    }
    // 实现系数计算
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::GpuArray<amrex::Real, num_terms>
        computeCoefficientImpl(
            const VarArray& prim,
            const ProbParm& pp,
            const amrex::GpuArray<amrex::Real, 3> norm) const
    {
        return { -prim[Alpha1], 1. };
    }
    // 计算节点 Flux
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::GpuArray<amrex::Real, num_terms>
        computeNodeFluxImpl(
            const VarArray& state,
            const ProbParm& pp,
            const amrex::GpuArray<amrex::Real, 3> norm) const
    {
        using namespace amrex;
        enum VelComponents {
            Ar1 = 0,
            Ar2,
            AMREX_D_DECL(U, V, W),
            Pressure,
            Alpha1
        };
        // 计算法向速度
        Real Vn = AMREX_D_TERM(norm[0] * state[U], +norm[1] * state[V], +norm[2] * state[W]);

        return { Vn, state[Alpha1] * Vn };
    }

    // 实现基于HLLC近似Riemann求解器的半节点通量计算
    template <bool first_order = false>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::GpuArray<amrex::Real, num_terms>
        computeFluxImpl(
            const VarArray& stateL, // 左侧状态变量数组
            const VarArray& stateR, // 右侧状态变量数组
            const ProbParm& pp, // 物理问题参数
            const amrex::GpuArray<amrex::Real, 3> norm, amrex::GpuArray<amrex::Real, BL_SPACEDIM> dx) const // 界面法向量
    {
        using namespace amrex;

        // 状态变量分量枚举定义
        enum VelComponents {
            Ar1 = 0, // 相1面积分数
            Ar2, // 相2面积分数
            AMREX_D_DECL(U, V, W), // 速度分量 (根据维度展开)
            Pressure, // 压力
            Alpha1 // 相1体积分数
        };
        enum { L,
            R }; // 左右状态标识符

        // 提取左右状态变量
        Real ar1l
            = stateL[Ar1],
            ar2l = stateL[Ar2];
        Real ar1r = stateR[Ar1], ar2r = stateR[Ar2];
        Real rl = ar1l + ar2l; // 左侧总密度
        Real rr = ar1r + ar2r; // 右侧总密度
                               // 提取压力参数
        Real pl = stateL[Pressure], pr = stateR[Pressure];
        Real alpha1l = stateL[Alpha1], alpha1r = stateR[Alpha1];

        //------------------------------------------------------
        Real kappa = 0.0;
        // 开关函数：仅在界面区激活
        Real alpha_avg = 0.5 * (alpha1l + alpha1r);
        Real S_switch = alpha_avg * (1.0 - alpha_avg); // ∈ [0,1]

        // --- 6. 计算压缩速度（受 CFL 限制）---
        Real u_comp = kappa;

        // 若你有 dt 和 dx，可进一步限幅（推荐）：
        // Real u_comp_max = C_comp * dx_norm / dt;
        // u_comp = amrex::min(u_comp, u_comp_max);

        // 确定压缩方向：从 alpha 小指向大（或用法向）
        Real n_dir = (alpha1r > alpha1l) ? 1.0 : -1.0;
        // 或更稳健：n_dir = (Sstar >= 0) ? 1.0 : -1.0; // 但推荐用 alpha 梯度
        //------------------------------------------------------

        // 提取速度分量（自动处理不同维度）
        Real AMREX_D_DECL(ul = stateL[U],
            vl = stateL[V],
            wl = stateL[W]);

        Real AMREX_D_DECL(uR = stateR[U],
            vR = stateR[V],
            wR = stateR[W]);

        // 计算法向速度分量
        GpuArray<Real, 2> Vn = {
            AMREX_D_TERM(ul * norm[0] + n_dir * u_comp, +vl * norm[1], +wl * norm[2]), // 左侧法向速度
            AMREX_D_TERM(uR * norm[0] + n_dir * u_comp, +vR * norm[1], +wR * norm[2]) // 右侧法向速度
        };

        // 计算声速
        Real cl = eos_prim<SoundSpeed>(stateL, pp); // 左侧声速
        Real cr = eos_prim<SoundSpeed>(stateR, pp); // 右侧声速

        // if (std::isnan(cl) || cl < 0) {
        //     std::string stateL_str = "StateL: [";
        //     for (int i = 0; i < NVar; ++i) {
        //         stateL_str += std::to_string(stateL[i]);
        //         if (i < NVar - 1)
        //             stateL_str += ", ";
        //     }
        //     stateL_str += "]";
        //     std::string error_msg = "Error: cl is invalid (cl = " + std::to_string(cl) + "). " + stateL_str;
        //     amrex::Abort(error_msg.c_str());
        // }
        // if (std::isnan(cr) || cr < 0) {
        //     std::string stateR_str = "StateR: [";
        //     for (int i = 0; i < NVar; ++i) {
        //         stateR_str += std::to_string(stateR[i]);
        //         if (i < NVar - 1)
        //             stateR_str += ", ";
        //     }
        //     stateR_str += "]";
        //     std::string error_msg = "Error: cr is invalid (cr = " + std::to_string(cr) + "). " + stateR_str;
        //     amrex::Abort(error_msg.c_str());
        // }

        // 计算平均密度的平方根
        Real sqrt_rl = std::sqrt(rl);
        Real sqrt_rr = std::sqrt(rr);
        // if (std::isnan(sqrt_rl)) {
        //     amrex::Abort("Error: sqrt_rl is NaN. Check the value of rl.");
        // }
        // if (std::isnan(sqrt_rr)) {
        //     amrex::Abort("Error: sqrt_rr is NaN. Check the value of rr.");
        // }

        // 计算加权平均速度
        Real VnBar = (Vn[L] * sqrt_rl + Vn[R] * sqrt_rr) / (sqrt_rl + sqrt_rr);

        // // 计算加权平均速度
        // Real VnBar = (Vn[L] * sqrt_rl + Vn[R] * sqrt_rr) / (sqrt_rl + sqrt_rr);

        // 计算加权平均声速平方
        Real cBar = sqrt((cl * cl * sqrt_rl + cr * cr * sqrt_rr) / (sqrt_rl + sqrt_rr)
            + 0.5 * (sqrt_rl * sqrt_rr) / ((sqrt_rl + sqrt_rr) * (sqrt_rl + sqrt_rr)) * (Vn[R] - Vn[L]) * (Vn[R] - Vn[L]));

        // 计算波速估计
        Real Ma_local = std::max(std::abs(Vn[L] / cl), std::abs(Vn[R] / cr)), Ma_limit = 0.1;
        Real phi = sin(std::min(1., Ma_local / Ma_limit) * M_PI_2);
        Real SL = std::min(Vn[L] - cl, VnBar - cBar); // 左行波速
        Real SR = std::max(Vn[R] + cr, VnBar + cBar); // 右行波速

        // 计算星区速度（接触间断速度）
        Real numerator = (pr - pl) + (rl * Vn[L] * (SL - Vn[L]) - rr * Vn[R] * (SR - Vn[R]));
        Real denominator = rl * (SL - Vn[L]) - rr * (SR - Vn[R]);
        Real Sstar = numerator / denominator; // 接触间断速度
        Real alpha1_half = Sstar > 0 ? alpha1l : alpha1r;
        Real u_star = Sstar;
        // if (SL >= 0) {
        //     u_star = Vn[L];
        // } else if (Sstar > 0) {
        //     Real Usf = (SL - Vn[L]) / (SL - Sstar);
        //     u_star = Vn[L] + SL * (Usf - 1.0);
        // } else if (SR > 0) {
        //     Real Usf = (SR - Vn[R]) / (SR - Sstar);
        //     u_star = Vn[R] + SR * (Usf - 1.0);
        // } else {
        //     u_star = Vn[R];
        // }

        return { u_star, u_star * alpha1_half };

        // if (SL >= 0) {
        //     return {
        //         Vn[L], Vn[L] * alpha1l
        //     };
        // } else if (SR <= 0) {
        //     return {
        //         Vn[R], Vn[R] * alpha1r
        //     };
        // }

        // Real UsfR = (SR - Vn[R]) / (SR - Sstar);
        // Real UsfL = (SL - Vn[L]) / (SL - Sstar);
        // phi = std::max(std::max(phi, (alpha1l + alpha1r) / 2), 1. - (alpha1l + alpha1r) / 2);
        // Real u_contact
        //     = (Vn[L] + Vn[R]) / 2 + (SL * phi * (UsfL - 1.) + std::abs(Sstar) * (UsfL - UsfR) + SR * phi * (UsfR - 1.)) / 2;
        // Real alpha1_half = u_contact > 0 ? alpha1l : alpha1r;

        // return { u_contact, u_contact * alpha1_half };
    }
};