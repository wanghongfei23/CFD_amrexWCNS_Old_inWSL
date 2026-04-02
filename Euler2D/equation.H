/**
 * @file
 * @brief 二维单相 Euler 参考实现的状态方程与初值定义。
 * @ingroup legacy_euler
 */

#pragma once

#include "Prob_Parm.H"
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
    constexpr Real gamma = 1.4; // 比热比
    AMREX_ASSERT(BL_SPACEDIM == 2); // 确保是二维问题
    const int* lo = box.loVect(); // 网格下界
    const int* hi = box.hiVect(); // 网格上界
    Array4<Real> const& state_array = state.array(); // 状态数组

#ifdef AMREX_USE_OMP
#pragma omp parallel for if (Gpu::notInLaunchRegion())
#endif
    for (int i = lo[0]; i <= hi[0]; ++i) {
        for (int j = lo[1]; j <= hi[1]; ++j) {
            // 计算当前网格中心的物理坐标
            Real x = problo[0] + (i + 0.5) * dx[0];
            Real y = problo[1] + (j + 0.5) * dx[1];

            // 根据初始条件设置变量
            Real rho, u, v, p;
            // if (x > 0.8 && y > 0.8) {
            //     rho = 1.5;
            //     u = 0.0;
            //     v = 0.0;
            //     p = 1.5;
            // } else if (x < 0.8 && y > 0.8) {
            //     rho = 0.5323;
            //     u = 1.206;
            //     v = 0.0;
            //     p = 0.3;
            // } else if (x < 0.8 && y < 0.8) {
            //     rho = 0.138;
            //     u = 1.206;
            //     v = 1.206;
            //     p = 0.029;
            // } else {
            //     rho = 0.5323;
            //     u = 0.0;
            //     v = 1.206;
            //     p = 0.3;
            // }
            // riemann2
            // if (x > 0.5 && y > 0.5) {
            //     rho = 1.0;
            //     u = 0.75;
            //     v = -0.5;
            //     p = 1.0;
            // } else if (x < 0.5 && y > 0.5) {
            //     rho = 2.0;
            //     u = 0.75;
            //     v = 0.5;
            //     p = 1.0;
            // } else if (x < 0.5 && y < 0.5) {
            //     rho = 1.0;
            //     u = -0.75;
            //     v = 0.5;
            //     p = 1.0;
            // } else {
            //     rho = 3.0;
            //     u = -0.75;
            //     v = -0.5;
            //     p = 1.0;
            // }
            // riemann3
            if (x > 0.5 && y > 0.5) {
                rho = 0.5313;
                u = 0.0;
                v = 0.0;
                p = 0.4;
            } else if (x < 0.5 && y > 0.5) {
                rho = 1.0;
                u = 0.7276;
                v = 0.0;
                p = 1.0;
            } else if (x < 0.5 && y < 0.5) {
                rho = 0.8;
                u = 0.0;
                v = 0.0;
                p = 1.0;
            } else {
                rho = 1.0;
                u = 0.0;
                v = 0.7276;
                p = 1.0;
            }

            // double Mach
            if (y >= std::sqrt(3.0) * (x - 1.0 / 6.0)) {
                rho = 8.0;
                u = 8.25 * cos(M_PI / 6);
                v = -8.25 * sin(M_PI / 6);
                p = 116.5;
            } else {
                rho = 1.4;
                u = 0;
                v = 0;
                p = 1.0;
            }

            // 根据欧拉方程计算守恒变量
            state_array(i, j, 0, 0) = rho; // 密度
            state_array(i, j, 0, 1) = rho * u; // x 方向动量
            state_array(i, j, 0, 2) = rho * v; // y 方向动量
            state_array(i, j, 0, 3) = p / (gamma - 1) + 0.5 * rho * (u * u + v * v); // 总能量
        }
    }
}

// 守恒变量转换为原始变量
// 输入：守恒变量cons，问题参数pp
// 输出：原始变量
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
VarArray cons_to_prim(VarArray cons, ProbParm const& pp)
{
    amrex::Real gamma = pp.gamma; // 气体绝热指数
    amrex::Real r, ru, rv, rE;
    r = cons[0]; // 密度
    ru = cons[1]; // x方向动量
    rv = cons[2]; // y方向动量
    rE = cons[3]; // 总能量
    amrex::Real u = ru / r; // x方向速度
    amrex::Real v = rv / r; // y方向速度
    amrex::Real E = rE / r; // 单位质量总能量
    amrex::Real q2 = (u * u + v * v) / 2; // 单位质量动能
    amrex::Real e = -q2 + E; // 单位质量内能
    amrex::Real RT = (gamma - 1) * e; // 气体常数乘温度
    amrex::Real p = r * RT; // 压力
    return { r, u, v, p };
}

// 原始变量转换为守恒变量
// 输入：原始变量prim，问题参数pp
// 输出：守恒变量
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
VarArray prim_to_cons(VarArray prim, ProbParm const& pp)
{
    // 从原始变量中提取值
    amrex::Real gamma = pp.gamma; // 气体绝热指数
    amrex::Real rho = prim[0]; // 密度
    amrex::Real u = prim[1]; // x 方向速度
    amrex::Real v = prim[2]; // y 方向速度
    amrex::Real p = prim[3]; // 压力

    // 计算动能
    amrex::Real ekin = 0.5 * (u * u + v * v); // 单位质量动能

    // 计算守恒变量
    amrex::Real rho_u = rho * u; // x 方向动量
    amrex::Real rho_v = rho * v; // y 方向动量
    amrex::Real E = p / (gamma - 1) + rho * ekin; // 总能量

    // 返回守恒变量
    return { rho, rho_u, rho_v, E };
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

// 原始变量转换为通量
// 输入：原始变量prim，问题参数pp，法向量norm
// 输出：通量
 template <int idim = 0>
AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        VarArray
        prim_to_flux(VarArray prim, ProbParm pp, amrex::GpuArray<amrex::Real, 3> norm = { 1, 0, 0 })
{
    // finish
    amrex::Real r = prim[0], u = prim[1], v = prim[2], p = prim[3]; // 密度、速度、压力
    amrex::Real gamma = pp.gamma; // 气体绝热指数

    amrex::Real Vn = (u * norm[0] + v * norm[1]); // 法向速度
    // amrex::Real H=GAMMA/(GAMMA-1)*p/r+(u*u+v*v)/2;
    VarArray res = {
        r * Vn, // 密度通量
        r * u * Vn + p * norm[0], // x方向动量通量
        r * v * Vn + p * norm[1], // y方向动量通量
        ((u * u + v * v) / 2 * r + gamma / (gamma - 1) * p) * Vn // 能量通量
    };
    // for (int i = 0; i < res.size(); ++i) {
    //     std::cout << "res[" << i << "] = " << res[i] << std::endl;
    // }
    return res;
}

// HLLC通量计算函数（二维）
// 输入：左侧密度rl，右侧密度rr，左侧x速度ul，右侧x速度ur，左侧y速度vl，右侧y速度vr，左侧压力pl，右侧压力pr，法向量norm，气体绝热指数gamma
// 输出：通量数组
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
VarArray HLLCFlux2D2(amrex::Real rl, amrex::Real rr,
    amrex::Real ul, amrex::Real ur,
    amrex::Real vl, amrex::Real vr,
    amrex::Real pl, amrex::Real pr,
    amrex::GpuArray<amrex::Real, 3> norm, amrex::Real gamma = 1.4)
{
    using namespace amrex;
    enum {
        L, R // 左侧和右侧
    };
    VarArray res; // 结果通量
    // if(pl<0 || pr<0 || rr<0 || rl<0)
    // {
    //     std::cout<<"fluxScheme error: positivity break\n";
    // }

    // 计算左右侧的总焓
    amrex::GpuArray<amrex::Real, 2> H;
    H[L] = (ul * ul + vl * vl) / 2 + pl / rl * gamma / (gamma - 1);
    H[R] = (ur * ur + vr * vr) / 2 + pr / rr * gamma / (gamma - 1);

    // 计算法向速度
    amrex::GpuArray<amrex::Real, 2> Vn = { (norm[1] > norm[0]) ? vl : ul, (norm[1] > norm[0]) ? vr : ur };
    amrex::Real cl = sqrt(pl / rl * gamma); // 左侧声速
    amrex::Real cr = sqrt(pr / rr * gamma); // 右侧声速

    // 计算波速
    amrex::Real rBar = (rl + rr) / 2, cBar = (cl + cr) / 2;
    amrex::Real ps = std::max(0.0, (pl + pr) / 2 - rBar * cBar * (Vn[R] - Vn[L]) / 2);
    amrex::Real ql = (ps <= pl) ? 1.0 : sqrt(1.0 + (gamma + 1.0) / (2.0 * gamma) * (ps / pl - 1.0));
    amrex::Real qr = (ps <= pr) ? 1.0 : sqrt(1.0 + (gamma + 1.0) / (2.0 * gamma) * (ps / pr - 1.0));
    amrex::Real SL = Vn[L] - cl * ql, SR = Vn[R] + cr * qr; // 左侧和右侧波速

    // 计算中间波速
    amrex::Real Sstar = ((pr - pl) + (rl * Vn[L] * (SL - Vn[L]) - rr * Vn[R] * (SR - Vn[R]))) / (rl * (SL - Vn[L]) - rr * (SR - Vn[R]));

    // 根据波速位置计算通量
    if (SL > 0) {
        // 全部在左侧
        res[0] = rl * Vn[L]; // 密度通量
        res[1] = rl * ul * Vn[L] + pl * norm[0]; // x方向动量通量
        res[2] = rl * vl * Vn[L] + pl * norm[1]; // y方向动量通量
        res[3] = rl * H[L] * Vn[L]; // 能量通量
    } else if (Sstar > 0) {
        amrex::Real U[4] = { rl,
            rl * ul,
            rl * vl,
            rl * H[L] - pl };
        amrex::Real F[4] = { rl * Vn[L], rl * ul * Vn[L] + pl * norm[0], rl * vl * Vn[L] + pl * norm[1], rl * H[L] * Vn[L] };
        amrex::Real Usf = (SL - Vn[L]) / (SL - Sstar);
        amrex::Real Ustar[4] = { Usf * rl,
            Usf * (Sstar * norm[0] + ul * norm[1]) * rl,
            Usf * (Sstar * norm[1] + vl * norm[0]) * rl,
            Usf * (U[3] + (Sstar - Vn[L]) * (rl * Sstar + pl / (SL - Vn[L]))) };
        for (int i = 0; i < 4; i++) {
            res[i] = F[i] + SL * (Ustar[i] - U[i]);
        }
    } else if (Sstar == 0) {
        // 在中间
        {
            amrex::Real U[4] = { rl,
                rl * ul,
                rl * vl,
                rl * H[L] - pl };
            amrex::Real F[4] = { rl * Vn[L], rl * ul * Vn[L] + pl * norm[0], rl * vl * Vn[L] + pl * norm[1], rl * H[L] * Vn[L] };
            amrex::Real Usf = (SL - Vn[L]) / (SL - Sstar);
            amrex::Real Ustar[4] = { Usf * rl,
                Usf * (Sstar * norm[0] + ul * norm[1]) * rl,
                Usf * (Sstar * norm[1] + vl * norm[0]) * rl,
                Usf * (U[3] + (Sstar - Vn[L]) * (rl * Sstar + pl / (SL - Vn[L]))) };
            for (int i = 0; i < 4; i++) {
                res[i] = (F[i] + SL * (Ustar[i] - U[i])) / 2;
            }
        }

        {
            amrex::Real U[4] = { rr, rr * ur, rr * vr, rr * H[R] - pr };

            amrex::Real F[4] = { rr * Vn[R], rr * ur * Vn[R] + pr * norm[0], rr * vr * Vn[R] + pr * norm[1], rr * H[R] * Vn[R] };
            amrex::Real Usf = (SR - Vn[R]) / (SR - Sstar);
            amrex::Real Ustar[4] = { Usf * rr,
                Usf * (Sstar * norm[0] + ur * norm[1]) * rr,
                Usf * (Sstar * norm[1] + vr * norm[0]) * rr,
                Usf * (U[3] + (Sstar - Vn[R]) * (rr * Sstar + pr / (SR - Vn[R]))) };
            for (int i = 0; i < 4; i++) {
                res[i] += (F[i] + SR * (Ustar[i] - U[i])) / 2;
            }
        }

    } else if (SR > 0) {
        amrex::Real U[4] = { rr, rr * ur, rr * vr, rr * H[R] - pr };

        amrex::Real F[4] = { rr * Vn[R], rr * ur * Vn[R] + pr * norm[0], rr * vr * Vn[R] + pr * norm[1], rr * H[R] * Vn[R] };
        amrex::Real Usf = (SR - Vn[R]) / (SR - Sstar);
        amrex::Real Ustar[4] = { Usf * rr,
            Usf * (Sstar * norm[0] + ur * norm[1]) * rr,
            Usf * (Sstar * norm[1] + vr * norm[0]) * rr,
            Usf * (U[3] + (Sstar - Vn[R]) * (rr * Sstar + pr / (SR - Vn[R]))) };
        for (int i = 0; i < 4; i++) {
            res[i] = F[i] + SR * (Ustar[i] - U[i]);
        }
    } else {
        // 全部在右侧
        res[0] = rr * Vn[R]; // 密度通量
        res[1] = rr * ur * Vn[R] + pr * norm[0]; // x方向动量通量
        res[2] = rr * vr * Vn[R] + pr * norm[1]; // y方向动量通量
        res[3] = rr * H[R] * Vn[R]; // 能量通量
    }
    for (int i = 0; i < res.size(); ++i) {
        if (isnan(res[i]))
            std::cout << "res[" << i << "] = " << res[i] << std::endl;
    }
    return res;
}

// Riemann求解器
// 输入：网格坐标i,j,k，原始变量数组prim_array，左侧状态stateL，右侧状态stateR，法向量norm，问题参数pp
// 输出：通量
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
VarArray
Riemann_solver_1(int /*i*/, int /*j*/, int /*k*/, const amrex::Array4<const amrex::Real>& /*prim_array*/, VarArray const& stateL, VarArray const& stateR, amrex::GpuArray<amrex::Real, 3> norm, ProbParm const& pp)
{
    return HLLCFlux2D2(stateL[0], stateR[0],
        stateL[1], stateR[1],
        stateL[2], stateR[2],
        stateL[3], stateR[3], norm, pp.gamma);
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
    amrex::Real rhou = state(i, j, k, 1); // x 方向动量
    amrex::Real rhov = state(i, j, k, 2); // y 方向动量
    amrex::Real rhoE = state(i, j, k, 3); // 总能量

    amrex::Real u = rhou / rho; // x 方向速度
    amrex::Real v = rhov / rho; // y 方向速度
    amrex::Real E = rhoE / rho; // 单位质量的总能量
    amrex::Real e = E - 0.5 * (u * u + v * v); // 单位质量的内能

    amrex::Real RT = (gamma - 1) * e; // 气体常数乘温度
    amrex::Real p = rho * RT; // 压力

    // 计算 x 方向和 y 方向的最大特征速度
    amrex::Real lambda_x = (std::sqrt(gamma * RT) + std::abs(u)) / dx[0];
    amrex::Real lambda_y = (std::sqrt(gamma * RT) + std::abs(v)) / dx[1];

    // 输出最大特征速度
    umax(i, j, k) = lambda_x + lambda_y;
}