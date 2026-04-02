/**
 * @file
 * @brief 二维单相 Euler 参考实现中的特征矩阵构造与变量变换。
 * @ingroup legacy_euler
 */

#include "eigen_system.H"
#include "equation.H"
#include "linear_algebra.H"

using namespace amrex;

/**
 * 构造函数，初始化特征系统
 * @param primL 左侧原始变量
 * @param primR 右侧原始变量
 * @param norm 法向量
 * @param pp 问题参数
 * @功能：计算平均状态并初始化左右特征矩阵
 */
AMREX_GPU_HOST_DEVICE
EigenSystem::EigenSystem(VarArray const& primL, VarArray const& primR, GpuArray<Real, 3> norm, ProbParm const& pp)
{
    // 使用 Roe 平均计算中间状态
    VarArray prim = average(primL, primR, norm, pp);
    enum { X, Y };

    // 提取物理量
    Real gamma = pp.gamma; // 气体绝热指数
    Real r = prim[0]; // 密度
    Real u = prim[1]; // x方向速度
    Real v = prim[2]; // y方向速度
    Real p = prim[3]; // 压力
    Real ek = (u * u + v * v) / 2; // 动能
    Real h = p / r * gamma / (gamma - 1); // 比焓
    Real c = sqrt(gamma * p / r); // 声速
    Real ht = h + ek; // 总焓
    Real Vn = u * norm[X] + v * norm[Y]; // 法向速度

    // 计算左特征矩阵
    leftEig = { (-norm[Y] * u + norm[X] * v),
        norm[Y],
        (-norm[X]),
        0,
        (h - ek),
        u,
        v,
        -1.0,
        (Vn / c + ek / h) / 2,
        (-norm[X] / c - u / h) / 2,
        (-norm[Y] / c - v / h) / 2,
        1.0 / (2.0 * h),
        (-Vn / c + ek / h) / 2,
        (norm[X] / c - u / h) / 2,
        (norm[Y] / c - v / h) / 2,
        1.0 / (2.0 * h) };

    // 计算右特征矩阵
    rightEig = { 0,
        1.0 / h,
        1,
        1,
        norm[Y],
        u / h,
        u - norm[X] * c,
        u + norm[X] * c,
        -norm[X],
        v / h,
        v - norm[Y] * c,
        v + norm[Y] * c,
        norm[Y] * u - norm[X] * v,
        ek / h,
        h + ek - Vn * c,
        h + ek + Vn * c };
}

/**
 * 将原始变量转换为特征变量
 * @param prim 原始变量
 * @param pp 问题参数
 * @return 特征变量charVars
 */
AMREX_GPU_HOST_DEVICE
VarArray
EigenSystem::prim_to_char(VarArray const& prim, ProbParm const& pp) const
{
    // 将原始变量转换为守恒变量
    VarArray cons = prim_to_cons(prim, pp);

    // 使用左特征矩阵将守恒变量转换为特征变量
    VarArray charVars = LA_M_X<1, NVar>(leftEig, cons);

    return charVars;
}

/**
 * 将特征变量转换为原始变量
 * @param charVars 特征变量
 * @param pp 问题参数
 * @return 原始变量prim
 */
AMREX_GPU_HOST_DEVICE
VarArray
EigenSystem::char_to_prim(VarArray const& charVars, ProbParm const& pp) const
{
    // 使用右特征矩阵将特征变量转换为守恒变量
    VarArray cons = LA_M_X<2, NVar>(rightEig, charVars);

    // 将守恒变量转换为原始变量
    VarArray prim = cons_to_prim(cons, pp);

    return prim;
}

/**
 * 将特征变量转换为守恒变量
 * @param charVars 特征变量
 * @param pp 问题参数
 * @return 守恒变量cons
 */
AMREX_GPU_HOST_DEVICE
VarArray
EigenSystem::char_to_cons(VarArray const& charVars, ProbParm const& pp) const
{
    // 使用右特征矩阵将特征变量转换为守恒变量
    VarArray cons = LA_M_X<2, NVar>(rightEig, charVars);

    return cons;
}

/**
 * 计算左右状态的平均
 * @param primL 左侧原始变量
 * @param primR 右侧原始变量
 * @param norm 法向量
 * @param pp 问题参数
 * @return 平均原始变量
 */
AMREX_GPU_HOST_DEVICE
VarArray EigenSystem::average(VarArray primL,
    VarArray primR,
    GpuArray<Real, 3> norm,
    ProbParm const& pp) const
{
    Real gamma = pp.gamma; // 气体绝热指数
    // 提取左右状态的物理量
    Real rl = primL[0], ul = primL[1], vl = primL[2], pl = primL[3]; // 左侧：密度、x速度、y速度、压力
    Real rr = primR[0], ur = primR[1], vr = primR[2], pr = primR[3]; // 右侧：密度、x速度、y速度、压力
    
    enum { LL, RR };
    GpuArray<Real, 2> H; // 总焓数组
    H[LL] = (ul * ul + vl * vl) / 2 + pl / rl * gamma / (gamma - 1); // 左侧总焓
    H[RR] = (ur * ur + vr * vr) / 2 + pr / rr * gamma / (gamma - 1); // 右侧总焓
    
    // 计算加权系数
    Real coef1 = sqrt(rl);
    Real coef2 = sqrt(rr);
    Real divisor = 1.0 / (coef1 + coef2);

    // 计算平均物理量
    Real r = sqrt(rl * rr); // 平均密度
    Real u = (coef1 * ul + coef2 * ur) * divisor; // 平均x速度
    Real v = (coef1 * vl + coef2 * vr) * divisor; // 平均y速度
    Real Vn = norm[0] * u + norm[1] * v; // 平均法向速度
    Real ht = (coef1 * H[LL] + coef2 * H[RR]) * divisor; // 平均总焓
    Real ek = (u * u + v * v) / 2; // 平均动能
    Real h = ht - ek; // 平均比焓
    Real c = sqrt((gamma - 1) * h); // 平均声速
    Real p = r * h * ((gamma - 1) / gamma); // 平均压力
    
    return { r, u, v, p }; // 返回平均原始变量
}