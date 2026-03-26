// 特征系统实现文件（一维）
// 包含一维欧拉方程特征值和特征向量的计算

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
    // 计算左右状态的平均
    VarArray prim = average(primL, primR, norm, pp);

    // 提取平均状态的物理量
    Real gamma = pp.gamma; // 气体绝热指数
    Real r = prim[0]; // 密度
    Real u = prim[1]; // 速度
    Real p = prim[2]; // 压力
    Real ek = (u * u) / 2; // 动能
    Real h = p / r * gamma / (gamma - 1); // 比焓
    Real c = sqrt(gamma * p / r); // 声速
    Real ht = h + ek; // 总焓

    // 计算左右特征矩阵
    Real gamma_1 = gamma - 1;
    // 左特征矩阵
    leftEig = { ht + c * (u - c) / gamma_1, -u - c / gamma_1, 1.0,
        -2.0 * ht + 4.0 * c * c / gamma_1, 2.0 * u, -2.0,
        ht - c * (u + c) / gamma_1, -u + c / gamma_1, 1.0 };
    Real factorEig = 0.5 * gamma_1 / (c * c); // 归一化因子

    // 归一化左特征矩阵
    for (int ii = 0; ii < 9; ii++)
        leftEig[ii] *= factorEig;

    // 右特征矩阵
    rightEig = { 1.0, 1.0, 1.0, u - c, u,
        u + c, ht - u * c, 0.5 * u * u, ht + u * c };
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
    Real rl = primL[0], ul = primL[1], pl = primL[2]; // 左侧：密度、速度、压力
    Real rr = primR[0], ur = primR[1], pr = primR[2]; // 右侧：密度、速度、压力
    
    enum { LL, RR };
    GpuArray<Real, 2> H; // 总焓数组
    H[LL] = (ul * ul) / 2 + pl / rl * gamma / (gamma - 1); // 左侧总焓
    H[RR] = (ur * ur) / 2 + pr / rr * gamma / (gamma - 1); // 右侧总焓
    
    // 计算加权系数
    Real coef1 = sqrt(rl);
    Real coef2 = sqrt(rr);
    Real divisor = 1.0 / (sqrt(rl) + sqrt(rr));

    // 计算平均物理量
    Real r = sqrt(rl * rr); // 平均密度
    Real u = (coef1 * ul + coef2 * ur) * divisor; // 平均速度
    Real ht = (coef1 * H[LL] + coef2 * H[RR]) * divisor; // 平均总焓
    Real ek = (u * u) / 2; // 平均动能
    Real h = ht - ek; // 平均比焓
    Real c = sqrt((gamma - 1) * h); // 平均声速
    Real p = r * h * ((gamma - 1) / gamma); // 平均压力
    
    return { r, u, p }; // 返回平均原始变量
}