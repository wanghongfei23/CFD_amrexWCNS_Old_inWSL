// 特征系统实现文件（两相流）
// 包含两相流欧拉方程特征值和特征向量的计算

#include "eigen_system.H"
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

    // 提取物理量
    Real ar1 = prim[0]; // 第一相密度
    Real ar2 = prim[1]; // 第二相密度
    Real c = prim[2]; // 声速
    Real r = prim[3]; // 总密度
    Real alpha1 = prim[4]; // 第一相体积分数
    Real nx = norm[0], ny = norm[1], nz = norm[2]; // 法向量分量

    Real factor = nx + ny + nz;

    // 初始化左特征矩阵
#if BL_SPACEDIM == 1
    Real ar2_over_c2r = ar2 / (c * c * r);
    Real ar1_over_c2r = ar1 / (c * c * r);
    Real half_over_cr = 1.0 / (2 * c * r);
    leftEig = {
        0, 0, 0, 0, 1.,
        0, 1., 0, -ar2_over_c2r, 0,
        1., 0, 0, -ar1_over_c2r, 0,
        0, 0, -1. / 2.0, half_over_cr, 0,
        0, 0, 1. / 2.0, half_over_cr, 0
    };
#elif BL_SPACEDIM == 2
    Real ar2_over_c2r = ar2 / (c * c * r);
    Real ar1_over_c2r = ar1 / (c * c * r);
    Real over_cr = 1.0 / (c * r);
    leftEig = {
        0, 0, 0, 0, 0, 1,
        0, 0, -ny, nx, 0, 0,
        0, 1., 0, 0, -ar2_over_c2r, 0,
        1.0, 0, 0, 0, -ar1_over_c2r, 0,
        0, 0, -nx, -ny, over_cr, 0,
        0, 0, nx, ny, over_cr, 0
    };
#elif BL_SPACEDIM == 3
    // 下面三行是ai给注释时候自己加的，有待商榷
    Real ar2_over_c2r = ar2 / (c * c * r);
    Real ar1_over_c2r = ar1 / (c * c * r);
    Real half_over_cr = 1.0 / (2 * c * r);
    leftEig = {
        0, 0, 0, 0, 0, 0, 1,
        0, 0, nx * (nz - ny) / factor, (-ny * ny + ny * nz + 1) / factor, (-ny * nz + nz * nz - 1) / factor, 0, 0,
        0, 0, (nx * nz + nz * nz + ny * ny) / factor, ny * (nz - nx) / factor, (-nx * nz + nz * nz - 1) / factor, 0, 0,
        0, 1, 0, 0, 0, -ar2_over_c2r, 0,
        1, 0, 0, 0, 0, -ar1_over_c2r, 0,
        0, 0, -nx / 2.0, -ny / 2.0, -nz / 2.0, half_over_cr, 0,
        0, 0, nx / 2.0, ny / 2.0, nz / 2.0, half_over_cr, 0
    };
#endif

    // 初始化右特征矩阵
#if BL_SPACEDIM == 1
    rightEig = {
        0, 0, 1., ar1 / c, ar1 / c,
        0, 1., 0, ar2 / c, ar2 / c,
        0, 0, 0, -1., 1.,
        0, 0, 0, c * r, c * r,
        1., 0, 0, 0, 0
    };
#elif BL_SPACEDIM == 2
    rightEig = {
        0, 0, 0, 1, ar1 / c / 2, ar1 / c / 2,
        0, 0, 1, 0, ar2 / c / 2, ar2 / c / 2,
        0, -ny, 0, 0, -nx / 2, nx / 2,
        0, nx, 0, 0, -ny / 2, ny / 2,
        0, 0, 0, 0, c * r / 2, c * r / 2,
        1, 0, 0, 0, 0, 0
    };
#elif BL_SPACEDIM == 3
    rightEig = {
        0, 0, 0, 0, 1, ar1 / c, ar1 / c,
        0, 0, 0, 1, 0, ar2 / c, ar2 / c,
        0, -ny, ny + nz, 0, 0, -nx, nx,
        0, nx + nz, -nx, 0, 0, -ny, ny,
        0, -ny, -nx, 0, 0, -nz, nz,
        0, 0, 0, 0, 0, c * r, c * r,
        1, 0, 0, 0, 0, 0, 0
    };
#endif
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
    // 提取左值
    Real ar1L = primL[0], ar2L = primL[1], cL = eos_prim<SoundSpeed>(primL, pp), rL = ar1L + ar2L;

    // 提取右值
    Real ar1R = primR[0], ar2R = primR[1], cR = eos_prim<SoundSpeed>(primR, pp), rR = ar1R + ar2R;

    // 计算算数平均
    Real ar1_avg = 0.5 * (ar1L + ar1R);
    Real ar2_avg = 0.5 * (ar2L + ar2R);
    Real c_avg = 0.5 * (cL + cR);
    Real r_avg = 0.5 * (rL + rR);
    Real alpha1_avg = 0.5 * (primL[Alpha1] + primR[Alpha1]);

    // 返回平均值
    return { ar1_avg, ar2_avg, c_avg, r_avg, alpha1_avg };
}

/**
 * 原始变量转为特征变量
 * @param prim 原始变量
 * @param pp 问题参数
 * @return 特征变量charVars
 */
AMREX_GPU_HOST_DEVICE
VarArray EigenSystem::prim_to_char(VarArray const& prim, ProbParm const& pp) const
{
    VarArray charVars = LA_M_X<0, NVar>(leftEig, prim);
    return charVars;
}

/**
 * 特征变量转为原始变量
 * @param charVars 特征变量
 * @param pp 问题参数
 * @return 原始变量prim
 */
AMREX_GPU_HOST_DEVICE
VarArray EigenSystem::char_to_prim(VarArray const& charVars, ProbParm const& pp) const
{
    VarArray prim = LA_M_X<0, NVar>(rightEig, charVars);
    return prim;
}

/**
 * 特征变量转为守恒变量
 * @param charVars 特征变量
 * @param pp 问题参数
 * @return 守恒变量cons
 */
AMREX_GPU_HOST_DEVICE
VarArray EigenSystem::char_to_cons(VarArray const& charVars, ProbParm const& pp) const
{
    VarArray prim = LA_M_X<0, NVar>(rightEig, charVars);
    VarArray cons = prim_to_cons(prim, pp);
    return cons;
}
