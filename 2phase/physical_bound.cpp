/**
 * @file
 * @brief 两相状态物理约束辅助函数的实现。
 * @ingroup two_phase_model
 */

#include "physical_bound.H"
#include "equation.H"

using namespace amrex;

/**
 * 计算每个维度的alpha值
 * @param dx 网格间距
 * @param primitiveVars 原始变量
 * @param pp 问题参数
 * @return 每个维度的alpha值
 * @功能：用于时间步长计算，根据速度和声速计算每个维度的权重
 */
AMREX_GPU_HOST_DEVICE
GpuArray<Real, BL_SPACEDIM> computeAlpha(const GpuArray<Real, BL_SPACEDIM>& dx, const VarArray& primitiveVars, const ProbParm& pp)
{
    // 提取速度分量
    GpuArray<Real, BL_SPACEDIM> velocity;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        velocity[i] = primitiveVars[XMom + i]; // 速度分量
    }

    // 计算声速
    Real c = prims_to_eos<SoundSpeed>(primitiveVars, pp);

    // 计算每个维度的 tau
    GpuArray<Real, BL_SPACEDIM> tau;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        tau[i] = (std::abs(velocity[i]) + c) / dx[i];
    }

    // 计算 alpha
    GpuArray<Real, BL_SPACEDIM> alpha;
    Real tau_sum = 0.0;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        tau_sum += tau[i];
    }
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        alpha[i] = tau[i] / tau_sum;
    }

    return alpha;
}

/**
 * 正量函数运算符重载
 * @param U 守恒变量
 * @param pp 问题参数
 * @return 需要保持为正的物理量数组
 * @功能：返回需要保持为正的物理量，用于正限制器
 */
AMREX_GPU_HOST_DEVICE
amrex::GpuArray<amrex::Real, 5>
PositiveQuantityFunctions::operator()(const VarArray& U, const ProbParm& pp) const
{
    amrex::GpuArray<amrex::Real, 5> result; // 输出需要保持为正的物理量

    // 确保 U 的索引范围有效
    AMREX_ASSERT(U.size() > 3 + BL_SPACEDIM);

    result[0] = U[AlphaRho1]; // 第一相密度
    result[1] = U[AlphaRho2]; // 第二相密度
    result[2] = U[Alpha1]; // 第一相体积分数
    result[3] = 1 - U[Alpha1]; // 第二相体积分数
    result[4] = cons_to_eos<SoundSpeedSqrDensity>(U, pp); // 声速平方密度

    return result;
}