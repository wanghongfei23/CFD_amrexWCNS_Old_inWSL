// 物理边界实现文件（一维）
// 包含物理边界相关的函数

#include "physical_bound.H"

using namespace amrex;

/**
 * 计算Alpha系数
 * @param dx 网格间距
 * @param primitiveVars 原始变量
 * @param pp 问题参数
 * @return Alpha系数数组
 * @功能：计算每个维度的alpha系数，用于时间步长计算
 */
AMREX_GPU_HOST_DEVICE
GpuArray<Real, BL_SPACEDIM> computeAlpha(const GpuArray<Real, BL_SPACEDIM>& dx, const VarArray& primitiveVars, const ProbParm& pp)
{
    // 提取原始变量
    Real rho = primitiveVars[0]; // 密度
    GpuArray<Real, BL_SPACEDIM> velocity; // 速度分量
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        velocity[i] = primitiveVars[i + 1];
    }
    Real p = primitiveVars[BL_SPACEDIM + 1]; // 压力

    // 计算声速
    Real c = std::sqrt(pp.gamma * p / rho);

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