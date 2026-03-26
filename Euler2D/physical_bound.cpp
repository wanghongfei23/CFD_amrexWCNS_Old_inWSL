// 物理边界实现文件
// 包含计算alpha值的函数，用于时间步进中的方向权重

#include "physical_bound.H"

using namespace amrex;

// 计算alpha值函数
// 输入：网格间距dx，原始变量primitiveVars，问题参数pp
// 输出：每个维度的alpha值数组
AMREX_GPU_HOST_DEVICE
GpuArray<Real, BL_SPACEDIM> computeAlpha(const GpuArray<Real, BL_SPACEDIM>& dx, const VarArray& primitiveVars, const ProbParm& pp)
{
    // 假设 primitiveVars 存储了密度 (rho), 速度 (u, v, w), 压力 (p)
    Real rho = primitiveVars[0];
    GpuArray<Real, BL_SPACEDIM> velocity;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        velocity[i] = primitiveVars[i + 1]; // 速度分量
    }
    Real p = primitiveVars[BL_SPACEDIM + 1];

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