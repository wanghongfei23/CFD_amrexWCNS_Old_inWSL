// 导数用户实现头文件
// 包含一阶导数项的用户实现

#pragma once
// #include "deriv.H"
#include "equation.H"

// 空AdvectionTerm实现
class FirstOrderDerivativeTermUsr {
public:
    // 编译期确定的元信息
    static constexpr int num_terms = 0; // 项数为0
    static constexpr amrex::GpuArray<int, num_terms> term_indices {}; // 项索引数组

    // 无参数构造函数
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE explicit FirstOrderDerivativeTermUsr()
    {
    }

    // 自动获取目标项索引（返回空数组）
    // 输出：目标项索引数组
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE static amrex::GpuArray<int, num_terms> getTermIndices()
    {
        return term_indices;
    }

    // 空实现（返回空数组）
    // 输入：原始变量，问题参数，法向量
    // 输出：系数数组
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::GpuArray<amrex::Real, num_terms>
        computeCoefficientImpl(
            const VarArray&,
            const ProbParm&,
            const amrex::GpuArray<amrex::Real, 3>) const
    {
        return {};
    }

    // 空实现（返回空数组）
    // 输入：原始变量，问题参数，法向量
    // 输出：节点通量数组
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::GpuArray<amrex::Real, num_terms>
        computeNodeFluxImpl(
            const VarArray&,
            const ProbParm&,
            const amrex::GpuArray<amrex::Real, 3>) const
    {
        return {};
    }

    // 空实现（返回空数组）
    // 输入：左侧原始变量，右侧原始变量，问题参数，法向量
    // 输出：通量数组
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::GpuArray<amrex::Real, num_terms>
        computeFluxImpl(
            const VarArray&,
            const VarArray&,
            const ProbParm&,
            const amrex::GpuArray<amrex::Real, 3>) const
    {
        return {};
    }
};