/**
 * @file
 * @brief 二维单相 Euler 参考实现中非守恒导数项的空策略。
 * @ingroup legacy_euler
 */

#pragma once
// #include "deriv.H"
#include "equation.H"
/**
 * @brief 不包含附加非守恒导数项的占位策略类。
 */
class FirstOrderDerivativeTermUsr {
public:
    // 编译期确定的元信息
    static constexpr int num_terms = 0;
    static constexpr amrex::GpuArray<int, num_terms> term_indices {};

    // 无参数构造函数
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE explicit FirstOrderDerivativeTermUsr()
    {
    }

    // 自动获取目标项索引（返回空数组）
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE static amrex::GpuArray<int, num_terms> getTermIndices()
    {
        return term_indices;
    }

    // 空实现（返回空数组）
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