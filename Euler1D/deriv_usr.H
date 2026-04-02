/**
 * @file
 * @brief 一维单相 Euler 参考实现中非守恒导数项的空策略。
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