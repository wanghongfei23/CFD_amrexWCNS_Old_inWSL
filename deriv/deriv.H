// 导数计算头文件
// 包含一阶导数项的模板类和类型定义

#pragma once
#include "deriv_usr.H"
#include "Prob_Parm.H"

// 基类模板：不再使用 CRTP，改为组合
// 用于计算一阶导数项
template <typename Derived>
class FirstOrderDerivativeTerm {
public:
    // 直接使用 Derived 的 num_terms
    static constexpr int num_terms = Derived::num_terms;

    // 构造函数：初始化 Derived 对象
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE explicit FirstOrderDerivativeTerm()
        : derived_obj()
        , // 默认构造 Derived 对象
        target_term_indices(Derived::getTermIndices())
    {
    }

    // 计算系数（通过 Derived 对象调用）
    // 输入：原始变量prim，问题参数pp，法向量norm
    // 输出：系数数组
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::GpuArray<amrex::Real, num_terms>
        computeCoefficient(
            const VarArray& prim,
            const ProbParm& pp,
            const amrex::GpuArray<amrex::Real, 3> norm) const
    {
        return derived_obj.computeCoefficientImpl(prim, pp, norm);
    }

    // 计算半节点 Flux（通过 Derived 对象调用）
    // 输入：左侧原始变量left_prim，右侧原始变量right_prim，问题参数pp，法向量norm，网格尺寸dx
    // 输出：通量数组
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::GpuArray<amrex::Real, num_terms>
        computeFlux(
            const VarArray& left_prim,
            const VarArray& right_prim,
            const ProbParm& pp,
            const amrex::GpuArray<amrex::Real, 3> norm, amrex::GpuArray<amrex::Real, BL_SPACEDIM> dx) const
    {
        return derived_obj.computeFluxImpl(left_prim, right_prim, pp, norm, dx);
    }

    // 计算节点 Flux（通过 Derived 对象调用）
    // 输入：原始变量prim，问题参数pp，法向量norm
    // 输出：节点通量数组
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::GpuArray<amrex::Real, num_terms>
        computeNodeFlux(
            const VarArray& prim,
            const ProbParm& pp,
            const amrex::GpuArray<amrex::Real, 3> norm) const
    {
        return derived_obj.computeNodeFluxImpl(prim, pp, norm);
    }

    // 获取目标项索引
    // 输出：目标项索引数组
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
        amrex::GpuArray<int, num_terms>
        getTargetTermIndices() const
    {
        return target_term_indices;
    }

private:
    Derived derived_obj; // 组合代替继承
    amrex::GpuArray<int, num_terms> target_term_indices;
};

// 对流项类型定义
using AdvectionTerm = FirstOrderDerivativeTerm<FirstOrderDerivativeTermUsr>;
// using AdvectionTerm = FirstOrderDerivativeTerm<FirstOrderDerivativeTermUsr2>;//植入压缩项