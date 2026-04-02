/**
 * @file
 * @brief 两相模型的特征矩阵与特征变量变换接口。
 * @ingroup reconstruction
 */

#pragma once

#include "Prob_Parm.H"
#include "equation.H"
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>

/**
 * @brief 在 Roe 平均状态上构建左右特征矩阵。
 */
class EigenSystem {

public:
    /**
     * @brief 使用界面左右原始变量构造特征系统。
     * @param primL 左状态原始变量。
     * @param primR 右状态原始变量。
     * @param norm 面法向量。
     * @param pp 两相物性参数。
     */
    explicit EigenSystem(VarArray const& primL,
        VarArray const& primR,
        amrex::GpuArray<amrex::Real, 3> norm,
        ProbParm const& pp);

    // 禁用默认构造函数
    EigenSystem() = delete;

    // 禁用拷贝构造函数
    EigenSystem(const EigenSystem&) = delete;

    // 禁用拷贝赋值运算符
    EigenSystem& operator=(const EigenSystem&) = delete;

    /**
     * @brief 将原始变量投影到特征空间。
     * @param prim 原始变量。
     * @param pp 两相物性参数。
     * @return 特征变量。
     */
    AMREX_GPU_HOST_DEVICE
    VarArray prim_to_char(VarArray const& prim, ProbParm const& pp) const;

    /**
     * @brief 将特征变量还原为原始变量。
     * @param charVars 特征变量。
     * @param pp 两相物性参数。
     * @return 原始变量。
     */
    AMREX_GPU_HOST_DEVICE
    VarArray char_to_prim(VarArray const& charVars, ProbParm const& pp) const;

    /**
     * @brief 将特征变量还原为守恒变量。
     * @param charVars 特征变量。
     * @param pp 两相物性参数。
     * @return 守恒变量。
     */
    AMREX_GPU_HOST_DEVICE
    VarArray char_to_cons(VarArray const& charVars, ProbParm const& pp) const;

    /**
     * @brief 访问左特征矩阵。
     * @return 左特征矩阵。
     */
    const VarMatrix& getLeftEig() const
    {
        return leftEig;
    }

    /**
     * @brief 访问右特征矩阵。
     * @return 右特征矩阵。
     */
    const VarMatrix& getRightEig() const
    {
        return rightEig;
    }

private:
    /**
     * @brief 计算 Roe 平均状态。
     * @param primL 左状态原始变量。
     * @param primR 右状态原始变量。
     * @param norm 面法向量。
     * @param pp 两相物性参数。
     * @return 平均后的辅助状态。
     */
    AMREX_GPU_HOST_DEVICE
    VarArray average(VarArray primL, VarArray primR, amrex::GpuArray<amrex::Real, 3> norm, ProbParm const& pp) const;

    // 左右特征矩阵
    VarMatrix leftEig;  // 左特征矩阵
    VarMatrix rightEig; // 右特征矩阵
};
