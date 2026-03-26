// 特征系统头文件
// 包含特征矩阵计算和特征变量转换的类定义

#pragma once

#include "Prob_Parm.H"
#include "equation.H"
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>

// 特征系统类，用于计算特征矩阵和处理特征变量转换
class EigenSystem {

public:
    // 构造函数，输入左右原始变量向量
    // 输入：左侧原始变量primL，右侧原始变量primR，法向量norm，问题参数pp
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

    // 原始变量转为特征变量
    // 输入：原始变量prim，问题参数pp
    // 输出：特征变量
    AMREX_GPU_HOST_DEVICE
    VarArray prim_to_char(VarArray const& prim, ProbParm const& pp) const;

    // 特征变量转为原始变量
    // 输入：特征变量charVars，问题参数pp
    // 输出：原始变量
    AMREX_GPU_HOST_DEVICE
    VarArray char_to_prim(VarArray const& charVars, ProbParm const& pp) const;

    // 特征变量转为守恒变量
    // 输入：特征变量charVars，问题参数pp
    // 输出：守恒变量
    AMREX_GPU_HOST_DEVICE
    VarArray char_to_cons(VarArray const& charVars, ProbParm const& pp) const;

    // 获取左特征矩阵
    // 输出：左特征矩阵
    const VarMatrix& getLeftEig() const
    {
        return leftEig;
    }

    // 获取右特征矩阵
    // 输出：右特征矩阵
    const VarMatrix& getRightEig() const
    {
        return rightEig;
    }

private:
    // Roe 平均方法
    // 输入：左侧原始变量primL，右侧原始变量primR，法向量norm，问题参数pp
    // 输出：平均后的变量
    AMREX_GPU_HOST_DEVICE
    VarArray average(VarArray primL, VarArray primR, amrex::GpuArray<amrex::Real, 3> norm, ProbParm const& pp) const;

    // 左右特征矩阵
    VarMatrix leftEig;  // 左特征矩阵
    VarMatrix rightEig; // 右特征矩阵
};
