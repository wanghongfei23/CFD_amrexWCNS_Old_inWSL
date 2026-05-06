/**
 * @file
 * @brief 基于特征分解的五阶重构实现。
 * @ingroup reconstruction
 */

#pragma once
#include "eigen_system.H"
#include "interpolation_schemes.H"
#include "linear_algebra.H"
#include "positive_preserving.H"
#include <AMReX_Array4.H>
#include <AMReX_BLassert.H>
// #include <vector>

// 模板函数 reconstruction：进行数值重构，计算特征值、特征向量，并返回拼接后的结果
// 参数：
//   i, j, k: 网格的索引位置
//   primin_array: 存储物理量（例如，压力、速度等）的 4D 数组
//   NUM_STATE: 物理状态量的数量（例如，速度、温度等）
//   norm: 用于计算方向的归一化向量，默认为 {1, 0, 0}，表示沿 i 方向的单位向量
template <int dir = 0>
AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        amrex::GpuArray<VarArray, 2>
        reconstruction(int i, int j, int k,
            amrex::Array4<amrex::Real const> const& primin_array,
            int NUM_STATE, ProbParm const& pp,
            amrex::GpuArray<amrex::Real, 3> norm,
            amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
{

    // 定义 offset 来控制某一维度的操作
    amrex::Dim3 offset = amrex::IntVect::TheDimensionVector(dir).dim3();

    // 初始化 primLR 向量，用于存储从 primin_array 中提取的物理量（大小为 NUM_STATE * 2）
    VarArray primCellL, primCellR;

    for (int ivar = 0; ivar < NUM_STATE; ivar++) {
        primCellL[ivar] = primin_array(i - 1 * offset.x, j - 1 * offset.y, k - 1 * offset.z, ivar); // 提取并存储到 primLR 中
        primCellR[ivar] = primin_array(i - 0 * offset.x, j - 0 * offset.y, k - 0 * offset.z, ivar); // 提取并存储到 primLR 中
    }

    // 使用 Roe 平均法对 primLR 进行处理，得到 primRef
    EigenSystem eig(primCellL, primCellR, norm, pp);

    // 初始化变量 it，qL 和 qR 用于存储左/右特征矩阵的变换结果
    amrex::GpuArray<amrex::GpuArray<amrex::Real, 5>, NVar> qL, qR;

    // 遍历指定方向，计算每个位置的物理量，并进行特征矩阵变换
    for (int ii = -3; ii <= 2; ii++) {
        VarArray primj;
        for (int ivar = 0; ivar < NUM_STATE; ivar++) {
            primj[ivar] = primin_array(i + ii * offset.x, j + ii * offset.y, k + ii * offset.z, ivar); // 获取当前点的物理量
        }

        int left_index = ii + 3;
        // 计算 charTemp 为特征矩阵的变换结果
        auto charTemp
            = eig.prim_to_char(primj, pp);

        // 将 charTemp 存入 qL 和 qR 中，用于后续计算
        if (left_index < 5) {
            for (int ivar = 0; ivar < NUM_STATE; ivar++) {
                qL[ivar][left_index] = charTemp[ivar];
            }
        }

        // 计算 qR（右特征矩阵），并存储变换结果
        int right_index = 5 - left_index;
        if (right_index < 5) {
            for (int ivar = 0; ivar < NUM_STATE; ivar++) {
                qR[ivar][right_index] = charTemp[ivar];
            }
        }
    }

    // 对左、右特征矩阵的数据进行 WENO5 重构
    // 【临时定位】格式选择
    VarArray charReconL, charReconR;
    for (int ivar = 0; ivar < NUM_STATE; ivar++) {
        charReconL[ivar] = TCNS_AS(qL[ivar]);
        charReconR[ivar] = TCNS_AS(qR[ivar]);
    }
    // TCNS_AS
    // weno_is

    // 使用 LA_M_X 进行矩阵-向量乘法，将计算结果存入 res
    VarArray consL = eig.char_to_cons(charReconL, pp); // 计算并将结果存储到 res 中
    VarArray consR = eig.char_to_cons(charReconR, pp); // 继续将结果存储到 res 中

    // 将原始变量 primCellL 和 primCellR 转换为守恒变量
    VarArray consCellL = prim_to_cons(primCellL, pp);
    VarArray consCellR = prim_to_cons(primCellR, pp);

    const PositiveQuantityFunctions pqf;

    // 对 consL 和 consR 应用保正限制器
    apply_positive_limiter(consCellL, consL, pqf, pp, NUM_STATE);
    apply_positive_limiter(consCellR, consR, pqf, pp, NUM_STATE);

    // 将保正后的守恒变量转换回原始变量
    VarArray primL = cons_to_prim(consL, pp);
    VarArray primR = cons_to_prim(consR, pp);

    // 返回拼接后的结果
    return { primL, primR };
}
