/**
 * @file
 * @brief 一维单相 Euler 参考实现的物性参数与变量编号。
 * @ingroup legacy_euler
 */

#ifndef PROB_PARM_H_
#define PROB_PARM_H_

#include <AMReX_Array.H>
#include <AMReX_REAL.H>

using namespace amrex::literals;

// 定义包含问题特定变量的结构体
struct ProbParm {
    amrex::Real gamma = 1.4; // 气体绝热指数
};

// 变量索引枚举
enum {
    Density = 0, // 密度
    XMom, // x方向动量
    Energy, // 能量
    NVar // 变量总数
};

// 变量数组类型定义
using VarArray = amrex::GpuArray<amrex::Real, NVar>;
// 变量矩阵类型定义
using VarMatrix = amrex::GpuArray<amrex::Real, NVar * NVar>;

#endif
