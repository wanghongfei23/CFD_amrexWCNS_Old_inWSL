/**
 * @file
 * @brief 两相求解器状态描述符、物理参数与边界组件注册实现。
 * @ingroup two_phase_model
 */

#include "AmrLevelCong.H"
#include <AMReX_GpuMemory.H>
#include <AMReX_ParmParse.H>
/**
 * 定义数据描述符。
 */
using namespace amrex;

/**
 * 变量设置函数，用于初始化问题参数和数据描述符
 * 功能：
 * 1. 初始化问题参数结构体
 * 2. 读取参数并设置物理边界条件
 * 3. 定义状态变量数量
 * 4. 添加状态描述符和组件
 */
void AmrLevelCong::variableSetUp()
{
    BL_ASSERT(desc_lst.size() == 0); // 确保描述符列表为空

    // 初始化包含问题特定变量的结构体
    h_prob_parm = new ProbParm {}; // 主机端问题参数
    d_prob_parm = (ProbParm*)The_Arena()->alloc(sizeof(ProbParm)); // 设备端问题参数
    // 将主机端参数复制到设备端
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_prob_parm, h_prob_parm + 1, d_prob_parm);

    // 读取参数，设置物理边界条件
    read_params();

    // 定义状态变量数量：2个密度 + 空间维度的动量 + 1个能量 + 1个体积分数
    NUM_STATE = 4 + BL_SPACEDIM;

    // 添加状态描述符
    desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
        StateDescriptor::Point, NUM_GROW, NUM_STATE,
        &cell_cons_interp); // 使用守恒变量插值方法

    // 使用ParmParse读取每个变量的边界条件
    ParmParse pp("amr");

    // 读取边界条件的辅助函数
    auto readBC = [&pp](const std::string& prefix) -> BCRec {
        // 初始化边界条件为外推边界
        int lo_bc[BL_SPACEDIM] = { AMREX_D_DECL(amrex::BCType::foextrap,
            amrex::BCType::foextrap,
            amrex::BCType::foextrap) };
        int hi_bc[BL_SPACEDIM] = { AMREX_D_DECL(amrex::BCType::foextrap,
            amrex::BCType::foextrap,
            amrex::BCType::foextrap) };
        for (int i = 0; i < BL_SPACEDIM; ++i) {
            std::string lo_bc_name = prefix + "_lo_bc_" + std::to_string(i);
            std::string hi_bc_name = prefix + "_hi_bc_" + std::to_string(i);
            pp.query(lo_bc_name.c_str(), lo_bc[i]);
            pp.query(hi_bc_name.c_str(), hi_bc[i]);
        }
        return BCRec(lo_bc, hi_bc);
    };

    // 读取每个变量的边界条件
    BCRec density1_bc = readBC("density1"); // 第一相密度边界条件
    BCRec density2_bc = readBC("density2"); // 第二相密度边界条件
    BCRec energy_bc = readBC("energy"); // 能量边界条件
    BCRec volumeFraction1_bc = readBC("volumeFraction1"); // 第一相体积分数边界条件

    // 设置边界函数
    StateDescriptor::BndryFunc bndryfunc(nullfill);
    bndryfunc.setRunOnGPU(true); // 确保边界函数在GPU上运行

    // 向状态描述符添加组件
    desc_lst.setComponent(State_Type, 0, "Density1", density1_bc, bndryfunc); // 第一相密度
    desc_lst.setComponent(State_Type, 1, "Density2", density2_bc, bndryfunc); // 第二相密度

    // 根据空间维度添加动量组件
    int comp_idx = 2; // 动量组件的起始索引
    const char* dir_char = "XYZ";
    for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
        std::string momentum_name = "Momentum_" + std::string(1, dir_char[dir]); // X, Y, Z方向动量
        BCRec momentum_bc = readBC(momentum_name);
        desc_lst.setComponent(State_Type, comp_idx++, momentum_name, momentum_bc, bndryfunc);
    }

    // 添加剩余组件
    desc_lst.setComponent(State_Type, comp_idx++, "TotalEnergy", energy_bc, bndryfunc); // 总能量
    desc_lst.setComponent(State_Type, comp_idx++, "VolumeFraction1", volumeFraction1_bc, bndryfunc); // 第一相体积分数

    BL_ASSERT(comp_idx == NUM_STATE); // 确保组件数量与定义一致
}

/**
 * 运行结束时清理数据描述符。
 * 功能：
 * 1. 清空描述符列表
 * 2. 释放主机端和设备端的问题参数
 */
void AmrLevelCong::variableCleanUp()
{
    desc_lst.clear(); // 清空描述符列表

    // 删除包含问题特定参数的结构体
    delete h_prob_parm; // 释放主机端参数
    The_Arena()->free(d_prob_parm); // 释放设备端参数
}