/**
 * @file
 * @brief 一维单相 Euler 参考实现的状态描述符与边界注册实现。
 * @ingroup legacy_euler
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
 * 2. 设置气体绝热指数
 * 3. 读取参数并设置物理边界条件
 * 4. 添加状态描述符和组件
 */
void AmrLevelCong::variableSetUp()
{
    BL_ASSERT(desc_lst.size() == 0); // 确保描述符列表为空

    // 初始化包含问题特定变量的结构体
    h_prob_parm = new ProbParm {}; // 主机端问题参数
    d_prob_parm = (ProbParm*)The_Arena()->alloc(sizeof(ProbParm)); // 设备端问题参数
    h_prob_parm->gamma = 1.4; // 设置气体绝热指数

    // 将主机端参数复制到设备端
    amrex::Gpu::copy(amrex::Gpu::hostToDevice, h_prob_parm, h_prob_parm + 1, d_prob_parm);

    // 读取参数，设置物理边界条件
    read_params();

    // 添加状态描述符
    desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
        StateDescriptor::Point, NUM_GROW, NUM_STATE,
        &cell_cons_interp); // 使用守恒变量插值方法

    // 设置边界条件
    int lo_bc[BL_SPACEDIM]; // 下界边界条件
    int hi_bc[BL_SPACEDIM]; // 上界边界条件
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        lo_bc[i] = hi_bc[i] = BCType::foextrap; // 外推边界条件
    }

    BCRec bc(lo_bc, hi_bc); // 创建边界条件记录

    // 设置边界函数
    StateDescriptor::BndryFunc bndryfunc(nullfill);
    bndryfunc.setRunOnGPU(true); // 确保边界函数在GPU上运行

    // 设置状态变量组件
    desc_lst.setComponent(State_Type, 0, "Density", bc, bndryfunc); // 密度
    desc_lst.setComponent(State_Type, 1, "Velocity", bc, bndryfunc); // 速度
    desc_lst.setComponent(State_Type, 2, "TotalEnergy", bc, bndryfunc); // 总能量
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