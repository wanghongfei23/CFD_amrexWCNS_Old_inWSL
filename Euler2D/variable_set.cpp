// 变量设置实现文件（二维）
// 包含变量初始化和清理相关的函数

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
 * 3. 添加状态描述符和组件
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

    // 添加状态描述符
    desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
        StateDescriptor::Point, NUM_GROW, NUM_STATE,
        &cell_cons_interp); // 使用守恒变量插值方法

    // 使用ParmParse读取每个变量的边界条件
    ParmParse pp("amr");

    // 读取边界条件的辅助函数
    auto readBC = [&pp](const std::string& prefix) -> BCRec {
        int lo_bc[BL_SPACEDIM]; // 下界边界条件
        int hi_bc[BL_SPACEDIM]; // 上界边界条件
        for (int i = 0; i < BL_SPACEDIM; ++i) {
            // 动态生成边界条件名称
            std::string lo_bc_name = prefix + "_lo_bc_" + std::to_string(i);
            std::string hi_bc_name = prefix + "_hi_bc_" + std::to_string(i);

            // 如果未指定，默认为外推边界条件
            pp.query(lo_bc_name.c_str(), lo_bc[i]);
            pp.query(hi_bc_name.c_str(), hi_bc[i]);

            if (!pp.contains(lo_bc_name.c_str())) {
                lo_bc[i] = BCType::foextrap;
            }
            if (!pp.contains(hi_bc_name.c_str())) {
                hi_bc[i] = BCType::foextrap;
            }
        }
        return BCRec(lo_bc, hi_bc);
    };

    // 读取每个变量的边界条件
    BCRec density_bc = readBC("density"); // 密度边界条件
    BCRec xvelocity_bc = readBC("xvelocity"); // x方向速度边界条件
    BCRec yvelocity_bc = readBC("yvelocity"); // y方向速度边界条件
    BCRec energy_bc = readBC("energy"); // 能量边界条件

    // 设置边界函数
    StateDescriptor::BndryFunc bndryfunc(nullfill);
    bndryfunc.setRunOnGPU(true); // 确保边界函数在GPU上运行

    // 设置状态变量组件
    desc_lst.setComponent(State_Type, 0, "Density", density_bc, bndryfunc); // 密度
    desc_lst.setComponent(State_Type, 1, "XVelocity", xvelocity_bc, bndryfunc); // x方向速度
    desc_lst.setComponent(State_Type, 2, "YVelocity", yvelocity_bc, bndryfunc); // y方向速度
    desc_lst.setComponent(State_Type, 3, "TotalEnergy", energy_bc, bndryfunc); // 总能量
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