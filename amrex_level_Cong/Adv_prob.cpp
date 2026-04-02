/**
 * @file
 * @brief AMReX 问题初始化入口的空实现。
 * @details 当前求解器不依赖 probin 文件额外参数，因此该文件仅保留 AMReX 约定接口。
 * @ingroup amr_core
 */

#include <AMReX_REAL.H>

extern "C" {
    // 问题初始化函数
    // 输入：初始化标志init，问题名称name，名称长度namelen，问题下界problo，问题上界probhi
    // 功能：从probin文件读取额外输入（此处无需操作）
    void amrex_probinit (const int* /*init*/,
                         const int* /*name*/,
                         const int* /*namelen*/,
                         const amrex::Real* /*problo*/,
                         const amrex::Real* /*probhi*/)
    {
        // 这里不需要做任何事情，
        // 因为没有需要从probin文件读取的额外输入
    }
}
