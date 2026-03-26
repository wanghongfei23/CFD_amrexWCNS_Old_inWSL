// 问题初始化实现文件
// 包含问题初始化相关的函数

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
