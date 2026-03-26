// 标记参数实现文件
// 包含网格标记相关的参数读取函数

#include <AMReX_GpuMemory.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Vector.H>

#include "AmrLevelCong.H"

// 获取标记参数函数
// 功能：从输入文件读取网格标记相关的参数
void AmrLevelCong::get_tagging_params()
{
    // 使用ParmParse从输入文件获取级别数
    amrex::ParmParse pp("tagging");
    pp.query("max_phierr_lev", max_phierr_lev); // 最大phi误差级别
    pp.query("max_phigrad_lev", max_phigrad_lev); // 最大phi梯度级别

    // 设置误差阈值的默认值，然后从输入文件读取
    if (max_phierr_lev != -1) {
        phierr.resize(max_phierr_lev, 1.0e+20); // 调整数组大小并设置默认值
        pp.queryarr("phierr", phierr); // 从输入文件读取phi误差阈值
    }
    if (max_phigrad_lev != -1) {
        phigrad.resize(max_phigrad_lev, 1.0e+20); // 调整数组大小并设置默认值
        pp.queryarr("phigrad", phigrad); // 从输入文件读取phi梯度阈值
    }
}
