/**
 * @file
 * @brief twoPhaseSolver 的主程序入口。
 * @details 负责初始化 AMReX、读取运行控制参数、执行主时间推进循环，并在运行前后导出调试数据。
 * @ingroup app_entry
 */

#include <iomanip>
#include <iostream>
#include <new>

#include <AMReX_Amr.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

#include "output_tmp.H"
#include <AMReX_MultiFab.H>
#include <fstream>
#include <iostream>

/**
 * @brief 返回全局 AMR 层级构建器。
 * @return 指向 LevelBld 实例的指针。
 */
amrex::LevelBld* getLevelBld();

#include <fstream> // 包含文件流头文件

/**
 * @brief 扫描所有 AMR 层级并报告 NaN 所在位置。
 * @param amr 当前求解器的 AMR 对象。
 * @return 发现任意 NaN 时返回 true。
 */
bool CheckForNaNAndPrint(Amr& amr)
{
    bool nanDetected = false;

    // 遍历所有层级
    for (int level = 0; level <= amr.finestLevel(); ++level) {
        const MultiFab& state = amr.getLevel(level).get_new_data(0); // 获取当前层级的状态数据
        const Geometry& geom = amr.Geom(level); // 获取几何对象

        // 使用 MFIter 遍历网格
        for (MFIter mfi(state); mfi.isValid(); ++mfi) {
            const Box& box = mfi.validbox();
            const FArrayBox& fab = state[mfi]; // 获取当前网格的数据

            // 遍历所有变量分量
            for (int comp = 0; comp < fab.nComp(); ++comp) {
                // 使用 ParallelFor 并行遍历 Box 中的单元
                amrex::LoopOnCpu(box, [&](int i, int j, int k) {
                    IntVect iv(AMREX_D_DECL(i, j, k));
                    if (std::isnan(fab(iv, comp))) { // 检查是否为 NaN
                        // 计算物理位置
                        amrex::Real physLoc[AMREX_SPACEDIM];
                        geom.CellCenter(iv, physLoc);

                        // 直接打印信息
                        std::cout << "NaN detected at the following location:\n";
                        std::cout << "  Simulation step: " << amr.levelSteps(0)
                                  << ", Simulation time: " << amr.cumTime() << "\n";
                        std::cout << "  Level: " << level << ", Index: " << iv
                                  << ", Component: " << comp << "\n";
                        std::cout << "  Physical Location: (";
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                            std::cout << physLoc[d];
                            if (d < AMREX_SPACEDIM - 1) {
                                std::cout << ", ";
                            }
                        }
                        std::cout << ")\n";
                        nanDetected = true;
                    }
                });
            }
        }
    }

    // 强制刷新输出缓冲区
    std::cout << std::flush;

    return nanDetected;
}

/**
 * @brief 运行 twoPhaseSolver 主流程。
 * @param argc 命令行参数个数。
 * @param argv 命令行参数数组。
 * @return 进程退出码。
 */
int main(int argc,
    char* argv[])
{
    amrex::Initialize(argc, argv); // 初始化AMReX
    amrex::Print() << "Dim= " << AMREX_SPACEDIM << '\n'; // 打印空间维度

    auto dRunTime1 = amrex::second(); // 记录开始时间

    int max_step; // 最大步数
    Real strt_time; // 开始时间
    Real stop_time; // 结束时间

    // 强制关闭GPU上的粒子分块
#ifdef AMREX_USE_GPU
    {
        ParmParse pp("particles");
        pp.add("do_tiling", 0);
    }
#endif

    // 读取参数
    {
        ParmParse pp;

        max_step = -1;
        strt_time = 0.0;
        stop_time = -1.0;

        pp.query("max_step", max_step);
        pp.query("strt_time", strt_time);
        pp.query("stop_time", stop_time);
    }

    // 检查参数有效性
    if (strt_time < 0.0) {
        amrex::Abort("MUST SPECIFY a non-negative strt_time");
    }

    if (max_step < 0 && stop_time < 0.0) {
        amrex::Abort("Exiting because neither max_step nor stop_time is non-negative.");
    }

    // 主模拟循环
    {
        Amr amr(getLevelBld()); // 创建Amr对象

        amr.init(strt_time, stop_time); // 初始化模拟

        // auto rho_start = SumFirstVariableOnLevel0(amr); // 记录初始密度
        // ExportAllLevelsToFiles(amr); // 导出所有层级到文件
        // volatile bool nanDetected = true;

        // 模拟主循环
        while (amr.okToContinue() && (amr.levelSteps(0) < max_step || max_step < 0) && (amr.cumTime() < stop_time || stop_time < 0.0)) {
            // 执行粗网格时间步，递归调用timeStep()
            amr.coarseTimeStep(stop_time);
            // bool nanDetected = CheckForNaNAndPrint(amr); // 检查NaN值
            // amrex::ParallelDescriptor::Barrier("Sync before Abort"); // 同步所有进程

            // // 如果检测到NaN值，终止模拟
            // if (nanDetected) {
            //     // 所有进程同步后调用 Abort
            //     amrex::Abort("NaN detected in the simulation.");
            //     break; // 防止后续代码执行（尽管 Abort 会直接退出）
            // }
        }

        // 写入最终的检查点和绘图文件
        if (amr.stepOfLastCheckPoint() < amr.levelSteps(0)) {
            amr.checkPoint();
        }

        if (amr.stepOfLastPlotFile() < amr.levelSteps(0)) {
            amr.writePlotFile();
        }
        // auto rho_end = SumFirstVariableOnLevel0(amr); // 记录最终密度

        // std::cout << "density loss= " << rho_start - rho_end << std::endl; // 打印密度损失

        // ExportAllLevelsToFiles(amr); // 导出所有层级到文件
    }

    auto dRunTime2 = amrex::second() - dRunTime1; // 计算运行时间

    ParallelDescriptor::ReduceRealMax(dRunTime2, ParallelDescriptor::IOProcessorNumber()); // 收集最大运行时间

    amrex::Print() << "Run time = " << dRunTime2 << std::endl; // 打印运行时间

    amrex::Finalize(); // 结束AMReX

    return 0;
}
