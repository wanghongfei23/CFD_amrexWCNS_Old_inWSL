#pragma once
#include <AMReX_MultiFab.H>
inline void Export1DEulerFieldToASCII(const amrex::MultiFab& mf, const amrex::Geometry& geom, const std::string& filename)
{
    // 打开输出文件
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        amrex::Abort("Error: Unable to open output file!");
    }

    // 获取网格的尺寸和物理范围
    const amrex::Real* dx = geom.CellSize(); // 网格单元的物理尺寸
    const amrex::Real* prob_lo = geom.ProbLo(); // 物理坐标起始点

    // 遍历 MultiFab 中的一维网格
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox(); // 当前网格的 Box
        auto const& arr = mf[mfi].array(); // 获取盒子内的 Array4 数据结构

        // 遍历一维网格点
        for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
            // 计算物理坐标
            amrex::Real x = prob_lo[0] + i * dx[0];

            // 输出物理坐标和流场变量
            outfile << x; // 输出 x 坐标

            // 遍历每个分量并输出
            for (int n = 0; n < mf.nComp(); ++n) {
                outfile << ' ' << arr(i, 0, 0, n); // 假设为 1D 情况，j=k=0
            }

            outfile << "\n"; // 换行
        }
    }

    // 关闭文件
    outfile.close();
}

inline void Export2DEulerFieldToASCII(const amrex::MultiFab& mf, const amrex::Geometry& geom, const std::string& filename)
{
    if (BL_SPACEDIM == 1) {
        Export1DEulerFieldToASCII(mf, geom, filename);
        return;
    }

    // 获取网格的尺寸和物理范围
    const amrex::Real* dx = geom.CellSize(); // 网格单元的物理尺寸
    const amrex::Real* prob_lo = geom.ProbLo(); // 物理坐标起始点

    // 获取网格的尺寸
    const amrex::Box& domain = geom.Domain();
    const int nx = domain.length(0); // x 方向网格点数
    const int ny = domain.length(1); // y 方向网格点数

    // 创建一个二维数组来存储数据
    std::vector<std::vector<std::vector<amrex::Real>>> data(nx, std::vector<std::vector<amrex::Real>>(ny, std::vector<amrex::Real>(mf.nComp())));

    // 遍历 MultiFab 中的二维网格
    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
        const amrex::Box& box = mfi.validbox(); // 当前网格的 Box
        auto const& arr = mf[mfi].array(); // 获取盒子内的 Array4 数据结构

        // 遍历当前分块中的网格点
        for (int j = box.smallEnd(1); j <= box.bigEnd(1); ++j) {
            for (int i = box.smallEnd(0); i <= box.bigEnd(0); ++i) {
                // 遍历每个分量并存储数据
                for (int n = 0; n < mf.nComp(); ++n) {
                    amrex::Real value = arr(i, j, 0, n); // 假设为 2D 情况，k=0

                    // 检测 NaN 并替换为 -1
                    if (std::isnan(value)) {
                        value = -1.0;
                        // amrex::Print() << i << ' ' << j << '\n';
                    }

                    data[i][j][n] = value; // 存储处理后的值
                }
            }
        }
    }

    // 打开输出文件
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        amrex::Abort("Error: Unable to open output file!");
    }

    // 写入 Tecplot 文件头
    outfile << "TITLE = \"2D Euler Field Data\"\n"; // 标题
    outfile << "VARIABLES = \"X\", \"Y\""; // 变量名：X 和 Y 坐标

    // 添加流场变量名
    for (int n = 0; n < mf.nComp(); ++n) {
        outfile << ", \"Var" << n + 1 << "\"";
    }
    outfile << "\n";

    // 写入区域信息
    outfile << "ZONE I=" << nx << ", J=" << ny << ", DATAPACKING=POINT\n";

    // 按照全局 (i, j) 顺序输出数据
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            // 计算物理坐标
            amrex::Real x = prob_lo[0] + i * dx[0]; // x 坐标
            amrex::Real y = prob_lo[1] + j * dx[1]; // y 坐标

            // 输出物理坐标
            outfile << x << ' ' << y;

            // 输出流场变量
            for (int n = 0; n < mf.nComp(); ++n) {
                outfile << ' ' << data[i][j][n];
            }

            outfile << "\n"; // 换行
        }
    }

    // 关闭文件
    outfile.close();
}

inline void ExportAllLevelsToFiles(amrex::Amr& amr)
{
    // 获取 AMR 的层级数量
    int num_levels = amr.finestLevel() + 1;

    // 遍历每个层级
    for (int lev = 0; lev < num_levels; ++lev) {
        // 获取当前层级的 MultiFab 和 Geometry
        amrex::MultiFab const& S_new = amr.getLevel(lev).get_new_data(0);
        amrex::Geometry const& geom = amr.Geom(lev);

        // 生成文件名，包含层级信息
        std::string filename = "euler_data_level_" + std::to_string(lev) + ".dat";

        // 调用导出函数，将数据写入文件
        Export2DEulerFieldToASCII(S_new, geom, filename);
    }
}

// 统计 Amr 对象中 level 0 的第一个变量的总和
inline amrex::Real SumFirstVariableOnLevel0(amrex::Amr& amr)
{
    amrex::Real sum = 0.0;

    // 获取 level 0 的 MultiFab
    const amrex::MultiFab& mf = amr.getLevel(0).get_new_data(0); // 假设第一个变量的索引为 0

    // 使用 Reduce 工具进行归约操作
    sum = amrex::ReduceSum(mf, 0, [=] AMREX_GPU_DEVICE(amrex::Box const& bx, amrex::Array4<const amrex::Real> const& arr) -> amrex::Real {
        amrex::Real local_sum = 0.0;

        // 遍历网格中的所有单元格
        amrex::Loop(bx, [&](int i, int j, int k) {
            local_sum += arr(i, j, k, 0); // 累加第一个变量的值
        });

        return local_sum;
    });

    return sum;
}

inline void WritePointDataToTecplot(amrex::Real time, const amrex::MultiFab& S_new, int i, int j, int k,
    const amrex::Real* var_values, const amrex::Real* residuals)
{
    // Generate filename based on point coordinates
    std::string filename = "point_data_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + ".dat";

    // Open file in append mode
    using namespace amrex;
    std::ofstream file(filename, std::ios::app);

    // Write header if the file is empty
    if (file.tellp() == 0) {
        file << "VARIABLES = \"time\", \"i\", \"j\", \"k\"";
        for (int n = 0; n < S_new.nComp(); ++n) {
            file << ", \"var" << n << "\""; // Variable names
        }
        for (int n = 0; n < S_new.nComp(); ++n) {
            file << ", \"res" << n << "\""; // Residual names
        }
        file << "\n";
    }

    // Write data for the specified point
    file << std::scientific << std::setprecision(6)
         << time << " " << i << " " << j << " " << k;

    // Write variable values
    for (int n = 0; n < S_new.nComp(); ++n) {
        file << " " << var_values[n];
    }

    // Write residuals
    for (int n = 0; n < S_new.nComp(); ++n) {
        file << " " << residuals[n];
    }

    file << "\n";

    // Close the file
    file.close();

    //     // Define points of interest
    // usage
    // int i1 = 59, j1 = 36, k1 = 0; // (i, j, k) = (59, 36, 0)
    // int i2 = 60, j2 = 36, k2 = 0; // (i, j, k) = (60, 36, 0)

    // // Get data for point 1
    // for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
    //     const Box& bx = mfi.validbox();
    //     if (bx.contains(IntVect(AMREX_D_DECL(i1, j1, k1)))) {
    //         const FArrayBox& fab = S_new[mfi];
    //         Real var_values[NUM_STATE];
    //         Real residuals[NUM_STATE];
    //         for (int n = 0; n < NUM_STATE; ++n) {
    //             var_values[n] = fab(IntVect(AMREX_D_DECL(i1, j1, k1)), n);
    //             residuals[n] = std::abs(var_values[n]); // Compute residual as absolute value
    //         }
    //         WritePointDataToTecplot(time, S_new, i1, j1, k1, var_values, residuals);
    //     }
    // }

    // // Get data for point 2
    // for (MFIter mfi(S_new); mfi.isValid(); ++mfi) {
    //     const Box& bx = mfi.validbox();
    //     if (bx.contains(IntVect(AMREX_D_DECL(i2, j2, k2)))) {
    //         const FArrayBox& fab = S_new[mfi];
    //         Real var_values[NUM_STATE];
    //         Real residuals[NUM_STATE];
    //         for (int n = 0; n < NUM_STATE; ++n) {
    //             var_values[n] = fab(IntVect(AMREX_D_DECL(i2, j2, k2)), n);
    //             residuals[n] = std::abs(var_values[n]); // Compute residual as absolute value
    //         }
    //         WritePointDataToTecplot(time, S_new, i2, j2, k2, var_values, residuals);
    //     }
    // }
}