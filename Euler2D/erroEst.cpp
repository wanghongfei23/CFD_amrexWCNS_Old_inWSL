/**
 * @file
 * @brief 二维单相 Euler 参考实现中的 AMR 标记逻辑。
 * @ingroup legacy_euler
 */

#include "AmrLevelCong.H"
#include "equation.H"
#include "teno_related.H"
/**
 * 用于网格重划分的误差估计。
 * 基于密度场的梯度计算误差，标记需要细化的网格区域
 */

/**
 * 计算标记值的函数，用于判断网格是否需要细化
 * @param i 网格点x坐标索引
 * @param j 网格点y坐标索引
 * @param k 网格点z坐标索引
 * @param S_new_array 状态数组
 * @return 是否需要标记为需要细化的布尔值
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE bool
computeTagValue(int i, int j, int k, amrex::Array4<amrex::Real const> const& S_new_array)
{
    bool res = false;

    // Lambda函数，用于计算给定方向上的标记
    auto computeFlag =
        [&](int di, int dj, int dk)
        -> AMREX_GPU_DEVICE bool {
        // 提取当前点及其周围点的密度值
        amrex::GpuArray<amrex::Real, 5> q = {
            S_new_array(i + di * (-2), j + dj * (-2), k + dk * (-2), 0),
            S_new_array(i + di * (-1), j + dj * (-1), k + dk * (-1), 0),
            S_new_array(i, j, k, 0),
            S_new_array(i + di * 1, j + dj * 1, k + dk * 1, 0),
            S_new_array(i + di * 2, j + dj * 2, k + dk * 2, 0)
        };
        // 使用TENO5格式计算标记
        return Teno5_CongZ_Flag(q);
    };

    // 计算各个维度上的标记
    if constexpr (SpaceDim >= 1) {
        res = res || computeFlag(1, 0, 0); // x方向
    }
    if constexpr (SpaceDim >= 2) {
        res = res || computeFlag(0, 1, 0); // y方向
    }
    if constexpr (SpaceDim >= 3) {
        res = res || computeFlag(0, 0, 1); // z方向
    }

    return res;
}

/**
 * 误差估计函数，用于标记需要细化的网格
 * @param tags 标记数组
 * @param clearval 清除值
 * @param tagval 标记值
 * @param time 当前时间
 * @param n_error_buf 误差缓冲区
 * @param ngrow 虚拟网格层数
 * @功能：根据密度场的梯度计算误差并标记需要细化的网格
 */
void AmrLevelCong::errorEst(amrex::TagBoxArray& tags,
    int /*clearval*/,
    int /*tagval*/,
    amrex::Real time,
    int /*n_error_buf*/,
    int /*ngrow*/)
{
    using namespace amrex;
    // 获取新的状态数据
    MultiFab& S_new = get_new_data(State_Type);

    // 为phi梯度检查正确填充补丁和虚拟细胞
    MultiFab density_tmp;
    if (level < max_phigrad_lev || true) {
        const Real cur_time = state[State_Type].curTime(); // 当前时间
        density_tmp.define(S_new.boxArray(), S_new.DistributionMap(), 1, 3); // 定义密度临时数组
        FillPatch(*this, density_tmp, 2, cur_time, State_Type, 0, 1); // 填充补丁
    }
    // 选择使用哪个密度数组
    MultiFab const& density = (level < max_phigrad_lev || true) ? density_tmp : S_new;

    const char tagval = TagBox::SET;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        // 遍历所有网格块
        for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& tilebx = mfi.tilebox(); // 当前瓦片的边界框
            Array4<Real const> const& S_new_array = S_new[mfi].const_array(); // 状态数组
            Array4<Real const> const& density_array = density[mfi].const_array(); // 密度数组
            auto tagarr = tags[mfi].array(); // 标记数组

            const auto& dx = geom.CellSize(); // 每个方向的网格尺寸
            const auto& problo = geom.ProbLo(); // 问题域的下界

            // 并行处理瓦片中的所有网格点
            amrex::ParallelFor(tilebx,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                    Real x = problo[0] + (i + 0.5) * dx[0]; // 计算当前网格中心的物理坐标
                    bool flag = computeTagValue<BL_SPACEDIM>(i, j, k, density_array); // 计算标记值
                    if (flag) //(x >= 0.6 && x <= 0.7) //(std::abs(density_array(i, j, k, 0) - 1.1) > 0.4) //
                        tagarr(i, j, k) = TagBox::SET; // 需要细化
                    else
                        tagarr(i, j, k) = TagBox::CLEAR; // 不需要细化
                    // std::cout << flag << '\n';
                });
        }
    }
    // std::cout << "error estimate finish \n";
}