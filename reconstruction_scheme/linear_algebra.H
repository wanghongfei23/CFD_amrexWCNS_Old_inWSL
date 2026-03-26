#pragma once
#include <AMReX_BLassert.H>
#include <AMReX_Math.H>
#include <AMReX_REAL.H>
#include <vector>

// 遍历指定方向（i、j 或 k）并对每个索引位置执行给定操作
// 参数：
//   dir: 遍历的方向 (0: i 方向, 1: j 方向, 2: k 方向)
//   i, j, k: 起始位置索引
//   minIdx, maxIdx: 遍历范围的最小值和最大值
//   operation: 传入的操作，接受遍历方向上的索引（i/j/k）并执行相应操作
template <int dir, typename T>
AMREX_GPU_HOST_DEVICE void traverse_in_direction(int i, int j, int k, int minIdx, int maxIdx, T operation)
{
    // 编译时判断，避免运行时开销
    if constexpr (dir == 0) { // 遍历 i 方向
        for (int ix = i + minIdx; ix <= i + maxIdx; ++ix) {
            operation(ix, j, k); // 在 i 方向遍历，固定 j 和 k
        }
    } else if constexpr (dir == 1) { // 遍历 j 方向
        for (int jy = j + minIdx; jy <= j + maxIdx; ++jy) {
            operation(i, jy, k); // 在 j 方向遍历，固定 i 和 k
        }
    } else if constexpr (dir == 2) { // 遍历 k 方向
        for (int kz = k + minIdx; kz <= k + maxIdx; ++kz) {
            operation(i, j, kz); // 在 k 方向遍历，固定 i 和 j
        }
    }
}

template <int dir, typename T>
AMREX_GPU_HOST_DEVICE void traverse_in_direction_with_ii(int i, int j, int k, int minIdx, int maxIdx, T operation)
{
    // 编译时判断，避免运行时开销
    if constexpr (dir == 0) { // 遍历 i 方向
        for (int ii = 0; ii <= maxIdx - minIdx; ++ii) {
            operation(ii + i + minIdx, j, k, ii); // 在 i 方向遍历，固定 j 和 k
        }
    } else if constexpr (dir == 1) { // 遍历 j 方向
        for (int jj = 0; jj <= maxIdx - minIdx; ++jj) {
            operation(i, jj + j + minIdx, k, jj); // 在 j 方向遍历，固定 i 和 k
        }
    } else if constexpr (dir == 2) { // 遍历 k 方向
        for (int kk = 0; kk <= maxIdx - minIdx; ++kk) {
            operation(i, j, kk + k + minIdx, kk); // 在 k 方向遍历，固定 i 和 j
        }
    }
}

template <int dir>
AMREX_GPU_HOST_DEVICE constexpr amrex::GpuArray<amrex::Real, 3> get_direction()
{
    if constexpr (dir == 0) {
        return { 1, 0, 0 }; // dir=0 返回 {1, 0, 0}
    } else if constexpr (dir == 1) {
        return { 0, 1, 0 }; // dir=1 返回 {0, 1, 0}
    } else if constexpr (dir == 2) {
        return { 0, 0, 1 }; // dir=2 返回 {0, 0, 1}
    }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real, 3>
get_direction(int dir)
{
    amrex::GpuArray<amrex::Real, 3> result = { 0, 0, 0 }; // 默认返回 {0, 0, 0}

    if (dir == 0) {
        result = { 1, 0, 0 }; // dir=0 返回 {1, 0, 0}
    } else if (dir == 1) {
        result = { 0, 1, 0 }; // dir=1 返回 {0, 1, 0}
    } else if (dir == 2) {
        result = { 0, 0, 1 }; // dir=2 返回 {0, 0, 1}
    }

    return result;
}

template <int mode = 0, int nvar>
AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
        amrex::GpuArray<amrex::Real, nvar>
        LA_M_X(const amrex::GpuArray<amrex::Real, nvar * nvar>& M, const amrex::GpuArray<amrex::Real, nvar>& x) noexcept
{
    int n = x.size();
    amrex::GpuArray<amrex::Real, nvar> res;
    res.fill(0.0);

    if constexpr (mode == 0) {
        // 正常计算
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                res[i] += M[i * n + j] * x[j];
            }
        }
    } else if constexpr (mode == 1) {
        // 优先计算 x 的第二和第三个元素（j=1 和 j=2）
        for (int j = 1; j <= 2; j++) {
            for (int i = 0; i < n; i++) {
                res[i] += M[i * n + j] * x[j];
            }
        }
        // 计算剩余元素
        for (int j = 0; j < n; j++) {
            if (j == 1 || j == 2)
                continue; // 跳过已经计算的列
            for (int i = 0; i < n; i++) {
                res[i] += M[i * n + j] * x[j];
            }
        }
    } else if constexpr (mode == 2) {
        // 优先计算 x 的第三和第四个元素（j=2 和 j=3）
        for (int j = 2; j <= 3; j++) {
            for (int i = 0; i < n; i++) {
                res[i] += M[i * n + j] * x[j];
            }
        }
        // 计算剩余元素
        for (int j = 0; j < n; j++) {
            if (j == 2 || j == 3)
                continue; // 跳过已经计算的列
            for (int i = 0; i < n; i++) {
                res[i] += M[i * n + j] * x[j];
            }
        }
    } else {
        // 默认模式，正常计算
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                res[i] += M[i * n + j] * x[j];
            }
        }
    }

    return res;
}

namespace amrex {

template <class FAB>
void LinComb3(FabArray<FAB>& dst,
    typename FAB::value_type a, const FabArray<FAB>& x, int xcomp,
    typename FAB::value_type b, const FabArray<FAB>& y, int ycomp,
    typename FAB::value_type c, const FabArray<FAB>& z, int zcomp,
    int dstcomp, int numcomp, const IntVect& nghost)
{
    AMREX_ASSERT(dst.boxArray() == x.boxArray());
    AMREX_ASSERT(dst.distributionMap == x.distributionMap);
    AMREX_ASSERT(dst.boxArray() == y.boxArray());
    AMREX_ASSERT(dst.distributionMap == y.distributionMap);
    AMREX_ASSERT(dst.boxArray() == z.boxArray());
    AMREX_ASSERT(dst.distributionMap == z.distributionMap);
    AMREX_ASSERT(dst.nGrowVect().allGE(nghost) && x.nGrowVect().allGE(nghost) && y.nGrowVect().allGE(nghost) && z.nGrowVect().allGE(nghost));

    BL_PROFILE("amrex::LinComb()");

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && dst.isFusingCandidate()) {
        auto const& dstma = dst.arrays();
        auto const& xma = x.const_arrays();
        auto const& yma = y.const_arrays();
        auto const& zma = z.const_arrays();
        ParallelFor(dst, nghost, numcomp,
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k, int n) noexcept {
                dstma[box_no](i, j, k, dstcomp + n) = a * xma[box_no](i, j, k, xcomp + n)
                    + b * yma[box_no](i, j, k, ycomp + n)
                    + c * zma[box_no](i, j, k, zcomp + n);
            });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.growntilebox(nghost);
            auto const& xfab = x.const_array(mfi);
            auto const& yfab = y.const_array(mfi);
            auto const& zfab = z.const_array(mfi);
            auto const& dfab = dst.array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, numcomp, i, j, k, n,
                {
                    dfab(i, j, k, dstcomp + n) = a * xfab(i, j, k, xcomp + n) + b * yfab(i, j, k, ycomp + n) + c * zfab(i, j, k, zcomp + n);
                });
        }
    }
}
template <class FAB>
void LinComb4(FabArray<FAB>& dst,
    typename FAB::value_type a, const FabArray<FAB>& x, int xcomp,
    typename FAB::value_type b, const FabArray<FAB>& y, int ycomp,
    typename FAB::value_type c, const FabArray<FAB>& z, int zcomp,
    typename FAB::value_type d, const FabArray<FAB>& w, int wcomp,
    int dstcomp, int numcomp, const IntVect& nghost)
{
    // 确保所有输入数组的分布和区域一致
    AMREX_ASSERT(dst.boxArray() == x.boxArray());
    AMREX_ASSERT(dst.distributionMap == x.distributionMap);
    AMREX_ASSERT(dst.boxArray() == y.boxArray());
    AMREX_ASSERT(dst.distributionMap == y.distributionMap);
    AMREX_ASSERT(dst.boxArray() == z.boxArray());
    AMREX_ASSERT(dst.distributionMap == z.distributionMap);
    AMREX_ASSERT(dst.boxArray() == w.boxArray()); // 新增对w的检查
    AMREX_ASSERT(dst.distributionMap == w.distributionMap); // 新增对w的检查
    AMREX_ASSERT(dst.nGrowVect().allGE(nghost) && x.nGrowVect().allGE(nghost) && y.nGrowVect().allGE(nghost) && z.nGrowVect().allGE(nghost) && w.nGrowVect().allGE(nghost)); // 新增对w的检查

    BL_PROFILE("amrex::LinComb4()");

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && dst.isFusingCandidate()) {
        auto const& dstma = dst.arrays();
        auto const& xma = x.const_arrays();
        auto const& yma = y.const_arrays();
        auto const& zma = z.const_arrays();
        auto const& wma = w.const_arrays(); // 新增对w的引用
        ParallelFor(dst, nghost, numcomp,
            [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k, int n) noexcept {
                dstma[box_no](i, j, k, dstcomp + n) = a * xma[box_no](i, j, k, xcomp + n) + b * yma[box_no](i, j, k, ycomp + n) + c * zma[box_no](i, j, k, zcomp + n) + d * wma[box_no](i, j, k, wcomp + n); // 新增第四项
            });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(dst, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.growntilebox(nghost);
            auto const& xfab = x.const_array(mfi);
            auto const& yfab = y.const_array(mfi);
            auto const& zfab = z.const_array(mfi);
            auto const& wfab = w.const_array(mfi); // 新增对w的引用
            auto const& dfab = dst.array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, numcomp, i, j, k, n,
                {
                    dfab(i, j, k, dstcomp + n) = a * xfab(i, j, k, xcomp + n) + b * yfab(i, j, k, ycomp + n) + c * zfab(i, j, k, zcomp + n) + d * wfab(i, j, k, wcomp + n); // 新增第四项
                });
        }
    }
}

} // namespace amrex

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::GpuArray<amrex::Real, 3>
compute_normal_vector(int i, int j, int k,
    amrex::Array4<amrex::Real const> const& primin_array,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx)
{
    if constexpr (BL_SPACEDIM == 1) {
        return { 1.0, 0.0, 0.0 };
    }
    // 使用数组存储梯度分量，初始化为 0.0
    amrex::Real dpsi[3] = { 0.0 };

    // 计算梯度分量，使用 AMREX_D_TERM 宏适配不同维度
    AMREX_D_TERM(
        dpsi[0] = (primin_array(i + 1, j, k, Alpha1) - primin_array(i - 1, j, k, Alpha1)) / (2.0 * dx[0]);,
        dpsi[1] = (primin_array(i, j + 1, k, Alpha1) - primin_array(i, j - 1, k, Alpha1)) / (2.0 * dx[1]);
        ,
        dpsi[2] = (primin_array(i, j, k + 1, Alpha1) - primin_array(i, j, k - 1, Alpha1)) / (2.0 * dx[2]););

    // least square
    //  dpsi[0] = ((primin_array(i + 1, j + 1, k, Alpha1) - primin_array(i - 1, j + 1, k, Alpha1))
    //                + (primin_array(i + 1, j - 1, k, Alpha1) - primin_array(i - 1, j - 1, k, Alpha1))
    //                + (primin_array(i + 1, j, k, Alpha1) - primin_array(i - 1, j, k, Alpha1)))
    //      / (6.0 * dx[0]);
    //  dpsi[1] = ((primin_array(i + 1, j + 1, k, Alpha1) - primin_array(i + 1, j - 1, k, Alpha1))
    //                + (primin_array(i, j + 1, k, Alpha1) - primin_array(i, j - 1, k, Alpha1))
    //                + (primin_array(i - 1, j + 1, k, Alpha1) - primin_array(i - 1, j - 1, k, Alpha1)))
    //      / (6.0 * dx[1]);

    // 计算法向量的模长
    amrex::Real norm_factor = 0.0;
    for (int dim = 0; dim < BL_SPACEDIM; ++dim) {
        norm_factor += dpsi[dim] * dpsi[dim];
    }
    norm_factor = std::sqrt(norm_factor);

    // 归一化法向量
    amrex::GpuArray<amrex::Real, 3> n_l = { 0.0 };
    if (norm_factor > 1e-8) { // 防止除以零
        for (int dim = 0; dim < BL_SPACEDIM; ++dim) {
            n_l[dim] = -dpsi[dim] / norm_factor;
        }
    }

    return n_l;
}