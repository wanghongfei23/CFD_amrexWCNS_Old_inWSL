/**
 * @file
 * @brief 二维单相 Euler 参考实现中的用户边界填充占位实现。
 * @ingroup legacy_euler
 */

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_PhysBCFunct.H>

using namespace amrex;

/**
 * 边界填充结构体
 * 用于处理计算域边界的虚拟细胞填充
 */
struct NullFill {
    /**
     * 运算符重载，用于执行边界填充操作
     * @param iv 单元索引
     * @param dest 目标数组
     * @param dcomp 目标组件起始索引
     * @param numcomp 组件数量
     * @param geom 几何信息
     * @param time 当前时间
     * @param bcr 边界条件记录
     * @param bcomp 边界组件
     * @param orig_comp 原始组件
     */
    AMREX_GPU_DEVICE
    void operator()(const IntVect& iv, Array4<Real> const& dest,
        const int dcomp, const int numcomp,
        GeometryData const& geom, const Real time,
        const BCRec* bcr, const int bcomp,
        const int orig_comp) const
    {
        // Lambda 函数：将原始变量转换为守恒变量
        const auto primitiveToConserved = [](Real rho, Real u, Real v, Real p, Real Gamma = 1.4) {
            return std::array<Real, 4> {
                rho, // 密度
                rho * u, // x 方向动量
                rho * v, // y 方向动量
                p / (Gamma - 1.0) + 0.5 * rho * (u * u + v * v) // 总能量
            };
        };

        auto ncomp = dest.nComp();

        // 获取计算域的边界索引
        const int ilo = geom.Domain().smallEnd(0); // x 方向下边界
        const int ihi = geom.Domain().bigEnd(0); // x 方向上边界
        const int jlo = geom.Domain().smallEnd(1); // y 方向下边界
        const int jhi = geom.Domain().bigEnd(1); // y 方向上边界
        const int klo = geom.Domain().smallEnd(2); // z 方向下边界
        const int khi = geom.Domain().bigEnd(2); // z 方向上边界

        // 获取当前单元的索引
        const auto [i, j, k] = iv.dim3();

        // 获取单元大小和计算域边界坐标
        const Real* dx = geom.CellSize();
        const Real* problo = geom.ProbLo();
        const Real* probhi = geom.ProbHi();

        // 计算当前单元的中心坐标
        Real x = problo[0] + (i + 0.5) * dx[0]; // x 坐标
        Real y = problo[1] + (j + 0.5) * dx[1]; // y 坐标
        Real z = problo[2] + (k + 0.5) * dx[2]; // z 坐标

        // 判断单元位于哪个边界
        bool is_left = x < problo[0]; // 左边界
        bool is_right = x > probhi[0]; // 右边界
        bool is_bottom = y < problo[1]; // 下边界
        bool is_top = y > probhi[1]; // 上边界
        bool is_back = z < problo[2]; // 后边界
        bool is_front = z > probhi[2]; // 前边界

        // 处理左边界
        if (is_left) {
            // 检查左边界是否设置为 ext_dir（外部边界条件）
            if (bcr[bcomp].lo(0) == BCType::ext_dir) {
                // 定义指定的原始变量
                constexpr Real rho = 8.0; // 密度
                constexpr Real u = 8.25 * cos(M_PI / 6); // x 方向速度
                constexpr Real v = -8.25 * sin(M_PI / 6); // y 方向速度
                constexpr Real p = 116.5; // 压力

                // 将原始变量转换为守恒变量
                auto conserved = primitiveToConserved(rho, u, v, p);

                // 将守恒变量赋值给 ghost cell
                for (int ivar = 0; ivar < ncomp; ++ivar) {
                    dest(i, j, k, ivar) = conserved[ivar];
                }
            }
        }

        // 处理右边界
        if (is_right) {
            // 处理右边界 ghost cell 的代码
        }

        // 处理下边界
        if (is_bottom) {
            // 检查下边界是否设置为 user_1（自定义边界条件）
            const Real* dx = geom.CellSize();
            const Real* problo = geom.ProbLo();

            // 检查 x 是否小于 1/6
            if (x < 1.0 / 6.0) {
                // 定义指定的原始变量
                constexpr Real rho = 8.0; // 密度
                constexpr Real u = 8.25 * cos(M_PI / 6); // x 方向速度
                constexpr Real v = -8.25 * sin(M_PI / 6); // y 方向速度
                constexpr Real p = 116.5; // 压力

                // 将原始变量转换为守恒变量
                auto conserved = primitiveToConserved(rho, u, v, p);

                // 将守恒变量赋值给虚拟细胞
                for (int ivar = 0; ivar < ncomp; ++ivar) {
                    dest(i, j, k, ivar) = conserved[ivar];
                }
            }
        }

        // 处理上边界
        if (is_top) {
            // 检查上边界是否设置为 user_1（自定义边界条件）
            if (bcr[bcomp].hi(1) == BCType::user_1) {
                // 计算动态阈值 gt
                Real gt = 1.0 / 6.0 + std::sqrt(3.0) / 3.0 * (1.0 + 20.0 * time);

                // 根据条件设置原始变量值
                std::array<Real, 4> exactValues;
                if (y - probhi[1] > std::sqrt(3.0) * (x - gt)) {
                    exactValues = { 8.0, 8.25 * cos(M_PI / 6), -8.25 * sin(M_PI / 6), 116.5 };
                } else {
                    exactValues = { 1.4, 0, 0, 1.0 };
                }

                // 将原始变量转换为守恒变量
                auto conserved = primitiveToConserved(exactValues[0], exactValues[1], exactValues[2], exactValues[3]);

                // 将守恒变量赋值给虚拟细胞
                for (int ivar = 0; ivar < ncomp; ++ivar) {
                    dest(i, j, k, ivar) = conserved[ivar];
                }
            }
        }

        // 处理后边界
        if (is_back) {
            // 处理后边界虚拟细胞的代码
        }

        // 处理前边界
        if (is_front) {
            // 处理前边界虚拟细胞的代码
        }
    }
};

/**
 * 边界填充函数
 * @param bx 边界框
 * @param data 数据数组
 * @param dcomp 目标组件起始索引
 * @param numcomp 组件数量
 * @param geom 几何信息
 * @param time 当前时间
 * @param bcr 边界条件记录
 * @param bcomp 边界组件
 * @param scomp 起始组件
 */
void nullfill(Box const& bx, FArrayBox& data,
    const int dcomp, const int numcomp,
    Geometry const& geom, const Real time,
    const Vector<BCRec>& bcr, const int bcomp,
    const int scomp)
{
    GpuBndryFuncFab<NullFill> gpu_bndry_func(NullFill {});
    gpu_bndry_func(bx, data, dcomp, numcomp, geom, time, bcr, bcomp, scomp);
}