#pragma once

#include "Prob_Parm.H"
#include "equation.H"

#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>

/**
 * 初始化一维两相算例的守恒量场
 * @param level 网格层级
 * @param time 当前时间
 * @param box 网格盒子
 * @param state 状态变量数组
 * @param dx 网格间距
 * @param problo 计算域左下角坐标
 * @param pp 两相物性参数
 */
#if BL_SPACEDIM == 1
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE void
initdata(int level, amrex::Real time,
    const amrex::Box& box, amrex::FArrayBox& state,
    const amrex::Real* dx, const amrex::Real* problo, const ProbParm& pp)
{
    using namespace amrex;

    // 原始变量数组
    VarArray prim;

    // 创建状态数组的 Array4 视图
    Array4<Real> state_arr = state.array();

    // 遍历网格框中的所有单元格（一维情况下仅 x 方向）
    const int i_lo = box.smallEnd(0);
    const int i_hi = box.bigEnd(0);

    for (int i = i_lo; i <= i_hi; ++i) {
        // 计算单元格中心在 x 方向的物理坐标
        Real x = problo[0] + (i + 0.5) * dx[0];

        // 根据选择的算例类型设置原始变量
        switch (pp.case_type) {
            case Case_1D_CenteredWater_Air_Moving: {
                // ============================================================================
                // 1D算例：基础一维两相流配置
                // ============================================================================
                // 编号: 1d_01
                // 描述: 中心区域为水，两侧为空气，整体以相同速度运动
                // 区域划分:
                //   - 水区域: 0.25 <= x < 0.75
                //     * 密度: 1000.0 kg/m^3
                //     * 速度: 100.0 m/s
                //     * 压力: 101325.0 Pa
                //     * 体积分数: alpha1 = 1.0
                //   - 空气区域: 其他区域
                //     * 密度: 1.204 kg/m^3
                //     * 速度: 100.0 m/s
                //     * 压力: 101325.0 Pa
                //     * 体积分数: alpha1 = 0.0
                // 激波: 无
                // 界面: 水与空气的分界面位于 x=0.25 和 x=0.75
                // ============================================================================
                constexpr Real epss = 1e-8;
                if (x >= 0.25 && x < 0.75) {
                    // 区域：0.25 <= x < 0.75（水）
                    prim[AlphaRho1] = 1000.0; // alpha1 * rho1 (kg/m^3)
                    prim[AlphaRho2] = epss; // alpha2 * rho2 (kg/m^3)
                    prim[XMom] = 100.0;       // u (m/s)
                    prim[Energy] = 101325.0;   // p (Pa)
                    prim[Alpha1] = 1.0 - epss; // alpha1
                } else {
                    // 区域：其他（空气）
                    prim[AlphaRho1] = epss;  // alpha1 * rho1 (kg/m^3)
                    prim[AlphaRho2] = 1.204;   // alpha2 * rho2 (kg/m^3)
                    prim[XMom] = 100.0;        // u (m/s)
                    prim[Energy] = 101325.0;    // p (Pa)
                    prim[Alpha1] = epss;     // alpha1
                }
                break;
            }
            case Case_1D_ShockTube_Water_Air: {
                // ============================================================================
                // 1D算例：一维两相激波管
                // ============================================================================
                // 编号: 1d_02
                // 描述: 左侧为高压水，右侧为低压空气，模拟激波管问题
                // 区域划分:
                //   - 高压水区域: x < 0.8
                //     * 密度: 1000.0 kg/m^3
                //     * 速度: 0.0 m/s
                //     * 压力: 1.0e9 Pa
                //     * 体积分数: alpha1 = 1.0
                //   - 低压空气区域: x >= 0.8
                //     * 密度: 20.0 kg/m^3
                //     * 速度: 0.0 m/s
                //     * 压力: 1.0e5 Pa
                //     * 体积分数: alpha1 = 0.0
                // 激波: 初始激波位于 x=0.8 处
                // 界面: 水与空气的分界面位于 x=0.8
                // ============================================================================
                constexpr Real epss = 1e-8;
                if (x < 0.8) {
                    // 区域：x < 0.8（高压水）
                    prim[AlphaRho1] = 1000.0; // alpha1 * rho1 (kg/m^3)
                    prim[AlphaRho2] = epss;   // alpha2 * rho2 (kg/m^3)
                    prim[XMom] = 0.0;         // u (m/s)
                    prim[Energy] = 1.0e9;     // p (Pa) - 高压
                    prim[Alpha1] = 1.0 - epss; // alpha1
                } else {
                    // 区域：x >= 0.8（低压空气）
                    prim[AlphaRho1] = epss;    // alpha1 * rho1 (kg/m^3)
                    prim[AlphaRho2] = 20.0;    // alpha2 * rho2 (kg/m^3)
                    prim[XMom] = 0.0;          // u (m/s)
                    prim[Energy] = 1.0e5;      // p (Pa) - 低压
                    prim[Alpha1] = epss;       // alpha1
                }
                break;
            }
            case Case_1D_ShockInterface_Water_Air: {
                // ============================================================================
                // 1D算例：激波与静止界面相互作用
                // ============================================================================
                // 编号: 1d_03
                // 描述: 三段区域配置，包含高速水流、静止水和空气
                // 区域划分:
                //   - 高速水流区域: 0.0 < x < 1.31
                //     * 密度: 1230.0 kg/m^3
                //     * 速度: 432.7 m/s
                //     * 压力: 1.0e9 Pa
                //     * 体积分数: alpha1 = 1.0
                //   - 静止水区域: 1.31 < x < 2.0
                //     * 密度: 1000.0 kg/m^3
                //     * 速度: 0.0 m/s
                //     * 压力: 1.0e5 Pa
                //     * 体积分数: alpha1 = 1.0
                //   - 空气区域: 其他区域
                //     * 密度: 1.0 kg/m^3
                //     * 速度: 0.0 m/s
                //     * 压力: 1.0e5 Pa
                //     * 体积分数: alpha1 = 0.0
                // 激波: 初始激波位于高速水流与静止水的界面处
                // 界面: 高速水流与静止水的分界面位于 x=1.31，静止水与空气的分界面位于 x=2.0
                // ============================================================================
                constexpr Real epss = 1e-8;
                if (x > 0.0 && x < 1.31) {
                    // 区域：0.0 < x < 1.31（高速水流）
                    prim[AlphaRho1] = 1230.0; // rho1 * alpha1 (kg/m^3)
                    prim[AlphaRho2] = epss;   // rho2 * alpha2 (kg/m^3)
                    prim[XMom] = 432.7;       // u (m/s) - 高速
                    prim[Energy] = 1.0e9;     // p (Pa) - 高压
                    prim[Alpha1] = 1.0 - epss; // alpha1
                } else if (x > 1.31 && x < 2.0) {
                    // 区域：1.31 < x < 2.0（静止水）
                    prim[AlphaRho1] = 1000.0; // rho1 * alpha1 (kg/m^3)
                    prim[AlphaRho2] = epss; // rho2 * alpha2 (kg/m^3)
                    prim[XMom] = 0.0;         // u (m/s) - 静止
                    prim[Energy] = 1.0e5;     // p (Pa) - 常压
                    prim[Alpha1] = 1.0 - epss; // alpha1
                } else {
                    // 区域：其他（空气）
                    prim[AlphaRho1] = epss;    // rho1 * alpha1 (kg/m^3)
                    prim[AlphaRho2] = 1.0;     // rho2 * alpha2 (kg/m^3)
                    prim[XMom] = 0.0;          // u (m/s) - 静止
                    prim[Energy] = 1.0e5;      // p (Pa) - 常压
                    prim[Alpha1] = epss;       // alpha1
                }
                break;
            }

        }

        // 将原始变量转换为守恒变量
        VarArray cons = prim_to_cons(prim, pp);
        // 使用 Array4 将守恒变量存储到状态数组中
        state_arr(i, 0, 0, AlphaRho1) = cons[AlphaRho1];
        state_arr(i, 0, 0, AlphaRho2) = cons[AlphaRho2];
        state_arr(i, 0, 0, XMom) = cons[XMom];
        state_arr(i, 0, 0, Energy) = cons[Energy];
        state_arr(i, 0, 0, Alpha1) = cons[Alpha1];
    }
}
#endif

#if BL_SPACEDIM == 2
/**
 * 初始化二维两相算例的守恒量场
 * @param level 网格层级
 * @param time 当前时间
 * @param box 网格盒子
 * @param state 状态变量数组
 * @param dx 网格间距
 * @param problo 计算域左下角坐标
 * @param pp 两相物性参数
 */
amrex::GpuArray<amrex::Real, 3>
compute_normal_vector1(int i, int j, int k,
    amrex::Array4<amrex::Real const> const& primin_array,
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx);

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE void
initdata(int level, amrex::Real time,
    const amrex::Box& box, amrex::FArrayBox& state,
    const amrex::Real* dx, const amrex::Real* problo, const ProbParm& pp)
{
    using namespace amrex;

    // 原始变量数组
    VarArray prim;

    // 创建状态数组的 Array4 视图
    Array4<Real> state_arr = state.array();

    // 遍历网格框中的所有单元格（二维情况）
    const int i_lo = box.smallEnd(0);
    const int i_hi = box.bigEnd(0);
    const int j_lo = box.smallEnd(1);
    const int j_hi = box.bigEnd(1);

    for (int j = j_lo; j <= j_hi; ++j) {
        for (int i = i_lo; i <= i_hi; ++i) {
            // 计算单元格中心在 x 和 y 方向的物理坐标
            Real x = problo[0] + (i + 0.5) * dx[0];
            Real y = problo[1] + (j + 0.5) * dx[1];

            // 根据选择的算例类型设置原始变量
            switch (pp.case_type) {
                case Case_2D_WaterColumn_Air_Stationary: {
                    // ============================================================================
                    // 2D算例：基础二维水柱配置
                    // ============================================================================
                    // 编号: 2d_01
                    // 描述: 简单的水柱和空气区域划分，无激波
                    // 区域划分:
                    //   - 水柱内部: 中心(0.04, 0.00) m，半径0.011 m
                    //     * 密度: 1000.0 kg/m^3
                    //     * 速度: (0.0, 0.0) m/s
                    //     * 压力: 1.01e5 Pa
                    //     * 体积分数: alpha1 = 1.0
                    //   - 空气区域: 水柱外部
                    //     * 密度: 1.17 kg/m^3
                    //     * 速度: (0.0, 0.0) m/s
                    //     * 压力: 1.01e5 Pa
                    //     * 体积分数: alpha1 = 0.0
                    // 激波: 无
                    // 界面: 水与空气的分界面为圆形，半径0.011 m
                    // ============================================================================
                    constexpr Real water_cylinder_x_center = 0.04; // 水柱中心 x 坐标 [m]
                    constexpr Real water_cylinder_y_center = 0.00; // 水柱中心 y 坐标 [m]
                    constexpr Real water_cylinder_radius = 0.011; // 水柱半径 [m]
                    constexpr Real epss = 1e-8;
                    
                    Real distance = std::sqrt((x - water_cylinder_x_center) * (x - water_cylinder_x_center) + 
                                           (y - water_cylinder_y_center) * (y - water_cylinder_y_center));
                                           
                    if (distance <= water_cylinder_radius) {
                        // 区域：水柱内部
                        prim[AlphaRho1] = 1000.0; // α1ρ1 (kg/m^3)
                        prim[AlphaRho2] = epss; // α2ρ2 (kg/m^3)
                        prim[XMom] = 0.0;         // u (m/s)
                        prim[YMom] = 0.0;         // v (m/s)
                        prim[Energy] = 1.01e5;     // p (Pa)
                        prim[Alpha1] = 1.0 - epss; // α1
                    } else { 
                        // 区域：水柱外部（空气）
                        prim[AlphaRho1] = epss;  // α1ρ1 (kg/m^3)
                        prim[AlphaRho2] = 1.17;    // α2ρ2 (kg/m^3)
                        prim[XMom] = 0.0;          // u (m/s)
                        prim[YMom] = 0.0;          // v (m/s)
                        prim[Energy] = 1.01e5;      // p (Pa)
                        prim[Alpha1] = epss;     // α1
                    }
                    break;
                }
                case Case_2D_ShockImpact_Water_Air_Mach2_4: {
                    // ============================================================================
                    // 2D算例：Mach 2.4 激波撞击水柱
                    // ============================================================================
                    // 编号: 2d_02
                    // 描述: 左侧为高速高压空气（激波后），中间为水柱，右侧为静止空气
                    // 区域划分:
                    //   - 激波后空气: x <= 0.029
                    //     * 密度: 3.7579 kg/m^3
                    //     * 速度: (574.57, 0.0) m/s
                    //     * 压力: 6.6189e5 Pa
                    //     * 体积分数: alpha1 = 0.0
                    //   - 水柱内部: 中心(0.04, 0.00) m，半径0.011 m
                    //     * 密度: 1000.0 kg/m^3
                    //     * 速度: (0.0, 0.0) m/s
                    //     * 压力: 1.01e5 Pa
                    //     * 体积分数: alpha1 = 1.0
                    //   - 静止空气: 水柱外部
                    //     * 密度: 1.17 kg/m^3
                    //     * 速度: (0.0, 0.0) m/s
                    //     * 压力: 1.01e5 Pa
                    //     * 体积分数: alpha1 = 0.0
                    // 激波: 初始激波位于 x=0.029 m 处
                    // 界面: 水与空气的分界面为圆形，半径0.011 m
                    // ============================================================================
                    constexpr Real water_cylinder_x_center = 0.04; // 水柱中心 x 坐标 [m]
                    constexpr Real water_cylinder_y_center = 0.00; // 水柱中心 y 坐标 [m]
                    constexpr Real water_cylinder_radius = 0.011; // 水柱半径 [m]
                    constexpr Real epss = 1e-8;

                    Real distance = std::sqrt((x - water_cylinder_x_center) * (x - water_cylinder_x_center) + 
                                           (y - water_cylinder_y_center) * (y - water_cylinder_y_center));
                                           
                    if (x <= 0.029) {
                        // 区域：x <= 0.029（激波后空气）
                        prim[AlphaRho1] = epss;    // α1ρ1 (kg/m^3)
                        prim[AlphaRho2] = 3.7579;  // α2ρ2 (kg/m^3)
                        prim[XMom] = 574.57;       // u (m/s) - 高速
                        prim[YMom] = 0.0;          // v (m/s)
                        prim[Energy] = 6.6189e5;    // p (Pa) - 高压
                        prim[Alpha1] = epss;       // α1
                    } else if (distance <= water_cylinder_radius + 1e-8) {
                        // 区域：水柱内部
                        prim[AlphaRho1] = 1000.0; // α1ρ1 (kg/m^3)
                        prim[AlphaRho2] = epss;   // α2ρ2 (kg/m^3)
                        prim[XMom] = 0.0;         // u (m/s) - 静止
                        prim[YMom] = 0.0;         // v (m/s) - 静止
                        prim[Energy] = 1.01e5;     // p (Pa) - 常压
                        prim[Alpha1] = 1.0 - epss; // α1
                    } else {
                        // 区域：水柱外部（静止空气）
                        prim[AlphaRho1] = epss;    // α1ρ1 (kg/m^3)
                        prim[AlphaRho2] = 1.17;    // α2ρ2 (kg/m^3)
                        prim[XMom] = 0.0;          // u (m/s) - 静止
                        prim[YMom] = 0.0;          // v (m/s) - 静止
                        prim[Energy] = 1.01e5;      // p (Pa) - 常压
                        prim[Alpha1] = epss;       // α1
                    }
                    break;
                }
                case Case_2D_ShockImpact_Water_Air_Mach10: {
                    // ============================================================================
                    // 2D算例：马赫 10 激波撞击水柱
                    // ============================================================================
                    // 编号: 2d_03
                    // 描述: 左侧为极高速度高压空气（激波后），中间为小水柱，右侧为静止空气
                    // 区域划分:
                    //   - 激波后空气: x <= 0.009
                    //     * 密度: 6.8571 kg/m^3
                    //     * 速度: (2817.9, 0.0) m/s
                    //     * 压力: 1.165e7 Pa
                    //     * 体积分数: alpha1 = 0.0
                    //   - 水柱内部: 中心(0.013, 0.00) m，半径0.004 m
                    //     * 密度: 1000.0 kg/m^3
                    //     * 速度: (0.0, 0.0) m/s
                    //     * 压力: 1.0e5 Pa
                    //     * 体积分数: alpha1 = 1.0
                    //   - 激波前空气: x > 0.009
                    //     * 密度: 1.2 kg/m^3
                    //     * 速度: (0.0, 0.0) m/s
                    //     * 压力: 1.0e5 Pa
                    //     * 体积分数: alpha1 = 0.0
                    // 激波: 初始激波位于 x=0.009 m 处
                    // 界面: 水与空气的分界面为圆形，半径0.004 m
                    // ============================================================================
                    constexpr Real water_cylinder_x_center = 0.013; // 水柱中心 x 坐标 [m]
                    constexpr Real water_cylinder_y_center = 0.00; // 水柱中心 y 坐标 [m]
                    constexpr Real water_cylinder_radius = 0.004;  // 水柱半径 [m]
                    constexpr Real epss = 1e-8;

                    Real distance = std::sqrt((x - water_cylinder_x_center) * (x - water_cylinder_x_center) + 
                                           (y - water_cylinder_y_center) * (y - water_cylinder_y_center));
                                           
                    if (distance <= water_cylinder_radius) {
                        // 区域：水柱内部
                        prim[AlphaRho1] = 1000.0; // α1ρ1 (kg/m^3)
                        prim[AlphaRho2] = epss;   // α2ρ2 (kg/m^3)
                        prim[XMom] = 0.0;         // u (m/s) - 静止
                        prim[YMom] = 0.0;         // v (m/s) - 静止
                        prim[Energy] = 1.0e5;      // p (Pa) - 常压
                        prim[Alpha1] = 1.0 - epss; // α1
                    } else if (x > 0.009 - 1e-8) {
                        // 区域：x > 0.009（激波前空气）
                        prim[AlphaRho1] = epss;    // α1ρ1 (kg/m^3)
                        prim[AlphaRho2] = 1.2;     // α2ρ2 (kg/m^3)
                        prim[XMom] = 0.0;          // u (m/s) - 静止
                        prim[YMom] = 0.0;          // v (m/s) - 静止
                        prim[Energy] = 1.0e5;      // p (Pa) - 常压
                        prim[Alpha1] = epss;       // α1
                    } else {
                        // 区域：x <= 0.009（激波后空气）
                        prim[AlphaRho1] = epss;     // α1ρ1 (kg/m^3)
                        prim[AlphaRho2] = 6.8571;   // α2ρ2 (kg/m^3)
                        prim[XMom] = 2.8179e3;      // u (m/s) - 极高速度
                        prim[YMom] = 0.0;           // v (m/s)
                        prim[Energy] = 1.165e7;     // p (Pa) - 极高压力
                        prim[Alpha1] = epss;        // α1
                    }
                    break;
                }
                case Case_2D_RMI_Air_SF6_Perturbed: {
                    // ============================================================================
                    // 2D算例：RMI_SF6 case1-半边截断
                    // ============================================================================
                    // 编号: 2d_04
                    // 描述: 2D Richtmyer–Meshkov 不稳定性（空气–SF6）
                    // 区域划分:
                    //   - SF6区域: x < x_interface（重流体，左侧）
                    //     * 密度: 5.04 kg/m^3
                    //     * 速度: (0.0, 0.0) m/s
                    //     * 压力: 1.0
                    //     * 体积分数: alpha1 = 0.0
                    //   - 激波后空气: x > x_shock（已受冲击空气，右侧）
                    //     * 密度: 1.411 kg/m^3
                    //     * 速度: (-0.39, 0.0) m/s
                    //     * 压力: 1.628
                    //     * 体积分数: alpha1 = 1.0
                    //   - 激波前空气: 介于界面与激波面之间
                    //     * 密度: 1.0 kg/m^3
                    //     * 速度: (0.0, 0.0) m/s
                    //     * 压力: 1.0
                    //     * 体积分数: alpha1 = 1.0
                    // 激波: 初始激波位于 x=3.2 处（激波从右向左传播）
                    // 界面: 扰动界面 x_interface = 2.9 - 0.2 * sin(2π(y + 0.25))
                    // ============================================================================
                    constexpr Real x_shock = 3.2;          // 激波面位置 x_shock
                    constexpr Real x_1 = 2.9;          // 基准位置 x_1
                    constexpr Real varepsilon = 0.2;          // 扰动振幅 varepsilon
                    constexpr Real epss = 1e-8;            // 小量用于避免浮点误差

                    // 计算扰动界面 x 坐标（作为 y 的函数）
                    Real x_interface = x_1 - varepsilon * std::sin(2.0 * M_PI * (y + 0.25));

                    // 判断空间位置 → 三区域划分
                    if (x < x_interface) {
                        // 区域 1: SF6（重流体，左侧）
                        prim[AlphaRho1] = epss;          // α₁ρ₁ ≈ 0 (空气几乎为0)
                        prim[AlphaRho2] = 5.04;          // α₂ρ₂ = ρ_SF6
                        prim[XMom]      = 0.0;           // u = 0
                        prim[YMom]      = 0.0;           // v = 0
                        prim[Energy]    = 1.0;           // p = 1.0
                        prim[Alpha1]    = epss;          // α₁ ≈ 0
                    } else if (x > x_shock + epss) {
                        // 区域 2: 已受冲击空气（激波后，右侧）
                        prim[AlphaRho1] = 1.411;         // α₁ρ₁ = ρ_post_air
                        prim[AlphaRho2] = epss;          // α₂ρ₂ ≈ 0
                        prim[XMom]      = -0.39;         // u = -0.39（向左）
                        prim[YMom]      = 0.0;           // v = 0
                        prim[Energy]    = 1.628;         // p = 1.628
                        prim[Alpha1]    = 1.0 - epss;    // α₁ ≈ 1
                    } else {
                        // 区域 3: 未受冲击空气（激波前，介于界面与激波面之间）
                        prim[AlphaRho1] = 1.0;           // α₁ρ₁ = ρ_air
                        prim[AlphaRho2] = epss;          // α₂ρ₂ ≈ 0
                        prim[XMom]      = 0.0;           // u = 0
                        prim[YMom]      = 0.0;           // v = 0
                        prim[Energy]    = 1.0;           // p = 1.0
                        prim[Alpha1]    = 1.0 - epss;    // α₁ ≈ 1
                    }
                    break;
                }
                case Case_2D_RMI_SF6_Air_Square: {
                    // ============================================================================
                    // 2D算例：RMI_SF6 方块
                    // ============================================================================
                    // 编号: 2d_05
                    // 描述: 2D Richtmyer–Meshkov 不稳定性（SF6方块–空气）
                    // 区域划分:
                    //   - SF6方块: 中心(0.5, 0.0) m，边长0.4 m
                    //     * 密度: 5.04 kg/m^3
                    //     * 速度: (0.0, 0.0) m/s
                    //     * 压力: 1.0
                    //     * 体积分数: alpha1 = 1.0
                    //   - 激波后空气: x < x_shock（左侧）
                    //     * 密度: 1.411 kg/m^3
                    //     * 速度: (0.39, 0.0) m/s
                    //     * 压力: 1.628
                    //     * 体积分数: alpha1 = 0.0
                    //   - 激波前空气: 激波前，介于激波面与方块之间
                    //     * 密度: 1.0 kg/m^3
                    //     * 速度: (0.0, 0.0) m/s
                    //     * 压力: 1.0
                    //     * 体积分数: alpha1 = 0.0
                    // 激波: 初始激波位于 x=0.2 处（激波马赫数：1.2）
                    // 界面: 正方形界面，中心(0.5, 0.0) m，边长0.4 m
                    // ============================================================================
                    constexpr Real x_shock = 0.2;     // 激波初始位置（左侧）
                    constexpr Real epss    = 1e-8;    // 小量

                    // -------------------------
                    // 理想方块界面参数
                    // -------------------------
                    constexpr Real x_c = 0.5;         // 方块中心 x 坐标
                    constexpr Real y_c = 0.0;         // 方块中心 y 坐标
                    constexpr Real L   = 0.4;         // 方块边长

                    // 理想方块内部判定
                    bool in_square = (std::abs(x - x_c) <= 0.5 * L) &&
                                    (std::abs(y - y_c) <= 0.5 * L);

                    // 分区：
                    // x < x_shock         -> 激波后空气（左侧）
                    // x_shock ~ square    -> 激波前空气
                    // in_square           -> SF6方块
                    if (in_square) {
                        // 区域 1: SF6方块内部（重流体）
                        prim[AlphaRho1] = 5.04;       // α₁ρ₁ = ρ_SF6 (相1现在是SF6)
                        prim[AlphaRho2] = epss;       // α₂ρ₂ ≈ 0
                        prim[XMom]      = 0.0;
                        prim[YMom]      = 0.0;
                        prim[Energy]    = 1.0;
                        prim[Alpha1]    = 1.0 - epss; // α₁ ≈ 1 (SF6体积分数接近1)
                    }
                    else if (x < x_shock - epss) {
                        // 区域 2: 已受冲击空气（激波后，左侧，轻流体）
                        prim[AlphaRho1] = epss;       // α₁ρ₁ ≈ 0 (SF6几乎为0)
                        prim[AlphaRho2] = 1.411;      // α₂ρ₂ = ρ_post_air (相2是空气)
                        prim[XMom]      = 0.39;       // 向右传播时，速度应为正
                        prim[YMom]      = 0.0;
                        prim[Energy]    = 1.628;
                        prim[Alpha1]    = epss;       // α₁ ≈ 0 (SF6体积分数接近0)
                    }
                    else {
                        // 区域 3: 未受冲击空气（激波前，轻流体）
                        prim[AlphaRho1] = epss;       // α₁ρ₁ ≈ 0
                        prim[AlphaRho2] = 1.0;        // α₂ρ₂ = ρ_air (相2是空气)
                        prim[XMom]      = 0.0;
                        prim[YMom]      = 0.0;
                        prim[Energy]    = 1.0;
                        prim[Alpha1]    = epss;       // α₁ ≈ 0
                    }
                    break;
                }
                case Case_2D_RMI_SF6_Air_Sector: {
                    // ============================================================================
                    // 2D算例：扇形区域内的RMI问题（参考RMI_SF6 case1-半边截断）
                    // ============================================================================
                    // 编号: 2d_06
                    // 描述: 90度扇形区域内的Richtmyer–Meshkov不稳定性（空气–SF6）
                    // 区域划分:
                    //   - 扇形内部:
                    //     * SF6区域: 靠近扇形中心，重流体（环形）
                    //     * 激波后空气: 靠近扇形边缘，已受冲击（环形）
                    //     * 激波前空气: 介于界面与激波面之间（环形）
                    //   - 扇形外部: 其他区域
                    // 激波: 初始激波位于扇形边缘（环形）
                    // 界面: 扰动界面（环形）
                    // ============================================================================
                    constexpr Real sector_center_x = 0.0;    // 扇形中心 x 坐标 [m]
                    constexpr Real sector_center_y = 0.0;    // 扇形中心 y 坐标 [m]
                    constexpr Real sector_radius_inner = 0.0; // 扇形内半径 [m] - 与计算域边界一致
                    constexpr Real sector_radius_outer = 3.5; // 扇形外半径 [m]
                    constexpr Real sector_angle_start = 0.0;   // 扇形起始角度 [弧度]
                    constexpr Real sector_angle_end = M_PI / 2; // 扇形结束角度 [弧度] (90度)
                    
                    // RMI参数（环形设置）
                    constexpr Real shock_radius = 3.2;          // 激波面半径
                    constexpr Real interface_radius = 2.9;      // 基准界面半径
                    constexpr Real varepsilon = 0.2;            // 扰动振幅
                    constexpr Real epss = 1e-8;                // 小量用于避免浮点误差

                    // 计算到中心的距离
                    Real distance = std::sqrt((x - sector_center_x) * (x - sector_center_x) + 
                                           (y - sector_center_y) * (y - sector_center_y));
                    
                    // 计算角度（相对于x轴正方向）
                    Real angle = std::atan2(y - sector_center_y, x - sector_center_x);
                    // 确保角度在[0, 2π)范围内
                    if (angle < 0) angle += 2 * M_PI;

                    // 判断是否在扇形区域内
                    bool in_sector = (distance >= sector_radius_inner) && 
                                    (distance <= sector_radius_outer) && 
                                    (angle >= sector_angle_start) && 
                                    (angle <= sector_angle_end);

                    if (in_sector) {
                        // 在扇形区域内，使用环形设置
                        // 计算扰动界面半径（作为角度的函数）
                        Real interface_radius_perturbed = interface_radius - varepsilon * std::sin(2.0 * M_PI * (angle + 0.25));

                        // 判断空间位置 → 三区域划分（环形）
                        if (distance < interface_radius_perturbed) {
                            // 区域 1: SF6（重流体，靠近中心）
                            prim[AlphaRho1] = epss;          // α₁ρ₁ ≈ 0 (空气几乎为0)
                            prim[AlphaRho2] = 5.04;          // α₂ρ₂ = ρ_SF6
                            prim[XMom]      = 0.0;           // u = 0
                            prim[YMom]      = 0.0;           // v = 0
                            prim[Energy]    = 1.0;           // p = 1.0
                            prim[Alpha1]    = epss;          // α₁ ≈ 0
                        } else if (distance > shock_radius + epss) {
                            // 区域 2: 已受冲击空气（激波后，靠近边缘）
                            prim[AlphaRho1] = 1.411;         // α₁ρ₁ = ρ_post_air
                            prim[AlphaRho2] = epss;          // α₂ρ₂ ≈ 0
                            // 激波向内传播，速度方向指向中心
                            prim[XMom]      = -0.39 * (x - sector_center_x) / distance; // 径向速度
                            prim[YMom]      = -0.39 * (y - sector_center_y) / distance; // 径向速度
                            prim[Energy]    = 1.628;         // p = 1.628
                            prim[Alpha1]    = 1.0 - epss;    // α₁ ≈ 1
                        } else {
                            // 区域 3: 未受冲击空气（激波前，介于界面与激波面之间）
                            prim[AlphaRho1] = 1.0;           // α₁ρ₁ = ρ_air
                            prim[AlphaRho2] = epss;          // α₂ρ₂ ≈ 0
                            prim[XMom]      = 0.0;           // u = 0
                            prim[YMom]      = 0.0;           // v = 0
                            prim[Energy]    = 1.0;           // p = 1.0
                            prim[Alpha1]    = 1.0 - epss;    // α₁ ≈ 1
                        }
                    } else {
                        // 区域：扇形外部（SF6）
                        prim[AlphaRho1] = epss;    // α1ρ1 (kg/m^3)
                        prim[AlphaRho2] = 5.04;    // α2ρ2 (kg/m^3)
                        prim[XMom] = 0.0;          // u (m/s)
                        prim[YMom] = 0.0;          // v (m/s)
                        prim[Energy] = 1.0;        // p = 1.0
                        prim[Alpha1] = epss;       // α1
                    }
                    break;
                }
                case Case_2D_Tin_Air_Implosion: {
                    // ============================================================================
                    // 2D算例：tin-air implosion
                    // ============================================================================
                    // 描述: 完整计算域上的锡-空气汇聚型内爆。
                    // 计算域: x in [-25, 25], y in [-25, 25]。
                    // 结构:
                    //   - r <= r_interface(theta): 空气腔体
                    //   - r_interface(theta) < r <= 24: 未受冲击的锡
                    //   - r > 24: 已受冲击的锡
                    // ============================================================================
                    constexpr int subcell_samples = 10;
                    constexpr Real cavity_base_radius = 20.0;
                    constexpr Real shock_radius = 24.0;
                    constexpr Real rho_tin_preshock = 7.28;
                    constexpr Real p_tin_preshock = 1.0;
                    constexpr Real rho_tin_postshock = 11.84;
                    constexpr Real p_tin_postshock = 1.0e6;
                    constexpr Real gamma_tin = 3.27;
                    constexpr Real p_inf_tin = 149500.0;
                    constexpr Real rho_air = 0.001;
                    constexpr Real p_air = 1.0;
                    constexpr Real gamma_air = 1.4;
                    constexpr Real p_inf_air = 0.0;
                    constexpr Real epss = 1.0e-8;
                    constexpr Real pi_over_2 = 0.5 * M_PI;

                    const Real e_tin_preshock = (p_tin_preshock + gamma_tin * p_inf_tin)
                        / ((gamma_tin - 1.0) * rho_tin_preshock);
                    const Real e_tin_postshock = (p_tin_postshock + gamma_tin * p_inf_tin)
                        / ((gamma_tin - 1.0) * rho_tin_postshock);
                    const Real e_air = (p_air + gamma_air * p_inf_air)
                        / ((gamma_air - 1.0) * rho_air);

                    const Real sample_dx = dx[0] / subcell_samples;
                    const Real sample_dy = dx[1] / subcell_samples;
                    const Real cell_x_lo = x - 0.5 * dx[0];
                    const Real cell_y_lo = y - 0.5 * dx[1];

                    int num_inside_shock = 0;
                    int num_inside_cavity = 0;

                    for (int sy = 0; sy < subcell_samples; ++sy) {
                        for (int sx = 0; sx < subcell_samples; ++sx) {
                            const Real sample_x = cell_x_lo + (sx + 0.5) * sample_dx;
                            const Real sample_y = cell_y_lo + (sy + 0.5) * sample_dy;
                            const Real radius = std::sqrt(sample_x * sample_x + sample_y * sample_y);

                            if (radius <= shock_radius) {
                                ++num_inside_shock;
                            }

                            const Real theta = std::atan2(sample_y, sample_x) - pi_over_2;
                            const Real r_interface = cavity_base_radius
                                + 0.4 * std::cos(22.0 * theta)
                                + 0.4 * std::cos(17.0 * theta)
                                + 0.3 * std::cos(29.0 * theta);

                            if (radius <= r_interface) {
                                ++num_inside_cavity;
                            }
                        }
                    }

                    const Real inv_total_samples = 1.0 / (subcell_samples * subcell_samples);
                    const Real shock_inside_fraction = num_inside_shock * inv_total_samples;

                    const Real rho1 = shock_inside_fraction * rho_tin_preshock
                        + (1.0 - shock_inside_fraction) * rho_tin_postshock;
                    const Real e1 = shock_inside_fraction * e_tin_preshock
                        + (1.0 - shock_inside_fraction) * e_tin_postshock;

                    Real alpha1 = 1.0 - num_inside_cavity * inv_total_samples;
                    if (alpha1 < epss) {
                        alpha1 = epss;
                    }
                    if (alpha1 > 1.0 - epss) {
                        alpha1 = 1.0 - epss;
                    }
                    const Real alpha2 = 1.0 - alpha1;
                    const Real rho2 = rho_air;

                    const Real total_internal_energy = alpha1 * rho1 * e1 + alpha2 * rho2 * e_air;
                    const Real gamma1_minus_1 = gamma_tin - 1.0;
                    const Real gamma2_minus_1 = gamma_air - 1.0;
                    const Real inv_gamma_minus_1 = alpha1 / gamma1_minus_1 + alpha2 / gamma2_minus_1;
                    const Real gamma_bar = 1.0 / inv_gamma_minus_1 + 1.0;
                    const Real gamma_pinf_bar = alpha1 * gamma_tin * p_inf_tin / gamma1_minus_1
                        + alpha2 * gamma_air * p_inf_air / gamma2_minus_1;
                    const Real p_inf_bar = gamma_pinf_bar * (gamma_bar - 1.0) / gamma_bar;
                    const Real pressure = (gamma_bar - 1.0) * total_internal_energy - gamma_bar * p_inf_bar;

                    prim[AlphaRho1] = alpha1 * rho1;
                    prim[AlphaRho2] = alpha2 * rho2;
                    prim[XMom] = 0.0;
                    prim[YMom] = 0.0;
                    prim[Energy] = pressure;
                    prim[Alpha1] = alpha1;
                    break;
                }

            }

            // 将原始变量转换为守恒变量
            VarArray cons = prim_to_cons(prim, pp);

            // 使用 Array4 将守恒变量存储到状态数组中
            state_arr(i, j, 0, AlphaRho1) = cons[AlphaRho1];
            state_arr(i, j, 0, AlphaRho2) = cons[AlphaRho2];
            state_arr(i, j, 0, XMom) = cons[XMom];
            state_arr(i, j, 0, YMom) = cons[YMom];
            state_arr(i, j, 0, Energy) = cons[Energy];
            state_arr(i, j, 0, Alpha1) = cons[Alpha1];
        }
    }
}
#endif

#if BL_SPACEDIM == 3
/**
 * 初始化三维两相算例的守恒量场
 * @param level 网格层级
 * @param time 当前时间
 * @param box 网格盒子
 * @param state 状态变量数组
 * @param dx 网格间距
 * @param problo 计算域左下角坐标
 * @param pp 两相物性参数
 */
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE void
initdata(int level, amrex::Real time,
    const amrex::Box& box, amrex::FArrayBox& state,
    const amrex::Real* dx, const amrex::Real* problo, const ProbParm& pp)
{
    using namespace amrex;

    // 原始变量数组
    VarArray prim;

    // 创建状态数组的 Array4 视图
    Array4<Real> state_arr = state.array();

    // 遍历网格框中的所有单元格（三维情况）
    const int i_lo = box.smallEnd(0);
    const int i_hi = box.bigEnd(0);
    const int j_lo = box.smallEnd(1);
    const int j_hi = box.bigEnd(1);
    const int k_lo = box.smallEnd(2);
    const int k_hi = box.bigEnd(2);

    for (int k = k_lo; k <= k_hi; ++k) {
        for (int j = j_lo; j <= j_hi; ++j) {
            for (int i = i_lo; i <= i_hi; ++i) {
                // 计算单元格中心的物理坐标
                Real x = problo[0] + (i + 0.5) * dx[0];
                Real y = problo[1] + (j + 0.5) * dx[1];
                Real z = problo[2] + (k + 0.5) * dx[2];

                // 根据选择的算例类型设置原始变量
                switch (pp.case_type) {
                    case Case_3D_RMI_Air_SF6_Perturbed: {
                        // ============================================================================
                        // 3D算例：三维Richtmyer–Meshkov不稳定性
                        // ============================================================================
                        // 编号: 3d_01
                        // 描述: 3D Richtmyer–Meshkov 不稳定性（空气–SF6）
                        // 区域划分:
                        //   - SF6区域: x < x_interface（重流体，左侧）
                        //     * 密度: 5.04 kg/m^3
                        //     * 速度: (0.0, 0.0, 0.0) m/s
                        //     * 压力: 1.0
                        //     * 体积分数: alpha1 = 0.0
                        //   - 激波后空气: x > x_shock（已受冲击空气，右侧）
                        //     * 密度: 1.411 kg/m^3
                        //     * 速度: (-0.39, 0.0, 0.0) m/s
                        //     * 压力: 1.628
                        //     * 体积分数: alpha1 = 1.0
                        //   - 激波前空气: 介于界面与激波面之间
                        //     * 密度: 1.0 kg/m^3
                        //     * 速度: (0.0, 0.0, 0.0) m/s
                        //     * 压力: 1.0
                        //     * 体积分数: alpha1 = 1.0
                        // 激波: 初始激波位于 x=3.2 处
                        // 界面: 扰动界面 x_interface = 2.9 - 0.2 * sin(2π(y + 0.25)) * sin(2π(z + 0.25))
                        // ============================================================================
                        constexpr Real x_shock = 3.2;          // 激波面位置 x_shock
                        constexpr Real x_1 = 2.9;          // 基准位置 x_1
                        constexpr Real varepsilon = 0.2;          // 扰动振幅 varepsilon
                        constexpr Real epss = 1e-8;            // 小量用于避免浮点误差
                        constexpr Real lamba = 1.0;            // 激波波速 lamba

                        // 计算扰动界面 x 坐标（作为 y 和 z 的函数，三维扰动）
                        // Real x_interface = x_1 - varepsilon * (std::sin(2.0 * M_PI * (y + 0.25)/lamba) * std::sin(2.0 * M_PI * (z + 0.25)/lamba));
                        Real x_interface = x_1 + varepsilon * (std::sin(2.0 * M_PI * (y + 0.25)/lamba) * std::sin(2.0 * M_PI * (z + 0.25)/lamba));

                        // 判断空间位置 → 三区域划分
                        if (x < x_interface) {
                            // 区域 1: SF6（重流体，左侧）
                            prim[AlphaRho1] = epss;          // α₁ρ₁ ≈ 0 (空气几乎为0)
                            prim[AlphaRho2] = 5.04;          // α₂ρ₂ = ρ_SF6
                            prim[XMom]      = 0.0;           // u = 0
                            prim[YMom]      = 0.0;           // v = 0
                            prim[ZMom]      = 0.0;           // w = 0
                            prim[Energy]    = 1.0;           // p = 1.0
                            prim[Alpha1]    = epss;          // α₁ ≈ 0
                        } else if (x > x_shock + epss) {
                            // 区域 2: 已受冲击空气（激波后，右侧）
                            prim[AlphaRho1] = 1.411;         // α₁ρ₁ = ρ_post_air
                            prim[AlphaRho2] = epss;          // α₂ρ₂ ≈ 0
                            prim[XMom]      = -0.39;         // u = -0.39（向左）
                            prim[YMom]      = 0.0;           // v = 0
                            prim[ZMom]      = 0.0;           // w = 0
                            prim[Energy]    = 1.628;         // p = 1.628
                            prim[Alpha1]    = 1.0 - epss;    // α₁ ≈ 1
                        } else {
                            // 区域 3: 未受冲击空气（激波前，介于界面与激波面之间）
                            prim[AlphaRho1] = 1.0;           // α₁ρ₁ = ρ_air
                            prim[AlphaRho2] = epss;          // α₂ρ₂ ≈ 0
                            prim[XMom]      = 0.0;           // u = 0
                            prim[YMom]      = 0.0;           // v = 0
                            prim[ZMom]      = 0.0;           // w = 0
                            prim[Energy]    = 1.0;           // p = 1.0
                            prim[Alpha1]    = 1.0 - epss;    // α₁ ≈ 1
                        }
                        break;
                    }

                }

                // 将原始变量转换为守恒变量
                VarArray cons = prim_to_cons(prim, pp);

                // 使用 Array4 将守恒变量存储到状态数组中
                state_arr(i, j, k, AlphaRho1) = cons[AlphaRho1];
                state_arr(i, j, k, AlphaRho2) = cons[AlphaRho2];
                state_arr(i, j, k, XMom) = cons[XMom];
                state_arr(i, j, k, YMom) = cons[YMom];
                state_arr(i, j, k, ZMom) = cons[ZMom];
                state_arr(i, j, k, Energy) = cons[Energy];
                state_arr(i, j, k, Alpha1) = cons[Alpha1];
            }
        }
    }
}
#endif
