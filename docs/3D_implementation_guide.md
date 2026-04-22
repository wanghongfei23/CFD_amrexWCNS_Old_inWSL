# 三维算例实现指南

## 项目结构分析

本项目是一个基于AMReX框架的两相流求解器，主要用于模拟多相流体动力学问题。当前配置为二维模拟，我们需要将其扩展到三维。

### 主要目录结构

- **2phase/**: 两相流求解器核心代码

  - `equation.H`: 状态方程、变量变换、Riemann求解与初值设置
  - `variable_set.cpp`: 变量定义和边界条件设置
  - `physical_bound.cpp`: 物理边界条件
  - `eigen_system.cpp`: 特征系统
  - `bc_nullfill.cpp`: 边界填充
  - `erroEst.cpp`: 误差估计
- **amrex_level_Cong/**: AMReX框架相关代码

  - `AmrLevelCong.cpp`: AMR级别实现
  - `AmrLevelCong_advance.cpp`: 时间推进
  - `Adv_prob.cpp`: 问题设置
- **reconstruction_scheme/**: 数值重构方案

  - `first_order_recon.H`: 一阶重构
  - `muscl_recon.H`: MUSCL重构
  - `fifth_order_recon.H`: 五阶重构
- **data/**: 输入参数文件

  - 包含各种二维算例的输入文件

### 当前配置

- 使用2D版本的AMReX库
- 代码中已有部分3D相关的宏定义和结构
- 缺少完整的3D初始条件设置和输入文件

## 实现三维算例的详细步骤

### 1. 修改编译配置

**目标**: 将AMReX库从2D版本切换到3D版本

**操作步骤**:

1. 打开 `CMakeLists.txt` 文件
2. 修改AMReX库路径，注释2D版本，启用3D版本:

```cmake
# 设置 AMReX 安装路径
# set(AMReX_DIR /opt/amrex_1d/lib/cmake/AMReX)
# set(AMReX_DIR /opt/amrex_2d/lib/cmake/AMReX)
set(AMReX_DIR /opt/amrex_3d/lib/cmake/AMReX)
```

### 2. 实现3D初始条件

**目标**: 在 `2phase/equation.H` 中添加3D版本的 `initdata` 函数

**操作步骤**:

1. 打开 `2phase/equation.H` 文件
2. 在文件末尾添加3D初始条件函数:

```cpp
#if BL_SPACEDIM == 3
/**
 * @brief 初始化三维算例的守恒量场。
 * @details 提供多种三维两相流算例配置。
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

                // ========================================================================
                //【算例】1：三维水柱配置
                // 描述：简单的水柱和空气区域划分
                // ========================================================================
                constexpr Real water_cylinder_x_center = 0.5; // 水柱中心 x 坐标 [m]
                constexpr Real water_cylinder_y_center = 0.5; // 水柱中心 y 坐标 [m]
                constexpr Real water_cylinder_z_center = 0.5; // 水柱中心 z 坐标 [m]
                constexpr Real water_cylinder_radius = 0.2; // 水柱半径 [m]
                constexpr Real epss = 1e-8;
              
                Real distance = std::sqrt(
                    (x - water_cylinder_x_center) * (x - water_cylinder_x_center) + 
                    (y - water_cylinder_y_center) * (y - water_cylinder_y_center) +
                    (z - water_cylinder_z_center) * (z - water_cylinder_z_center)
                );
                                  
                if (distance <= water_cylinder_radius) {
                    // 区域：水柱内部
                    prim[AlphaRho1] = 1000.0; // α1ρ1 (kg/m^3)
                    prim[AlphaRho2] = epss; // α2ρ2 (kg/m^3)
                    prim[XMom] = 0.0;         // u (m/s)
                    prim[YMom] = 0.0;         // v (m/s)
                    prim[ZMom] = 0.0;         // w (m/s)
                    prim[Energy] = 1.01e5;     // p (Pa)
                    prim[Alpha1] = 1.0 - epss; // α1
                } else { 
                    // 区域：水柱外部（空气）
                    prim[AlphaRho1] = epss;  // α1ρ1 (kg/m^3)
                    prim[AlphaRho2] = 1.17;    // α2ρ2 (kg/m^3)
                    prim[XMom] = 0.0;          // u (m/s)
                    prim[YMom] = 0.0;          // v (m/s)
                    prim[ZMom] = 0.0;          // w (m/s)
                    prim[Energy] = 1.01e5;      // p (Pa)
                    prim[Alpha1] = epss;       // α1
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
```

### 3. 验证3D数值方法

**目标**: 确保所有数值函数正确处理3D情况

**检查要点**:

1. **变量变换函数** (`cons_to_prim`, `prim_to_cons`):

   - 已使用 `AMREX_D_DECL` 宏，应自动支持3D
   - 检查是否正确处理Z方向动量
2. **通量计算** (`prim_to_flux`):

   - 已使用 `AMREX_D_TERM` 宏，应自动支持3D
   - 检查法向速度计算是否包含Z分量
3. **Riemann求解器** (`Riemann_solver_1`, `Riemann_solver_2`):

   - 已使用 `AMREX_D_DECL` 宏，应自动支持3D
   - 检查是否正确处理Z方向速度
4. **界面法向量计算** (`compute_normal_vector1`):

   - 已包含Z方向梯度计算
   - 检查是否正确归一化3D向量

### 4. 创建3D输入文件

**目标**: 在 `data/` 目录下创建3D输入文件

**操作步骤**:

1. 复制现有的2D输入文件，例如 `data/inputs2phase2d`
2. 修改为3D配置，保存为 `data/inputs2phase3d`:

```
# 三维两相流算例输入文件

amrex.spacedim = 3

# 网格设置
amrex.geometry.prob_lo = 0.0 0.0 0.0
amrex.geometry.prob_hi = 1.0 1.0 1.0

amrex.amr.n_cell = 64 64 64
amrex.amr.max_level = 0

# 时间步进
amrex.integration.cfl = 0.3
amrex.integration.stop_time = 0.1
amrex.integration.dt = 1e-6

# 输出设置
amrex.amr.plot_file = plt
amrex.amr.plot_int = 10

# 边界条件
amrex.density1_lo_bc_0 = 1
amrex.density1_hi_bc_0 = 1
amrex.density1_lo_bc_1 = 1
amrex.density1_hi_bc_1 = 1
amrex.density1_lo_bc_2 = 1
amrex.density1_hi_bc_2 = 1

amrex.density2_lo_bc_0 = 1
amrex.density2_hi_bc_0 = 1
amrex.density2_lo_bc_1 = 1
amrex.density2_hi_bc_1 = 1
amrex.density2_lo_bc_2 = 1
amrex.density2_hi_bc_2 = 1

amrex.energy_lo_bc_0 = 1
amrex.energy_hi_bc_0 = 1
amrex.energy_lo_bc_1 = 1
amrex.energy_hi_bc_1 = 1
amrex.energy_lo_bc_2 = 1
amrex.energy_hi_bc_2 = 1

amrex.VolumeFraction1_lo_bc_0 = 1
amrex.VolumeFraction1_hi_bc_0 = 1
amrex.VolumeFraction1_lo_bc_1 = 1
amrex.VolumeFraction1_hi_bc_1 = 1
amrex.VolumeFraction1_lo_bc_2 = 1
amrex.VolumeFraction1_hi_bc_2 = 1

amrex.Momentum_X_lo_bc_0 = 1
amrex.Momentum_X_hi_bc_0 = 1
amrex.Momentum_X_lo_bc_1 = 1
amrex.Momentum_X_hi_bc_1 = 1
amrex.Momentum_X_lo_bc_2 = 1
amrex.Momentum_X_hi_bc_2 = 1

amrex.Momentum_Y_lo_bc_0 = 1
amrex.Momentum_Y_hi_bc_0 = 1
amrex.Momentum_Y_lo_bc_1 = 1
amrex.Momentum_Y_hi_bc_1 = 1
amrex.Momentum_Y_lo_bc_2 = 1
amrex.Momentum_Y_hi_bc_2 = 1

amrex.Momentum_Z_lo_bc_0 = 1
amrex.Momentum_Z_hi_bc_0 = 1
amrex.Momentum_Z_lo_bc_1 = 1
amrex.Momentum_Z_hi_bc_1 = 1
amrex.Momentum_Z_lo_bc_2 = 1
amrex.Momentum_Z_hi_bc_2 = 1

# 物理参数
prob.gamma1 = 1.4
prob.pInf1 = 0.0
prob.gamma2 = 1.0
prob.pInf2 = 0.0
```

### 5. 编译和运行

**目标**: 编译3D版本并运行算例

**操作步骤**:

1. 清理并重新编译项目:

cd build && rm -rf * && cmake .. && make -j4

```bash
cd build
rm -rf *
cmake ..
make -j4
```

2. 运行3D算例:

```bash
./twoPhaseSolver ../data/inputs2phase3d
```

### 6. 结果可视化

**目标**: 查看3D模拟结果

**操作步骤**:

1. 使用ParaView打开输出的plot文件:

```bash
paraview plt00000/
```

2. 在ParaView中设置合适的可视化参数，如等值面、切片等，查看3D流场结构

## 常见问题和解决方案

### 1. 编译错误

**问题**: 找不到3D版本的AMReX库

**解决方案**:

- 确保3D版本的AMReX库已正确安装
- 检查CMakeLists.txt中的路径是否正确

### 2. 运行时错误

**问题**: 数组越界或内存访问错误

**解决方案**:

- 检查3D初始条件函数中的数组访问
- 确保所有3D相关的索引都正确处理

### 3. 物理结果不正确

**问题**: 3D模拟结果与预期不符

**解决方案**:

- 检查3D初始条件设置
- 验证Riemann求解器在3D情况下的正确性
- 确保边界条件在3D情况下正确应用

### 4. 性能问题

**问题**: 3D模拟运行速度过慢

**解决方案**:

- 调整网格分辨率，从粗网格开始测试
- 考虑使用并行计算
- 优化数值算法，减少计算量

## 扩展建议

1. **实现更多3D算例**:

   - 3D Richtmyer-Meshkov不稳定性
   - 3D激波-气泡相互作用
   - 3D多相湍流
2. **添加3D可视化工具**:

   - 开发专门的3D结果后处理脚本
   - 集成ParaView的Python接口
3. **性能优化**:

   - 利用GPU加速3D计算
   - 实现自适应网格加密(AMR)策略
   - 优化内存访问模式
4. **物理模型扩展**:

   - 添加表面张力模型
   - 实现多组分模型
   - 考虑化学反应

## 总结

将两相流求解器从2D扩展到3D需要以下关键步骤:

1. **切换到3D AMReX库**
2. **实现3D初始条件**
3. **验证3D数值方法**
4. **创建3D输入文件**
5. **编译和运行3D算例**
6. **可视化和分析结果**

项目代码已经使用了AMReX的维度无关宏，因此大部分代码可以自动适应3D情况。主要工作是添加3D特定的初始条件和确保编译配置正确。

通过以上步骤，您应该能够成功实现三维算例的数值模拟，并探索更复杂的多相流现象。
