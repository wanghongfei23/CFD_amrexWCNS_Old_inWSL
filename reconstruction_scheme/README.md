# 重构方案（Reconstruction Scheme）目录介绍

本文档介绍 `reconstruction_scheme` 目录中的插值格式文件及其功能。这些文件实现了不同阶数和类型的数值重构方法，用于计算流体力学中的界面状态。

## 目录结构

```
reconstruction_scheme/
├── 1th_01_recon.H         # 一阶 Godunov 型重构
├── 2th_01_muscl_recon.H   # MUSCL 限幅重构
├── 5th_01_recon.H         # 基于特征分解的五阶重构
├── 5th_02_thinc_recon.H   # 五阶 THINC/WENO 混合重构
├── eigen_system.H          # 特征矩阵与特征变量变换
├── interpolation_schemes.H # WENO、TENO 等插值方案
├── linear_algebra.H        # 线性代数与方向遍历工具
├── positive_preserving.H   # 保正限制器
└── teno_related.H          # TENO 相关功能
```

## 重构方法介绍

### 1. 一阶重构（1th_01_recon.H）

**功能**：实现一阶 Godunov 型重构，是最简单的重构方法。

**特点**：
- 直接使用单元格中心值作为界面状态
- 计算简单，耗散较大
- 适用于初始测试或作为基准方法

**核心函数**：
```cpp
template <int dir = 0>
amrex::GpuArray<VarArray, 2> reconstruction(int i, int j, int k, ...)
```

### 2. MUSCL 限幅重构（2th_01_muscl_recon.H）

**功能**：实现二阶 MUSCL（Monotonic Upwind Scheme for Conservation Laws）重构。

**特点**：
- 使用 minmod 限制器确保单调性
- 二阶精度，比一阶重构更准确
- 计算效率较高

**核心函数**：
- `musclInterpolation`：MUSCL 插值计算
- `reconstruction`：重构主函数

### 3. 基于特征分解的五阶重构（5th_01_recon.H）

**功能**：实现基于特征空间的五阶 WENO 重构。

**特点**：
- 在特征空间进行重构，提高激波捕捉能力
- 五阶精度，精度高
- 包含保正限制器，确保物理量非负

**核心函数**：
```cpp
template <int dir = 0>
amrex::GpuArray<VarArray, 2> reconstruction(int i, int j, int k, ...)
```

### 4. 五阶 THINC/WENO 混合重构（5th_02_thinc_recon.H）

**功能**：实现 THINC（Tangent of Hyperbola for INterface Capturing）与 WENO 混合的五阶重构。

**特点**：
- 结合 THINC 界面捕捉和 WENO 高精度特性
- 特别适合多相流模拟中的界面处理
- 包含多种重构策略（原始变量重构、特征变量重构、混合重构）

**核心函数**：
- `phi_THINC`：THINC 界面重构
- `primitive_reconstruction`：原始变量重构
- `characteristic_reconstruction`：特征变量重构
- `mixture_reconstruction`：混合重构
- `reconstruction`：重构主函数

## 辅助模块

### 1. 特征系统（eigen_system.H）

**功能**：实现特征矩阵与特征变量变换。

**核心类**：
- `EigenSystem`：在 Roe 平均状态上构建左右特征矩阵

**主要方法**：
- `prim_to_char`：原始变量到特征变量的转换
- `char_to_prim`：特征变量到原始变量的转换
- `char_to_cons`：特征变量到守恒变量的转换

### 2. 插值方案（interpolation_schemes.H）

**功能**：提供各种高级插值方案。

**主要函数**：
- `weno5_Z`：五阶 WENO-Z 插值
- `weno5_JSchen`：另一种 WENO 实现
- `weno_is`：改进的 WENO 插值
- `WENO3`：三阶 WENO 插值
- `Teno5_CongZ`：五阶 TENO 插值

### 3. 线性代数（linear_algebra.H）

**功能**：提供轻量级线性代数工具和方向遍历功能。

**主要函数**：
- `traverse_in_direction`：方向遍历
- `traverse_in_direction_with_ii`：带索引的方向遍历
- `LA_M_X`：矩阵-向量乘法
- `compute_normal_vector`：法向量计算
- `LinComb3`/`LinComb4`：线性组合操作

### 4. 保正限制器（positive_preserving.H）

**功能**：确保物理量的非负性。

**主要功能**：
- 提供 `apply_positive_limiter` 函数，确保密度、体积分数等物理量非负

## 使用方法

### 选择重构方法

在代码中，可以通过包含不同的头文件来选择不同的重构方法：

```cpp
// 使用一阶重构
#include "reconstruction_scheme/1th_01_recon.H"

// 或使用 MUSCL 重构
#include "reconstruction_scheme/2th_01_muscl_recon.H"

// 或使用五阶特征重构
#include "reconstruction_scheme/5th_01_recon.H"

// 或使用五阶 THINC/WENO 混合重构
#include "reconstruction_scheme/5th_02_thinc_recon.H"
```

### 调用重构函数

所有重构方法都提供了相同的 `reconstruction` 模板函数接口：

```cpp
// 沿 x 方向（dir=0）重构
auto [left_state, right_state] = reconstruction<0>(i, j, k, primin_array, NUM_STATE, pp, norm, dx);

// 沿 y 方向（dir=1）重构
auto [left_state, right_state] = reconstruction<1>(i, j, k, primin_array, NUM_STATE, pp, norm, dx);

// 沿 z 方向（dir=2）重构
auto [left_state, right_state] = reconstruction<2>(i, j, k, primin_array, NUM_STATE, pp, norm, dx);
```

## 性能考虑

- **一阶重构**：计算最快，但精度最低
- **MUSCL 重构**：平衡了计算速度和精度
- **五阶重构**：精度最高，但计算成本也最高
- **THINC/WENO 混合重构**：在多相流模拟中特别有效，但计算复杂度较高

## 适用场景

| 重构方法 | 适用场景 | 优势 |
|---------|---------|------|
| 一阶重构 | 初始测试、简单流动 | 计算快，实现简单 |
| MUSCL 重构 | 一般流体模拟 | 平衡速度和精度 |
| 五阶特征重构 | 高精度模拟、激波捕捉 | 精度高，激波处理好 |
| THINC/WENO 混合重构 | 多相流模拟、界面捕捉 | 界面处理效果好 |

## 扩展与定制

如果需要扩展或定制重构方法，可以：

1. 在 `interpolation_schemes.H` 中添加新的插值方案
2. 在 `eigen_system.H` 中修改特征系统实现
3. 创建新的重构文件，实现特定需求的重构方法

## 注意事项

- 不同重构方法的计算成本差异较大，应根据具体问题选择合适的方法
- 高阶重构在边界附近需要特殊处理，确保精度
- 在多相流模拟中，界面捕捉的质量对结果影响很大，建议使用 THINC/WENO 混合重构
