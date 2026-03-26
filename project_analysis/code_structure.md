# 代码结构详解

## 根目录文件

| 文件名 | 功能描述 | 重要性 |
|--------|---------|--------|
| main.cpp | 主程序入口，包含AMR初始化和主循环 | 核心 |
| CmakeLists.txt | CMake构建配置文件 | 构建 |
| Make.package | AMReX包配置文件 | 构建 |
| testEigenSystem.cpp | 特征系统测试 | 测试 |
| testEquation.cpp | 方程测试 | 测试 |

## 2phase 目录

| 文件名 | 功能描述 | 重要性 |
|--------|---------|--------|
| equation.H | 状态方程和物理模型 | 核心 |
| equation.cpp | 方程实现 | 核心 |
| eigen_system.cpp | 特征系统计算 | 核心 |
| deriv_usr.H | HLLC Riemann 求解器 | 核心 |
| variable_set.cpp | 变量定义和转换 | 核心 |
| bc_nullfill.cpp | 边界条件处理 | 重要 |
| physical_bound.cpp | 物理边界条件 | 重要 |
| erroEst.cpp | 误差估计（AMR） | 重要 |
| Prob_Parm.H | 问题参数定义 | 重要 |

## amrex_level_Cong 目录

| 文件名 | 功能描述 | 重要性 |
|--------|---------|--------|
| AmrLevelCong.H | AMR级别类声明 | 核心 |
| AmrLevelCong.cpp | 类实现 | 核心 |
| AmrLevelCong_advance.cpp | 时间推进实现 | 核心 |
| LevelBldCong.cpp | 级别构建器 | 重要 |
| Tagging_params.cpp | 网格加密参数 | 重要 |
| output_tmp.H | 输出处理 | 重要 |
| Adv_prob.cpp | 问题初始化 | 重要 |

## deriv 目录

| 文件名 | 功能描述 | 重要性 |
|--------|---------|--------|
| deriv.H | 一阶导数项模板 | 核心 |

## reconstruction_scheme 目录

| 文件名 | 功能描述 | 重要性 |
|--------|---------|--------|
| positive_preserving.cpp | 保正重构算法 | 核心 |

## data 目录

| 文件名 | 功能描述 | 重要性 |
|--------|---------|--------|
| inputs2phase1d | 1D双相流配置 | 重要 |
| inputs2phase2d | 2D双相流配置 | 重要 |
| inputs2phase2dMach10 | 高速流动配置 | 重要 |
| euler_data_level_3.dat | 欧拉方程数据 | 数据 |

## 核心类和函数

### AmrLevelCong 类
- **功能**: 处理单个AMR级别的计算
- **主要方法**:
  - `advance()`: 时间推进
  - `errorEst()`: 误差估计
  - `initData()`: 初始化数据
  - `compute_dSdt()`: 计算时间导数

### FirstOrderDerivativeTerm 模板
- **功能**: 处理一阶导数项
- **主要方法**:
  - `computeCoefficient()`: 计算系数
  - `computeFlux()`: 计算通量
  - `computeNodeFlux()`: 计算节点通量

### HLLC 求解器
- **功能**: 计算界面通量
- **核心算法**:
  1. 计算左右状态
  2. 估计波速 SL, SR
  3. 计算接触间断速度 S*
  4. 计算界面通量

### 状态方程 (EOS)
- **功能**: 计算物理状态
- **支持的计算**:
  - 压力
  - 声速
  - 总能量
  - 声速平方×密度

## 代码执行流程

1. **初始化阶段**:
   - 读取配置文件
   - 建立AMR网格
   - 初始化物理参数

2. **时间推进阶段**:
   - 计算时间步长 dt
   - 执行Runge-Kutta时间积分
   - 重构状态变量
   - 计算通量
   - 更新保守变量

3. **AMR阶段**:
   - 误差估计
   - 网格加密/粗化
   - 数据传递

4. **输出阶段**:
   - 写入checkpoint文件
   - 写入plotfile文件
   - 导出数据到文件

## 依赖关系

- **AMReX**:
  - 网格管理
  - 并行计算
  - I/O操作
  - 数据结构

- **MPI**:
  - 分布式计算
  - 进程间通信

## 编译选项

- **CMAKE_BUILD_TYPE**: Debug/Release
- **AMReX_SPACEDIM**: 1/2/3 (空间维度)
- **优化选项**: -O3
- **C++标准**: C++23

## 代码特点

1. **模板元编程**: 利用C++模板实现通用算法
2. **GPU支持**: 利用AMReX的GPU加速能力
3. **并行计算**: 支持多核心并行
4. **模块化设计**: 清晰的代码结构
5. **性能优化**: 内存访问优化，计算密集部分的优化