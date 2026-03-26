# 快速上手指南

## 环境准备

### 1. 系统要求

- **操作系统**: Linux (WSL2 或原生 Linux)
- **编译器**: GCC 9.0+ 或 Clang
- **依赖**: 
  - AMReX 框架
  - MPI (OpenMPI 或 MPICH)
  - CMake 3.12+

### 2. AMReX 安装

```bash
# 克隆 AMReX 仓库
git clone https://github.com/AMReX-Codes/amrex.git

# 编译 AMReX
cd amrex
mkdir build && cd build
cmake .. -DAMReX_SPACEDIM=2 -DAMReX_MPI=ON
make -j4

# 安装 (可选)
sudo make install
```

### 3. 项目配置

编辑 `CmakeLists.txt` 文件，设置 AMReX 路径:

```cmake
# 设置 AMReX 安装路径
set(AMReX_DIR /opt/amrex_2d/lib/cmake/AMReX)
```

## 编译项目

```bash
# 进入项目目录
cd /home/archwanghongfei/Documents/GitHubWSL/CFD_inWSL_MF_initial

# 创建构建目录
mkdir -p build && cd build

# 配置 CMake
cmake ..

# 编译
make -j4
```

## 运行示例

### 1. 2D 双相流算例

```bash
# 运行 2D 双相流算例
./twoPhaseSolver ../data/inputs2phase2d
```

### 2. 1D 激波管问题

```bash
# 运行 1D 激波管算例
./twoPhaseSolver ../data/inputs2phase1d
```

### 3. 高速流动算例

```bash
# 运行高速流动算例
./twoPhaseSolver ../data/inputs2phase2dMach10
```

## 配置文件说明

### inputs2phase2d 主要参数

```bash
# 网格设置
amr.n_cell = 64 64       # 网格数量
amr.max_level = 2         # 最大网格级别
amr.ref_ratio = 2         # 网格加密比例

# 时间控制
amr.max_step = 100        # 最大时间步数
amr.stop_time = 0.1       # 停止时间

# 输出控制
amr.plot_file = plt       # 输出前缀
amr.plot_int = 10         # 输出间隔
amr.check_file = chk      # 检查点前缀
amr.check_int = 50        # 检查点间隔

# 物理参数
amr.gamma = 1.4 1.4       # 两种流体的比热比
amr.pInf = 0.0 0.0        # 两种流体的参考压力
```

## 结果查看

### 1. 输出文件

- **plt00000/**: 包含可视化数据
- **chk00000/**: 包含检查点数据
- **euler_data_level_3.dat**: 导出的数据文件

### 2. 使用 ParaView 查看

```bash
# 打开 ParaView
paraview

# 选择 File → Open → plt00000
# 选择 AMReX PlotFile
```

### 3. 数据导出

项目会自动导出数据到 `euler_data_level_3.dat` 文件，可以用 Python 或 MATLAB 分析:

```python
import numpy as np
import matplotlib.pyplot as plt

# 读取数据
data = np.loadtxt('euler_data_level_3.dat')

# 绘制密度分布
plt.plot(data[:, 0], data[:, 1])
plt.xlabel('X')
plt.ylabel('Density')
plt.show()
```

## 常见问题

### 1. 编译错误

**问题**: 找不到 AMReX
**解决**: 检查 `AMReX_DIR` 是否正确设置

**问题**: 缺少 MPI
**解决**: 安装 OpenMPI 或 MPICH

### 2. 运行错误

**问题**: NaN 检测到
**解决**: 检查初始条件，调整 CFL 数

**问题**: 内存不足
**解决**: 减少网格数量或最大级别

**问题**: 并行错误
**解决**: 检查 MPI 配置，尝试使用 `mpirun -np 1` 运行

### 3. 结果异常

**问题**: 相界面不稳定
**解决**: 增加网格分辨率，调整数值参数

**问题**: 激波捕捉不准确
**解决**: 检查 Riemann 求解器设置

## 调优建议

### 1. 性能优化

- **网格设置**: 根据问题复杂度调整 `amr.n_cell` 和 `amr.max_level`
- **并行计算**: 使用 `mpirun -np N` 进行并行计算
- **编译优化**: 使用 `-O3` 优化级别

### 2. 数值参数调整

- **CFL 数**: 默认为 0.8，可根据稳定性调整
- **重构格式**: 可在 `reconstruction_scheme` 中选择不同的重构方法
- **边界条件**: 根据具体问题调整边界条件

### 3. 内存管理

- **输出频率**: 减少 `amr.plot_int` 以减少 I/O 操作
- **网格加密**: 合理设置 `amr.ref_ratio` 和 `amr.max_level`
- **数据类型**: 检查是否使用了合适的数据类型

## 高级用法

### 1. 自定义初始条件

修改 `Adv_prob.cpp` 文件中的初始条件设置:

```cpp
// 设置初始条件
for (int i = 0; i < ncells; i++) {
    if (x[i] < 0.5) {
        // 左侧状态
        data[i][AlphaRho1] = 1.0;
        data[i][AlphaRho2] = 0.0;
        data[i][XMOM] = 0.0;
        data[i][YMOM] = 0.0;
        data[i][Energy] = 1.0;
        data[i][Alpha1] = 1.0;
    } else {
        // 右侧状态
        data[i][AlphaRho1] = 0.0;
        data[i][AlphaRho2] = 1.0;
        data[i][XMOM] = 0.0;
        data[i][YMOM] = 0.0;
        data[i][Energy] = 0.1;
        data[i][Alpha1] = 0.0;
    }
}
```

### 2. 自定义物理模型

修改 `equation.H` 文件中的状态方程:

```cpp
// 自定义状态方程
template <EOSType eType>
amrex::Real eos_cons(VarArray cons, ProbParm const& pp) {
    // 实现自定义状态方程
}
```

### 3. 自定义数值格式

修改 `deriv_usr.H` 文件中的通量计算:

```cpp
// 自定义通量计算
amrex::GpuArray<amrex::Real, num_terms>
computeFluxImpl(const VarArray& stateL, const VarArray& stateR, ...) {
    // 实现自定义通量计算
}
```

## 故障排除

### 1. 检查 AMReX 版本

```bash
# 检查 AMReX 版本
amrex-config --version
```

### 2. 检查编译配置

```bash
# 检查 CMake 配置
cd build
cmake .. -LA
```

### 3. 运行调试版本

```bash
# 编译调试版本
cmake .. -DCMAKE_BUILD_TYPE=Debug
make -j4

# 运行调试版本
gdb --args ./twoPhaseSolver ../data/inputs2phase2d
```

### 4. 查看日志文件

检查运行时生成的日志文件，了解详细的错误信息。

## 下一步

1. **熟悉代码结构**: 查看 `code_structure.md` 了解代码组织
2. **学习可视化**: 参考 `visualization_guide.md` 学习如何可视化结果
3. **尝试不同算例**: 运行不同的输入文件，了解不同问题的模拟
4. **扩展功能**: 根据需要修改代码，添加新的物理模型或数值方法