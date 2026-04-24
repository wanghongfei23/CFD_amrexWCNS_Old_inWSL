# CFD_inWSL_MF_initial

基于 AMReX 的多相可压缩流动求解器实验仓库。当前主线构建目标是 `twoPhaseSolver`，代码入口在 `main.cpp`，求解器主体由 `2phase/`、`amrex_level_Cong/`、`reconstruction_scheme/` 和 `deriv/` 组成。

当前仓库更适合按“正在使用的源码快照”来理解，而不是按历史设计目标来理解。以当前 `CMakeLists.txt` 为准，默认构建会链接 `/opt/amrex_2d/lib/cmake/AMReX`，因此默认产物是二维版本的 `twoPhaseSolver`。

## 当前状态

- 可执行文件：`twoPhaseSolver`
- 主程序入口：`main.cpp`
- 默认 AMReX 维度：2D
- 构建系统：CMake
- 并行方式：MPI
- 运行输入：`data/` 目录下的 AMReX inputs 文件

需要注意的是，仓库里虽然保留了 1D、2D 和 3D 的输入文件与部分实现痕迹，但当前 `CMakeLists.txt` 默认启用的是双相求解器主线，并且默认链接的是二维 AMReX。想跑 1D 或 3D，通常需要先切换 `AMReX_DIR` 再重新编译。

## 依赖环境

建议运行环境：

- Linux 或 WSL
- CMake 3.12 及以上
- 支持 C++23 的编译器
- MPI
- AMReX

当前工程要求 `AMReX_DIR` 指向已安装的 AMReX CMake 配置目录。默认配置如下：

```cmake
set(AMReX_DIR /opt/amrex_2d/lib/cmake/AMReX)
```

如果你要切换到 1D 或 3D，请先修改根目录 `CMakeLists.txt` 中的 `AMReX_DIR`，再重新配置和编译。

## 编译方法

在项目根目录执行：

```bash
cmake -S . -B build
cmake --build build -j
```

编译成功后，可执行文件通常位于：

```bash
build/twoPhaseSolver
```

如果你修改了以下内容之一，需要重新编译：

- 求解器实现文件，例如 `2phase/*.cpp`、`amrex_level_Cong/*.cpp`
- 头文件中的核心逻辑，例如 `2phase/equation.H`、`2phase/Prob_Parm.H`
- 重构与通量相关实现，例如 `reconstruction_scheme/` 下文件

如果你只修改 `data/` 下的输入文件，一般不需要重新编译。

## 运行方法

建议从 `build/` 目录启动，这样输出文件会落在构建目录，不会把仓库根目录弄乱。

串行运行示例：

```bash
cd build
./twoPhaseSolver ../data/inputs2phase2d
```

MPI 运行示例：

```bash
cd build
mpiexec -n 8 ./twoPhaseSolver ../data/inputs2phase2d
```

如果只是做快速检查，可以直接在命令行覆盖 inputs 参数：

```bash
./twoPhaseSolver ../data/inputs2phase2d max_step=10
./twoPhaseSolver ../data/inputs2phase2d max_step=10 amr.plot_files_output=0 amr.checkpoint_files_output=0
```

## 输入文件说明

输入文件位于 `data/` 目录，例如：

- `data/inputs2phase2d`
- `data/inputs2phase2dMach10`
- `data/inputs2phase3dRMI_SF6_2half`

这些文件负责控制：

- 时间推进参数，例如 `max_step`、`stop_time`、`adv.cfl`
- 几何区域与网格，例如 `geometry.prob_lo`、`geometry.prob_hi`、`amr.n_cell`
- AMR 参数，例如 `amr.max_level`、`amr.ref_ratio`、`amr.regrid_int`
- 边界条件映射相关参数

但当前双相算例并不是“完全由 inputs 文件定义”。下面两点尤其重要：

1. 初始条件有相当一部分硬编码在 `2phase/equation.H` 的 `initdata` 相关逻辑中。
2. 部分物性参数默认值定义在 `2phase/Prob_Parm.H` 中，而不是由 inputs 完整驱动。

如果你改了 inputs 却发现结果几乎不变，优先检查这两个位置。

## 输出结果

程序运行后通常会产生两类输出：

1. AMReX 标准输出，例如 plotfile 和 checkpoint。
2. 项目自定义 ASCII 输出，例如 `euler_data_level_0.dat`。

根据当前源码，自定义 ASCII 输出由 `amrex_level_Cong/output_tmp.H` 导出。即使关闭 plotfile 或 checkpoint，这类 ASCII 文件仍然可能继续生成。

## 目录结构

项目根目录下常用部分如下：

- `2phase/`：双相模型、边界、特征系统、变量定义等核心实现
- `amrex_level_Cong/`：AMReX 层级推进、参数读取、时间推进和输出逻辑
- `reconstruction_scheme/`：重构方法、正性保持及相关线性代数实现
- `deriv/`：导数与派生量相关实现
- `data/`：算例输入文件
- `docs/`：当前源码行为说明、3D 实现说明等补充文档
- `tests/`：测试目录，目前内容较少
- `build/`：本地构建输出目录

仓库中还保留了 `Euler1D/`、`Euler2D/` 等目录，它们更像历史分支或备用实现，不是当前 `twoPhaseSolver` 默认构建路径的一部分。

## 常见修改入口

如果你的目标是调整算例，通常从下面几个位置入手：

- 改网格、终止条件、AMR 参数：修改 `data/` 下的 inputs 文件
- 改双相初始条件：修改 `2phase/equation.H`
- 改物性默认参数：修改 `2phase/Prob_Parm.H`
- 改 AMR 推进与读参逻辑：修改 `amrex_level_Cong/`
- 改重构或正性保持策略：修改 `reconstruction_scheme/`

## 相关文档

如果需要更细的源码级说明，可以继续看：

- `docs/current_source_usage_guide.md`
- `docs/3D_implementation_guide.md`
- `docs/time_advancement_function_call_chain.md`
- `docs/using.md`

## 快速开始

如果你只是想确认当前仓库能正常跑起来，建议直接使用下面的顺序：

```bash
cmake -S . -B build
cmake --build build -j
cd build
./twoPhaseSolver ../data/inputs2phase2d max_step=10
```

这个流程最适合先验证编译链、AMReX 安装和输入文件解析是否正常。
