@page runtime-guide 当前源码版本使用指南
@ingroup analysis_docs

@tableofcontents

## 1. 文档适用范围

这份文档只说明当前仓库源码的真实行为，而不是历史文档中的设计目标。

写这份文档的原因很直接：仓库里部分旧文档仍然保留了早期项目描述，和现在真正编译运行的求解器并不一致。下面的说明以当前源码、当前 CMake 配置以及一次实际启动验证为准。

当前确认过的事实：

- 当前启用的主可执行文件是 `twoPhaseSolver`。
- 当前构建配置链接的是二维 AMReX 安装，因此实际运行维度是 2D。
- 当前算例的运行参数主要来自 `data/` 目录下的 inputs 文件。
- 当前双相算例的初始条件不是由 inputs 文件完全控制，而是部分硬编码在源码中。

## 2. 当前仓库实际在运行什么

### 2.1 主程序和构建入口

- 主程序入口：`main.cpp`
- CMake 构建入口：`CMakeLists.txt`
- 当前启用的目标：`twoPhaseSolver`

当前 `CMakeLists.txt` 收集的是下面这些源码：

- `main.cpp`
- `2phase/*.cpp`
- `amrex_level_Cong/*.cpp`
- `reconstruction_scheme/*.cpp`
- `deriv/*.cpp`

这说明当前主线不是旧 README 里描述的简单二维不可压 Navier-Stokes 求解器，而是一个基于 AMReX 的双相可压缩求解器。

### 2.2 当前运行维度

当前 CMake 里把 `AMReX_DIR` 指向了 `/opt/amrex_2d/lib/cmake/AMReX`，因此当前构建出的 `twoPhaseSolver` 是二维版本。

这有一个直接后果：

- `data/inputs2phase2d`
- `data/inputs2phase2dMach10`

这两个二维算例可以直接用于当前 build。

而下面这些一维算例文件：

- `data/inputs2phase1d`
- `data/inputs2phase1dforRiemannProb`
- `data/inputs2phase1dforShockInterface`

更适合在一维 AMReX 版本下重新构建后使用，不能默认当作当前二维二进制的主用法。

## 3. 环境要求

当前源码版本至少需要下面这些依赖：

- Linux 或 WSL
- CMake
- C++ 编译器
- MPI
- AMReX

其中最关键的是 AMReX 安装维度必须和你想运行的求解器一致。当前仓库默认写的是二维 AMReX。

## 4. 编译方法

如果你已经有可用的 `build/` 目录，可以直接在现有目录中重新配置和编译。

推荐流程如下：

1. 进入项目根目录。
2. 确认 `CMakeLists.txt` 中的 `AMReX_DIR` 指向正确的安装路径。
3. 创建或进入 `build/` 目录。
4. 运行 CMake。
5. 编译 `twoPhaseSolver`。

示例命令：

    cd /path/to/CFD_inWSL_MF_initial
    mkdir -p build
    cd build
    cmake ..
    cmake --build . -j

如果你只是修改了 inputs 文件，一般不需要重新编译。

如果你改了下面这些源码文件，就需要重新编译：

- `2phase/equation.H`
- `2phase/Prob_Parm.H`
- `2phase/bc_nullfill.cpp`
- `2phase/erroEst.cpp`
- 以及其他求解器实现文件

## 5. 运行方法

### 5.1 推荐运行位置

建议从 `build/` 目录运行，这样输出文件会集中落在 `build/` 目录下，不会把项目根目录弄乱。

### 5.2 当前已验证可用的启动方式

下面这个命令已经实际验证可以启动当前二维双相算例：

    cd build
    ./twoPhaseSolver ../data/inputs2phase2d

如果只是想试跑少量步数，可以在命令后面覆盖参数：

    ./twoPhaseSolver ../data/inputs2phase2d max_step=10

也可以临时关闭 plotfile 和 checkpoint：

    ./twoPhaseSolver ../data/inputs2phase2d max_step=10 amr.plot_files_output=0 amr.checkpoint_files_output=0

### 5.3 输出文件会写到哪里

当前程序有两类输出：

1. AMReX 标准输出

- plotfile，例如 `plt00000`
- checkpoint，例如 `chk00000`

2. 当前项目自定义 ASCII 输出

- `euler_data_level_0.dat`
- `euler_data_level_1.dat`
- `euler_data_level_2.dat`

第二类输出是当前项目在 `output_tmp.H` 里主动写出的，因此即使你把 plotfile 和 checkpoint 关掉，ASCII 文件仍然可能继续生成。

## 6. 当前可直接使用的算例文件

### 6.1 当前二维 build 直接可跑

#### `data/inputs2phase2d`

这是当前最适合入门和修改的二维双相算例，实际可以启动。

特点：

- 二维区域
- 启用了 AMR
- 边界条件比较完整
- 对应当前源码中的二维激波-水柱型初始条件

#### `data/inputs2phase2dMach10`

这是另一个二维双相算例模板，适合在现有二维 build 下继续改。

### 6.2 当前仓库里有但不建议直接拿当前二维二进制跑

#### `data/inputs2phase1d`
#### `data/inputs2phase1dforRiemannProb`
#### `data/inputs2phase1dforShockInterface`

这些文件更像是一维版本的输入模板。

如果你想真正跑一维问题，建议：

1. 把构建改成一维 AMReX。
2. 重新编译求解器。
3. 再使用这些一维 inputs。

## 7. 当前算例是怎么被配置进去的

当前双相求解器的算例配置分成三层：

1. inputs 文件控制运行参数、网格、边界、AMR 选项。
2. `2phase/equation.H` 控制初始条件。
3. `2phase/Prob_Parm.H` 控制两相物性参数。

这三层必须区分清楚，否则很容易出现“明明改了参数，但结果没有变化”的情况。

## 8. 如何修改算例

### 8.1 改时间推进和停止条件

直接修改 inputs 文件中的这些参数：

- `max_step`
- `stop_time`
- `adv.cfl`
- `adv.do_reflux`
- `adv.v`

其中当前源码里真正从 `adv` 段读取的主要是：

- `v`
- `cfl`
- `do_reflux`

所以这些项改完会直接影响实际运行。

### 8.2 改计算域和网格

直接修改 inputs 文件中的这些项：

- `geometry.prob_lo`
- `geometry.prob_hi`
- `geometry.is_periodic`
- `amr.n_cell`
- `amr.max_level`
- `amr.ref_ratio`
- `amr.regrid_int`
- `amr.blocking_factor`
- `amr.max_grid_size`

这是最适合放在 inputs 文件里的内容，也是当前源码真正支持的配置方式。

### 8.3 改初始条件

这是当前仓库最容易误解的地方。

当前双相初始条件不是完全由 inputs 文件定义，而是硬编码在 `2phase/equation.H` 的 `initdata` 函数里。

#### 一维初值位置

在 `2phase/equation.H` 里的一维 `initdata` 分支。

#### 二维初值位置

在 `2phase/equation.H` 里二维 `initdata` 分支。

当前二维分支大致把区域分成三类：

- 左侧激波后区域
- 圆柱内部液相区域
- 外部气相区域

你通常会在这里修改下面这些量：

- 激波位置
- 圆柱中心坐标
- 圆柱半径
- 各区域的密度
- 各区域的速度
- 各区域的压力
- 各区域的体积分数

也就是说，如果你想把“激波打水柱”改成“两个圆斑”或者“分层界面”，主要工作不在 inputs，而在 `2phase/equation.H`。

### 8.4 改边界条件

当前边界条件是通过 inputs 文件中的 `amr.xxx_lo_bc_dir` 和 `amr.xxx_hi_bc_dir` 这类参数传入的。

双相算例当前会读取这些变量的边界：

- `density1`
- `density2`
- `Momentum_X`
- `Momentum_Y`
- `energy`
- `volumeFraction1`

例如在二维问题里，你会看到这样的字段：

- `amr.density1_lo_bc_0`
- `amr.Momentum_X_hi_bc_0`
- `amr.Momentum_Y_lo_bc_1`
- `amr.energy_hi_bc_1`

#### 当前输入文件里边界类型常用数值

- `-1`：反射奇对称，常用于法向速度
- `1`：反射偶对称，常用于密度、压力、切向速度
- `2`：一阶外推
- `3`：外部 Dirichlet，需要用户自己填 ghost cell

#### 当前源码下边界修改的实际建议

如果你只是想改成常见的外推或反射边界，可以直接改 inputs 文件中的边界类型数值。

如果你想做真正的定值入口、定值出口或者用户自定义边界，就不能只改数字，因为当前 `2phase/bc_nullfill.cpp` 基本是空实现。也就是说：

- 反射、外推、周期这类 AMReX 内建行为可以直接用。
- `ext_dir`、`user_1` 这类需要用户填 ghost cell 的边界，当前仓库还没有真正补全实现。

如果你要做这类边界，需要同时修改 `2phase/bc_nullfill.cpp`。

### 8.5 改两相物性参数

当前两相物性参数定义在 `2phase/Prob_Parm.H` 中，包括：

- 两相 `gamma`
- 两相 `pInf`
- `alphaEps`

当前源码并没有把这些量从 inputs 文件读进来，所以：

- 改 inputs 文件里的同名概念，当前通常不会生效。
- 真正要改，应该改 `2phase/Prob_Parm.H` 后重新编译。

如果你后面希望把物性参数也改成“只改 inputs 即可”，那需要继续改源码，把它们接入 `ParmParse`。

### 8.6 改 AMR 加密行为

当前 AMR 的几何级参数仍然建议在 inputs 里改：

- `amr.max_level`
- `amr.ref_ratio`
- `amr.regrid_int`
- `amr.blocking_factor`
- `amr.max_grid_size`

但具体“哪里加密”的判据主要在 `2phase/erroEst.cpp` 里。

当前实现更偏向于基于总密度变化来做 TENO 标记，不是一个纯粹由 `tagging.phierr` 或 `tagging.phigrad` 完全控制的简单阈值器。

所以如果你只是改：

- `tagging.phierr`
- `tagging.max_phierr_lev`
- `tagging.phigrad`
- `tagging.max_phigrad_lev`

不一定能得到你预期中的细化变化。

如果你想让 AMR 明确跟着相界面、压力梯度或者激波走，通常还要直接修改 `2phase/erroEst.cpp`。

## 9. 当前源码里哪些参数容易让人误判

### 9.1 `adv.num_state`

inputs 文件里有这个参数，但当前运行时会报它是未使用变量。

原因是当前状态量个数由源码里的 `NUM_STATE = 4 + BL_SPACEDIM` 决定，不是由 inputs 驱动。

### 9.2 `adv.do_tracers`

inputs 文件里也有这个参数，但当前运行时同样会显示为未使用变量。

### 9.3 物性参数不是 inputs 驱动

旧文档里可能会让人以为双相 `gamma` 和 `pInf` 可以直接在 inputs 里配，但当前源码并没有这样实现。

### 9.4 一维 inputs 不是当前二维二进制的主用法

仓库里同时放着 1D 和 2D inputs，但当前默认 CMake 指向二维 AMReX，因此当前 build 更适合直接跑二维算例。

## 10. 推荐的实际改算例流程

如果你要基于现有二维算例做自己的问题，建议按这个顺序来：

1. 复制 `data/inputs2phase2d` 为一个新的 inputs 文件。
2. 先只修改计算域、网格、时间步、输出频率和边界类型。
3. 运行一次，确认网格和边界行为正常。
4. 再修改 `2phase/equation.H` 中的二维 `initdata`，把初始条件改成你的目标问题。
5. 如果问题涉及不同材料参数，再修改 `2phase/Prob_Parm.H`。
6. 重新编译。
7. 再次运行并观察 `euler_data_level_*.dat`、plotfile 或 checkpoint 输出。

这是当前源码版本下最省事、最不容易误判的一套工作流。

## 11. 如果你后面想继续改进这个仓库

从“便于用”的角度看，当前最值得优先补的两项是：

1. 把初始条件从 `2phase/equation.H` 挪一部分到 inputs 文件中。
2. 把 `gamma`、`pInf` 等物性参数接入 `ParmParse`。

这样后面换算例时，就不用每次都改源码再重新编译。