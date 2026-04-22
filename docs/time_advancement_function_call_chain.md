@page time-advance-chain 2phase 模型流场变量时间推进函数调用关系详解
@ingroup analysis_docs

@tableofcontents

## 1. 概述

本文档详细说明基于 AMReX 框架的两相流 CFD 求解器中，流场变量从当前时刻 $t^n$ 推进到下一时刻 $t^{n+1}$ 所涉及的完整函数调用关系。

**核心物理模型**：Allaire 五方程两相流模型  
**数值方法**：有限体积法 (FVM) + Runge-Kutta 时间积分  
**空间离散**：高阶重构格式 (THINC/WENO) + HLLC Riemann 求解器  
**网格架构**：自适应网格加密 (AMR)

---

## 2. 总体架构图

```
┌─────────────────────────────────────────────────────────┐
│  AmrLevelCong::advance(time, dt, iteration, ncycle)    │
│                     【时间推进主循环】                    │
└────────────────┬────────────────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────────────────────┐
│  1. swapTimeLevels(dt)  // 交换时间层数据                │
│  2. get_new_data()      // 获取 S_new (n+1 时刻)          │
│  3. get_old_data()      // 获取 S_old (n 时刻)            │
└────────────────┬────────────────────────────────────────┘
                 │
                 ▼
┌─────────────────────────────────────────────────────────┐
│  RK(RKOrder=3, state, time, dt,                         │
│     compute_dSdt<RKOrder>, ...)                          │
│           【三阶 Runge-Kutta 多阶段循环】                   │
└────────────────┬────────────────────────────────────────┘
                 │
        ┌────────┴────────┬────────────────┐
        │  阶段 1          │  阶段 2          │  阶段 3
        ▼                 ▼                 ▼
┌─────────────────────────────────────────────────────────┐
│  compute_dSdt(stage, Sborder, dSdt, dt_sub, ...)        │
│              【计算时间导数 ∂U/∂t】                        │
└────────────────┬────────────────────────────────────────┘
                 │
                 ├──────────────────┬──────────────────┐
                 ▼                  ▼                  ▼
        ┌─────────────────┐ ┌──────────────┐ ┌─────────────────┐
        │守恒变量→原始变量│ │  idim=0 (x)  │ │  idim=1 (y)     │
        │cons_to_prim     │ │flux_solve    │ │flux_solve       │
        └─────────────────┘ └──────────────┘ └─────────────────┘
                                              (2D 情况)
                 
                 │
                 ▼
┌─────────────────────────────────────────────────────────┐
│  flux_solve_in_dim<idim>                                │
│          【空间离散：计算 -∇·F(U)】                        │
│  ┌──────────────────────────────────────────┐           │
│  │ 阶段 1: 节点通量计算                       │           │
│  │   prim_to_flux → F_i                     │           │
│  │                                          │           │
│  │ 阶段 2: Riemann 求解器                      │           │
│  │   reconstruction → L,R 状态               │           │
│  │   Riemann_solver → F_{i+1/2}             │           │
│  │                                          │           │
│  │ 阶段 3: 通量差分更新 RHS                   │           │
│  │   ΔF/Δx → ∂U/∂t                          │           │
│  └──────────────────────────────────────────┘           │
└────────────────┬────────────────────────────────────────┘
                 │
                 ▼
          返回 dS/dt
                 │
                 ▼
┌─────────────────────────────────────────────────────────┐
│  RK 组合：S_new = Σ(ω_k * S_stage_k)                    │
│          LinComb3(...)                                   │
└────────────────┬────────────────────────────────────────┘
                 │
                 ▼
        S_new (n+1 时刻守恒变量)
```

---

## 3. 详细函数调用链

### 3.1 时间推进入口

#### 📍 `AmrLevelCong::advance()`

**文件位置**: `amrex_level_Cong/AmrLevelCong_advance.cpp:30-130`

```cpp
Real AmrLevelCong::advance(Real time, Real dt, int iteration, int ncycle)
```

**功能描述**:  
时间推进主函数，负责协调 Runge-Kutta 多阶段积分过程，管理新旧数据交换和通量回流。

**详细步骤**:

##### 步骤 1: 打印统计信息（仅 level 0）
```cpp
if (level == 0) {
    MultiFab& S_mm = get_new_data(State_Type);
    for (int ivar = 0; ivar < NVar; ++ivar) {
        Real maxval = S_mm.max(ivar);
        Real minval = S_mm.min(ivar);
        amrex::Print() << "Variable " << ivar
                       << ": max = " << maxval
                       << ", min = " << minval << std::endl;
    }
}
```

##### 步骤 2: 交换时间层
```cpp
for (int k = 0; k < NUM_STATE_TYPE; k++) {
    state[k].allocOldData();      // 分配旧数据内存
    state[k].swapTimeLevels(dt);  // 交换新旧时间层
}
```

##### 步骤 3: 获取数据指针
```cpp
MultiFab& S_new = get_new_data(State_Type);  // n+1 时刻状态
MultiFab& S_old = get_old_data(State_Type);  // n 时刻状态

const Real prev_time = state[State_Type].prevTime();
const Real cur_time = state[State_Type].curTime();
const Real ctr_time = 0.5 * (prev_time + cur_time);
```

##### 步骤 4: 准备通量寄存器
```cpp
FluxRegister* fine = 0;
FluxRegister* current = 0;

int finest_level = parent->finestLevel();

if (do_reflux && level < finest_level) {
    fine = &getFluxReg(level + 1);
    fine->setVal(0.0);
}

if (do_reflux && level > 0) {
    current = &getFluxReg(level);
}
```

##### 步骤 5: 调用 Runge-Kutta 方法
```cpp
constexpr int RKOrder = 3;

RK(RKOrder, State_Type, time, dt, iteration, ncycle,
    // Lambda 1: 计算时间导数
    [&](int stage, MultiFab& dSdt, MultiFab const& S, Real t, Real dtsub) {
        compute_dSdt<RKOrder>(stage, S, dSdt, dtsub, current, fine);
    },
    // Lambda 2: 每阶段后处理（空操作）
    [&](int /*stage*/, MultiFab& S) {}
);
```

**返回值**: 实际使用的时间步长 `dt`

---

### 3.2 计算时间导数

#### 📍 `compute_dSdt<RKOrder>()`

**文件位置**: `amrex_level_Cong/AmrLevelCong_advance.cpp:178-392`

```cpp
template <int RKOrder>
void AmrLevelCong::compute_dSdt(int stage, 
                                 const MultiFab& Sborder, 
                                 MultiFab& dSdt, 
                                 Real dt,
                                 FluxRegister* current, 
                                 FluxRegister* fine)
```

**功能描述**:  
计算状态变量的时间导数 $\frac{\partial U}{\partial t}$，包括空间离散和通量计算。

#### 3.2.1 数据准备

##### 守恒变量→原始变量转换
```cpp
// 获取旧数据（用于保正限制器）
const MultiFab& old_data = get_old_data(State_Type);

// 创建原始变量数组（带鬼细胞）
MultiFab prim(grids, dmap, NUM_STATE, NUM_GROW, MFInfo(), Factory());

// 对每个网格块执行转换
for (MFIter mfi(Sborder, false); mfi.isValid(); ++mfi) {
    Box ghost_box = amrex::grow(mfi.validbox(), NUM_GROW);
    Array4<Real> const& prim_array = prim[mfi].array();
    Array4<Real const> const& statein_array = Sborder[mfi].const_array();
    ProbParm const& pp = *d_prob_parm;
    
    // 🔑 核心转换函数
    cons_to_prim_for_array4(ghost_box, prim_array, statein_array, pp);
}
```

##### 初始化时间导数和通量数组
```cpp
// 初始化 dS/dt = 0
dSdt.setVal(0.0);

// 创建通量数组
MultiFab fluxes[BL_SPACEDIM];
GpuArray<Real, BL_SPACEDIM> dx = geom.CellSizeArray();

// 如需回流，定义通量寄存器
if (do_reflux) {
    for (int j = 0; j < BL_SPACEDIM; j++) {
        BoxArray ba = prim.boxArray();
        ba.surroundingNodes(j);
        fluxes[j].define(ba, dmap, NUM_STATE, 0);
    }
}
```

---

#### 3.2.2 核心转换函数详解

##### 📍 `cons_to_prim_for_array4()`

**文件位置**: `2phase/equation.H:176-197`

```cpp
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void cons_to_prim_for_array4(const amrex::Box& ghost_box,
    amrex::Array4<amrex::Real> const& prim_array,
    amrex::Array4<amrex::Real const> const& statein_array,
    ProbParm const& pp = ProbParm())
{
    ParallelFor(ghost_box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        int ncomp = NVar;
        VarArray input;  // 输入守恒变量
        
        // 读取守恒变量
        for (int ivar = 0; ivar < ncomp; ivar++) {
            input[ivar] = statein_array(i, j, k, ivar);
        }
        
        // 🔑 转换为原始变量
        VarArray primi = cons_to_prim(input, pp);
        
        // 写入原始变量数组
        for (int ivar = 0; ivar < ncomp; ivar++) {
            prim_array(i, j, k, ivar) = primi[ivar];
        }
    });
}
```

**参数说明**:
- `[in] ghost_box`: 包含鬼细胞的计算区域
- `[out] prim_array`: 原始变量输出数组
- `[in] statein_array`: 守恒变量输入数组
- `[in] pp`: 问题参数（包含状态方程参数）

---

##### 📍 `cons_to_prim()`

**文件位置**: `2phase/equation.H:113-137`

```cpp
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
VarArray cons_to_prim(VarArray cons, ProbParm const& pp)
{
    using namespace amrex;

    // 提取守恒变量
    Real rho1 = cons[AlphaRho1];   // 相 1 密度
    Real rho2 = cons[AlphaRho2];   // 相 2 密度
    Real AMREX_D_DECL(rhou = cons[XMom],   // x 动量
                      rhov = cons[YMom],   // y 动量
                      rhow = cons[ZMom]);  // z 动量
    Real rhoE = cons[Energy];      // 总能量
    Real alpha1 = cons[Alpha1];    // 相 1 体积分数

    // 计算总密度
    Real rho = rho1 + rho2;

    // 计算速度分量
    Real AMREX_D_DECL(u = rhou / rho,   // x 速度
                      v = rhov / rho,   // y 速度
                      w = rhow / rho);  // z 速度

    // 🔑 通过状态方程计算压力
    Real pressure = cons_to_eos<Pressure>(cons, pp);

    // 返回原始变量 {ρ₁, ρ₂, u, v, w, p, α₁}
    return { rho1, rho2, AMREX_D_DECL(u, v, w), pressure, alpha1 };
}
```

**变量映射关系**:
| 守恒变量 | 符号 | 原始变量 | 符号 |
|---------|------|---------|------|
| AlphaRho1 | $\alpha_1\rho_1$ | AlphaRho1 | $\alpha_1\rho_1$ |
| AlphaRho2 | $\alpha_2\rho_2$ | AlphaRho2 | $\alpha_2\rho_2$ |
| XMom | $\rho u$ | XMom | $u$ |
| YMom | $\rho v$ | YMom | $v$ |
| Energy | $\rho E$ | Energy | $p$ |
| Alpha1 | $\alpha_1$ | Alpha1 | $\alpha_1$ |

---

##### 📍 `cons_to_eos<EOSType::Pressure>()`

**文件位置**: `2phase/equation.H:23-64`

```cpp
template <EOSType eType>
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real cons_to_eos(VarArray cons, ProbParm const& pp)
{
    using namespace amrex;
    
    // 提取守恒变量
    Real rho1 = cons[AlphaRho1], rho2 = cons[AlphaRho2];
    Real AMREX_D_DECL(rhou = cons[XMom], rhov = cons[YMom], rhow = cons[ZMom]);
    Real rhoE = cons[Energy], alpha1 = cons[Alpha1];

    // 计算动能
    Real halfRhoU_sqr = (AMREX_D_TERM(rhou*rhou, +rhov*rhov, +rhow*rhow)) / 2;

    // 计算混合区 gamma 和 pInf
    amrex::GpuArray<amrex::Real, 2> gamma_pInf = pp.cal_gamma(alpha1);
    Real gamma_bar = gamma_pInf[0];    // 混合绝热指数
    Real pInf_bar = gamma_pInf[1];     // 混合参考压力

    // 总密度
    Real rho = rho1 + rho2;

    // Lambda 函数：计算状态量
    auto computeState = [&](Real gamma, Real pInf, Real rho) -> Real {
        // 内能 = 总能 - 动能
        Real internalEnergy = rhoE - halfRhoU_sqr / rho;
        
        // 压力 = (γ-1)*内能 - γ*pInf
        Real pressure = (gamma - 1) * internalEnergy - gamma * pInf;
        
        // 声速² = γ*(γ-1)*(内能-pInf)/ρ
        Real soundSpeedSqrDensity = gamma * (gamma - 1) * (internalEnergy - pInf);
        
        // 声速
        Real soundSpeed = sqrt(soundSpeedSqrDensity / rho);

        // 根据模板参数返回相应状态量
        if constexpr (eType == Pressure) {
            return pressure;
        } else if constexpr (eType == SoundSpeed) {
            return soundSpeed;
        } else if constexpr (eType == SoundSpeedSqrDensity) {
            return soundSpeedSqrDensity;
        }
    };

    return computeState(gamma_bar, pInf_bar, rho);
}
```

**状态方程类型** (`EOSType` 枚举):
```cpp
enum EOSType {
    Pressure = 0,           // 压力 p
    SoundSpeed,             // 声速 c
    SoundSpeedSqrDensity,   // 声速平方×密度 ρc²
    TotalEnergy             // 总能量 ρE
};
```

---

### 3.3 分维度通量计算

#### 📍 `flux_solve_in_dim<idim>()`

**文件位置**: `amrex_level_Cong/AmrLevelCong_advance.cpp:395-590`

```cpp
template <int idim>
void AmrLevelCong::flux_solve_in_dim(const FluxSolveParams& params)
```

**功能描述**:  
在指定空间维度上计算通量散度，实现空间离散化。采用**三阶段流程**。

**参数结构体** (`FluxSolveParams`):
```cpp
struct FluxSolveParams {
    const amrex::Box& bx;                                    // 计算区域
    amrex::GpuArray<amrex::FArrayBox*, AMREX_SPACEDIM>& flux;  // 通量数组
    amrex::GpuArray<amrex::FArrayBox*, AMREX_SPACEDIM>& ncflux;// 非守恒通量
    amrex::FArrayBox& cellFlux;                              // 单元格通量
    amrex::FArrayBox& cellNcFlux;                            // 单元格非守恒通量
    const amrex::FArrayBox& primin;                          // 原始变量
    const amrex::FArrayBox& consin;                          // 守恒变量
    amrex::FArrayBox& rhs;                                   // 右端项 (dS/dt)
    const amrex::GpuArray<amrex::Box, BL_SPACEDIM>& nbx;      // 邻居边界框
    amrex::Real dt;                                          // 时间步长
    amrex::Real RKTimeCoef;                                  // RK 时间系数
};
```

---

#### 阶段 1: 节点通量计算

```cpp
// 扩展通量计算区域（包含虚拟网格）
Box fluxBx = amrex::grow(bx, num_flux_ghost_cells);

// 初始化节点通量为 0
cellFlux.setVal<RunOn::Device>(0.0);

// 获取数组视图
Array4<Real const> const& primin_array = primin.const_array();
Array4<Real> const& nodeFlux_array = cellFlux.array();

// ----------------------------
// 阶段 1：计算节点通量 F_i
// ----------------------------
ParallelFor(fluxBx, [=] AMREX_GPU_DEVICE(int xIdx, int yIdx, int zIdx) {
    // 1. 提取原始状态变量
    VarArray primitiveState;
    for (int stateIdx = 0; stateIdx < NUM_STATE; stateIdx++) {
        primitiveState[stateIdx] = primin_array(xIdx, yIdx, zIdx, stateIdx);
    }

    // 2. 🔑 计算 Euler 通量
    VarArray fluxState = prim_to_flux(primitiveState, pp, norm);
    
    // 3. 存储到节点通量数组
    for (int stateIdx = 0; stateIdx < NUM_STATE; stateIdx++) {
        nodeFlux_array(xIdx, yIdx, zIdx, stateIdx) = fluxState[stateIdx];
    }
});
```

##### 📍 `prim_to_flux()`

**文件位置**: `2phase/equation.H:200-237`

```cpp
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
VarArray prim_to_flux(const VarArray& prim, 
                      ProbParm const& pp, 
                      amrex::GpuArray<amrex::Real, 3> norm)
{
    using namespace amrex;

    // 提取原始变量
    Real rho1 = prim[AlphaRho1], rho2 = prim[AlphaRho2];
    Real AMREX_D_DECL(u = prim[XMom], v = prim[YMom], w = prim[ZMom]);
    Real pressure = prim[Energy], alpha1 = prim[Alpha1];

    // 计算总密度
    Real rho = rho1 + rho2;

    // 计算法向速度 Vn = u·n
    Real Vn = AMREX_D_TERM(u*norm[0], +v*norm[1], +w*norm[2]);

    // 🔑 计算总能量 ρE
    Real rhoE = prims_to_eos<TotalEnergy>(prim, pp);

    // 构造通量向量 F = [ρ₁Vn, ρ₂Vn, ρuVn+p·n, ρvVn+p·n, ρwVn+p·n, (ρE+p)Vn, 0]ᵀ
    VarArray flux = {
        rho1 * Vn,                          // 相 1 密度通量
        rho2 * Vn,                          // 相 2 密度通量
        AMREX_D_DECL(
            rho*u*Vn + pressure*norm[0],    // x 动量通量
            rho*v*Vn + pressure*norm[1],    // y 动量通量
            rho*w*Vn + pressure*norm[2]),   // z 动量通量
        (rhoE + pressure) * Vn,             // 能量通量
        0.                                  // 相分数通量 (非守恒)
    };

    return flux;
}
```

**通量向量物理意义** (以 x 方向为例，norm=[1,0,0]):
$$
F_x = \begin{bmatrix}
\alpha_1\rho_1 u \\
\alpha_2\rho_2 u \\
\rho u^2 + p \\
\rho uv \\
\rho uw \\
(\rho E + p)u \\
0
\end{bmatrix}
$$

---

#### 阶段 2: Riemann 求解器（半节点通量）

```cpp
// 获取法向量（如 x 方向：norm=[1,0,0]）
ProbParm const& pp = *d_prob_parm;
amrex::GpuArray<amrex::Real, 3> norm = get_direction<idim>();
GpuArray<Real, BL_SPACEDIM> dx = geom.CellSizeArray();

// 获取数组视图
Array4<Real> const& fluxRef = flux[idim]->array();
Array4<Real> const& ncfluxRef = ncflux[idim]->array();

// ----------------------------
// 阶段 2：Riemann 求解器计算 F_{i+1/2}
// ----------------------------
ParallelFor(nbx[idim], [=] AMREX_GPU_DEVICE(int i, int j, int k) {
    
    // 🔑 1. 空间重构获取左右状态 U_L, U_R
    amrex::GpuArray<VarArray, 2> stateLR = reconstruction<idim>(
        i, j, k, primin_array, NUM_STATE, pp, norm, dx);

    // 🔑 2. Riemann 求解器计算界面通量
    VarArray iFlux = Riemann_solver_1<idim>(
        i, j, k, primin_array, stateLR[0], stateLR[1], norm, pp, dx);

    // 3. 提取左右单元守恒变量和通量
    VarArray consL, consR, fluxCellL, fluxCellR;
    Dim3 offset = amrex::IntVect::TheDimensionVector(idim).dim3();
    
    for (int ivar = 0; ivar < NUM_STATE; ivar++) {
        consL[ivar] = consin_array(i-offset.x, j-offset.y, k-offset.z, ivar);
        consR[ivar] = consin_array(i, j, k, ivar);
        fluxCellL[ivar] = nodeFlux_array(i-offset.x, j-offset.y, k-offset.z, ivar);
        fluxCellR[ivar] = nodeFlux_array(i, j, k, ivar);
    }
    
    // 4. 🔑 四阶精度通量修正
    for (int ivar = 0; ivar < NUM_STATE; ivar++) {
        iFlux[ivar] = 4.0/3.0 * iFlux[ivar]
            - 1.0/6.0 * (nodeFlux_array(i-offset) + nodeFlux_array(i));
    }

    // 5. 计算非守恒项通量（α·du/dx 中的 u_{1/2}）
    AdvectionTerm advTerm;
    auto flux_var = advTerm.computeFlux(stateLR[0], stateLR[1], pp, norm, dx);
    auto indices = advTerm.getTargetTermIndices();

    // 6. 计算左右单元的系数 α
    VarArray primL, primR;
    for (int ivar = 0; ivar < NUM_STATE; ivar++) {
        primL[ivar] = primin_array(i-offset.x, j-offset.y, k-offset.z, ivar);
        primR[ivar] = primin_array(i, j, k, ivar);
    }
    auto coefL = advTerm.computeCoefficient(primL, pp, norm);
    auto coefR = advTerm.computeCoefficient(primR, pp, norm);

    // 7. 计算增量通量
    VarArray incFluxL, incFluxR;
    incFluxL.fill(0.0);
    incFluxR.fill(0.0);
    
    for (int idx = 0; idx < AdvectionTerm::num_terms; idx++) {
        incFluxL[indices[idx]] += coefL[idx]
            * (4.0/3.0 * flux_var[idx]
                - 1.0/6.0 * (nodencFlux_array(i-offset, idx) + nodencFlux_array(i, idx)));
        incFluxR[indices[idx]] += coefR[idx]
            * (4.0/3.0 * flux_var[idx]
                - 1.0/6.0 * (nodencFlux_array(i-offset, idx) + nodencFlux_array(i, idx)));
    }

    // 8. 计算中间状态（用于保正检测）
    Real multL = 2. * dt * RKTimeCoef / dx[idim] / computeAlpha(dx, primL, pp)[idim];
    Real multR = -2. * dt * RKTimeCoef / dx[idim] / computeAlpha(dx, primR, pp)[idim];

    VarArray consStarL, consStarR;
    PositiveQuantityFunctions pqf;
    for (int ivar = 0; ivar < NUM_STATE; ivar++) {
        consStarL[ivar] = consL[ivar] - multL * (iFlux[ivar] + incFluxL[ivar] - fluxCellL[ivar]);
        consStarR[ivar] = consR[ivar] - multR * (iFlux[ivar] + incFluxR[ivar] - fluxCellR[ivar]);
    }

    // 9. 🔑 保正限制器（如出现负密度/压力）
    bool leftFlag = containsInvalidValues(consStarL, pp, pqf);
    bool rightFlag = containsInvalidValues(consStarR, pp, pqf);
    
    if (leftFlag || rightFlag) {
        // 切换到一阶格式（HLL）
        VarArray foFlux = Riemann_solver_2(i, j, k, primin_array, primL, primR, norm, pp);
        
        // 计算限制器参数 θ
        Real thetaL = compute_theta(consStarL, consFoL, pp, pqf);
        Real thetaR = compute_theta(consStarR, consFoR, pp, pqf);
        Real theta = min(thetaL, thetaR);
        
        // 凸组合修正：F_final = θ*F_HLL + (1-θ)*F_HLLC
        apply_convex_combination(foFlux, iFlux, theta, 0, NVar);
        apply_convex_combination(foncFluxL, incFluxL, theta, 0, NVar);
        apply_convex_combination(foncFluxR, incFluxR, theta, 0, NVar);
    }

    // 10. 存储最终通量
    for (int ivar = 0; ivar < NVar; ivar++) {
        fluxRef(i, j, k, ivar) = iFlux[ivar];        // 守恒通量
        ncfluxRef(i, j, k, ivar) = incFluxL[ivar];   // 左侧非守恒通量
        ncfluxRef(i, j, k, ivar+NVar) = incFluxR[ivar]; // 右侧非守恒通量
    }
});
```

---

##### 📍 `reconstruction<idim>()`

**文件位置**: `reconstruction_scheme/fifth_order_thinc_recon.H:98-...`

```cpp
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::GpuArray<VarArray, 2>
primitive_reconstruction(int i, int j, int k,
    amrex::Array4<amrex::Real const> const& primin_array,
    int NUM_STATE, ProbParm const& pp,
    amrex::GpuArray<amrex::Real, 3> norm,
    amrex::GpuArray<amrex::Real, BL_SPACEDIM> dx)
```

**功能**: 使用五阶 THINC/WENO 格式在网格界面重构左右状态。

**重构流程** (以标量 φ 为例):

1. **构建 stencil**: $\{\phi_{i-2}, \phi_{i-1}, \phi_i, \phi_{i+1}, \phi_{i+2}\}$

2. **THINC 插值** (双曲正切重构):
```cpp
AMREX_GPU_HOST_DEVICE
amrex::Real phi_THINC(amrex::Real phi_im1, amrex::Real phi_i, 
                      amrex::Real phi_ip1, amrex::Real x, 
                      amrex::Real beta = 1.8)
{
    auto sgn = [](amrex::Real val) -> amrex::Real {
        return val > 0 ? 1.0 : (val < 0 ? -1.0 : 0.0);
    };

    amrex::Real min_phi = std::min(phi_im1, phi_ip1);
    amrex::Real max_phi = std::max(phi_im1, phi_ip1);
    amrex::Real phi_bar = (max_phi + min_phi) / 2;
    amrex::Real phi_delta = (max_phi - min_phi) / 2;
    amrex::Real theta = sgn(phi_ip1 - phi_im1);

    // 单调性检测
    if ((phi_im1 - phi_i) * (phi_i - phi_ip1) > 0) {
        // 计算界面位置
        amrex::Real di = atanh(1.0 - 1e-8 - (2.0 - 2e-8) * phi_i) / beta;
        
        // THINC 重构公式
        return 0.5 + (0.5 - 1e-8) * theta * tanh(beta * (x - theta * di));
    } else {
        return phi_i;  // 非单调区域用中心值
    }
}
```

3. **WENO 非线性加权** (未在代码片段中显示完整实现):
   - 计算多个子 stencil 的平滑度指示器
   - 计算非线性权重
   - 加权得到最终重构值

4. **返回左右状态**:
```cpp
return { stateL, stateR };  // stateL: 左状态，stateR: 右状态
```

---

##### 📍 `Riemann_solver_1<idim>()`

**文件位置**: `2phase/equation.H:295-510`

```cpp
template <int idim>
AMREX_FORCE_INLINE
AMREX_GPU_HOST_DEVICE
VarArray
Riemann_solver_1(int i, int j, int k, 
                 const amrex::Array4<const amrex::Real>& prim_array,
                 VarArray const& stateL, VarArray const& stateR,
                 amrex::GpuArray<amrex::Real, 3> norm,
                 ProbParm const& pp,
                 amrex::GpuArray<amrex::Real, BL_SPACEDIM> dx)
```

**功能**: HLLC 近似 Riemann 求解器，计算界面通量。

**详细算法**:

##### 步骤 1: 提取左右状态
```cpp
enum VelComponents { Ar1=0, Ar2, AMREX_D_DECL(U,V,W), Pressure, Alpha1 };

// 左状态
Real ar1l = stateL[Ar1], ar2l = stateL[Ar2];
Real pl = stateL[Pressure], alpha1l = stateL[Alpha1];
Real AMREX_D_DECL(ul = stateL[U], vl = stateL[V], wl = stateL[W]);

// 右状态
Real ar1r = stateR[Ar1], ar2r = stateR[Ar2];
Real pr = stateR[Pressure], alpha1r = stateR[Alpha1];
Real AMREX_D_DECL(uR = stateR[U], vR = stateR[V], wR = stateR[W]);

// 总密度
Real rl = ar1l + ar2l;
Real rr = ar1r + ar2r;
```

##### 步骤 2: 计算声速
```cpp
Real cl = prims_to_eos<SoundSpeed>(stateL, pp);
Real cr = prims_to_eos<SoundSpeed>(stateR, pp);
```

##### 步骤 3: Roe 平均计算
```cpp
// 计算法向速度
GpuArray<Real, 2> Vn = { 
    AMREX_D_TERM(ul*norm[0] + v*norm[1] + w*norm[2]),
    AMREX_D_TERM(uR*norm[0] + vR*norm[1] + wR*norm[2])
};

// Roe 平均密度平方根
Real sqrt_rl = std::sqrt(rl);
Real sqrt_rr = std::sqrt(rr);

// Roe 平均速度
Real VnBar = (Vn[L]*sqrt_rl + Vn[R]*sqrt_rr) / (sqrt_rl + sqrt_rr);

// Roe 平均声速
Real cBar = sqrt((cl*cl*sqrt_rl + cr*cr*sqrt_rr) / (sqrt_rl + sqrt_rr)
    + 0.5*(sqrt_rl*sqrt_rr)/((sqrt_rl+sqrt_rr)^2) * (Vn[R]-Vn[L])^2);
```

##### 步骤 4: 估计波速
```cpp
// Davis 估计
Real SL = std::min(Vn[L] - cl, VnBar - cBar);  // 左行波速
Real SR = std::max(Vn[R] + cr, VnBar + cBar);  // 右行波速
```

##### 步骤 5: 计算中间波速 S*
```cpp
Real Sstar = ((pr - pl) + (rl*Vn[L]*(SL-Vn[L]) - rr*Vn[R]*(SR-Vn[R]))) 
             / (rl*(SL-Vn[L]) - rr*(SR-Vn[R]));
```

##### 步骤 6: HLLC 分段通量
```cpp
VarArray res;

if (SL >= 0) {
    // 情况 1: 所有波向左，取左通量
    res = prim_to_flux(stateL, pp, norm);
    
} else if (Sstar > 0) {
    // 情况 2: 接触间断向左，计算左星区通量
    VarArray UL = prim_to_cons(stateL, pp);
    VarArray FL = prim_to_flux(stateL, pp, norm);
    
    Real Usf = (SL - Vn[L]) / (SL - Sstar);
    
    VarArray Ustar = { 
        Usf * ar1l,
        Usf * ar2l,
        AMREX_D_DECL(
            Usf * (ul + norm[0]*(Sstar-Vn[L])) * rl,
            Usf * (vl + norm[1]*(Sstar-Vn[L])) * rl,
            Usf * (wl + norm[2]*(Sstar-Vn[L])) * rl),
        Usf * (UL[Pressure] + (Sstar-Vn[L])*(rl*Sstar + pl/(SL-Vn[L]))),
        alpha1l 
    };
    
    for (int ivar = 0; ivar < NVar-1; ivar++) {
        res[ivar] = FL[ivar] + SL * (Ustar[ivar] - UL[ivar]);
    }
    
} else if (SR > 0) {
    // 情况 3: 接触间断向右，计算右星区通量
    VarArray UR = prim_to_cons(stateR, pp);
    VarArray FR = prim_to_flux(stateR, pp, norm);
    
    Real Usf = (SR - Vn[R]) / (SR - Sstar);
    
    VarArray Ustar = { 
        Usf * ar1r,
        Usf * ar2r,
        AMREX_D_DECL(
            Usf * (uR + norm[0]*(Sstar-Vn[R])) * rr,
            Usf * (vR + norm[1]*(Sstar-Vn[R])) * rr,
            Usf * (wR + norm[2]*(Sstar-Vn[R])) * rr),
        Usf * (UR[Pressure] + (Sstar-Vn[R])*(rr*Sstar + pr/(SR-Vn[R]))),
        alpha1r 
    };
    
    for (int ivar = 0; ivar < NVar-1; ivar++) {
        res[ivar] = FR[ivar] + SR * (Ustar[ivar] - UR[ivar]);
    }
    
} else {
    // 情况 4: 所有波向右，取右通量
    res = prim_to_flux(stateR, pp, norm);
}

res[Alpha1] = 0;  // 相分数方程非守恒
return res;
```

**HLLC 波结构**:
```
         SL        S*        SR
         ↓         ↓         ↓
    ----●---------●---------●----
   U_L   |   U*_L  |  U*_R   |  U_R
         |         |         |
    左状态  左星区    右星区    右状态
```

---

#### 阶段 3: 通量差分更新 RHS

```cpp
// ----------------------------
// 阶段 3：更新右手边（RHS）∂U/∂t
// ----------------------------
ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
    for (int ivar = 0; ivar < NUM_STATE; ivar++) {
        // 1. 守恒项通量差分：-∂F/∂x
        Real flux_diff = (fluxRef(i+offset) - fluxRef(i)) / dx[idim];
        stateout_array(i, j, k, ivar) -= flux_diff;

        // 2. 非守恒项差分：-∂(αu)/∂x
        // 注意：用右面的左值减去左面的右值
        Real ncflux_diff = (ncfluxRef(i+offset, ivar) - ncfluxRef(i, ivar+NVar)) / dx[idim];
        stateout_array(i, j, k, ivar) -= ncflux_diff;
    }
});
```

**数学形式**:
$$
\frac{\partial U_i}{\partial t} = -\frac{F_{i+1/2} - F_{i-1/2}}{\Delta x} - \left(\frac{\alpha u_{i+1/2} - \alpha u_{i-1/2}}{\Delta x}\right)
$$

---

## 4. 核心物理过程

### 4.1 对流项离散

采用**有限体积法**离散对流项：
$$
\frac{\partial U}{\partial t} + \nabla \cdot F(U) = 0
$$

空间离散格式：
$$
\left.\nabla \cdot F\right|_i \approx \frac{F_{i+1/2} - F_{i-1/2}}{\Delta x}
$$

### 4.2 界面捕捉

**THINC/WENO 混合格式**:
- **THINC**: 在界面区使用双曲正切函数重构，保证锐利界面
- **WENO**: 在光滑区使用五阶 WENO 格式，减少数值耗散

### 4.3 正性保持

**限制器策略**:
1. 预测中间状态 $U^*$
2. 检测是否出现负密度/压力
3. 如无效，切换到一阶 HLL 格式
4. 使用凸组合限制器：
   $$
   F_{final} = \theta F_{HLL} + (1-\theta) F_{HLLC}
   $$
   其中 $\theta \in [0,1]$ 由保正条件确定

### 4.4 时间积分

**三阶 Runge-Kutta 方法** (Shu-Osher 格式):

```
阶段 1: U^(1) = U^n + Δt * L(U^n)
阶段 2: U^(2) = 0.75*U^n + 0.25*U^(1) + 0.25*Δt * L(U^(1))
阶段 3: U^{n+1} = (2/3)*U^n + (1/3)*U^(2) + (2/3)*Δt * L(U^(2))
```

等价于：
$$
U^{n+1} = U^n + \frac{\Delta t}{6}(k_1 + k_2 + 4k_3)
$$

其中 $k_i = L(U^{(i-1)})$。

### 4.5 AMR 耦合

**通量回流 (Reflux)**:
- 粗网格通量：$F_{crse}$
- 细网格通量：$F_{fine}$
- 在粗细网格界面进行通量校正：
  $$
  F_{crse}^{corrected} = F_{crse} + \frac{1}{r^d}(F_{fine} - F_{crse})
  $$
  其中 $r$ 为细化比，$d$ 为空间维度。

---

## 5. 完整调用序列

```
┌────────────────────────────────────────────────────────────┐
│ 1. AmrLevelCong::advance(time, dt, iteration, ncycle)      │
│    文件：amrex_level_Cong/AmrLevelCong_advance.cpp:30      │
└─────────────┬──────────────────────────────────────────────┘
              │
              ├─ state[k].swapTimeLevels(dt)
              ├─ get_new_data(State_Type) → S_new
              ├─ get_old_data(State_Type) → S_old
              │
              └─ RK(RKOrder=3, ..., 
                    compute_dSdt<RKOrder>(stage, S, dSdt, dt_sub),
                    ...)
                   │
                   └─ [RK 阶段循环：stage=1,2,3]
                       │
                       └─ 2. compute_dSdt(stage, Sborder, dSdt, dt, ...)
                           文件：amrex_level_Cong/AmrLevelCong_advance.cpp:178
                           │
                           ├─ 2.1 守恒变量→原始变量转换
                           │   └─ cons_to_prim_for_array4(...)
                           │       文件：2phase/equation.H:176
                           │       │
                           │       └─ cons_to_prim(cons, pp)
                           │           文件：2phase/equation.H:113
                           │           │
                           │           └─ cons_to_eos<Pressure>(...)
                           │               文件：2phase/equation.H:23
                           │
                           ├─ 2.2 [对每个维度 idim=0,1,(2)]
                           │   └─ flux_solve_in_dim<idim>(params)
                           │       文件：amrex_level_Cong/AmrLevelCong_advance.cpp:395
                           │       │
                           │       ├─ 阶段 1: 节点通量计算
                           │       │   └─ prim_to_flux(prim, pp, norm)
                           │       │       文件：2phase/equation.H:200
                           │       │       │
                           │       │       └─ prims_to_eos<TotalEnergy>(...)
                           │       │           文件：2phase/equation.H:67
                           │       │
                           │       ├─ 阶段 2: Riemann 求解器
                           │       │   │
                           │       │   ├─ reconstruction<idim>(...)
                           │       │   │   文件：reconstruction_scheme/fifth_order_thinc_recon.H:98
                           │       │   │   │
                           │       │   │   └─ phi_THINC(...)
                           │       │   │       文件：fifth_order_thinc_recon.H:16
                           │       │   │
                           │       │   └─ Riemann_solver_1<idim>(...)
                           │       │       文件：2phase/equation.H:295
                           │       │       │
                           │       │       ├─ prims_to_eos<SoundSpeed>(...)
                           │       │       │   文件：2phase/equation.H:67
                           │       │       ├─ prim_to_cons(...)
                           │       │       │   文件：2phase/equation.H:140
                           │       │       └─ prim_to_flux(...)
                           │       │           文件：2phase/equation.H:200
                           │       │
                           │       └─ 阶段 3: 通量差分更新 RHS
                           │           (直接计算，无函数调用)
                           │
                           └─ [回流处理] reflux()
                               └─ FineAdd/CrseInit(...)
                           
              │
              └─ [RK 组合] LinComb/LinComb3(...) → S_new (n+1 时刻)
```

---

## 6. 关键数据结构

### 6.1 变量数组 (`VarArray`)

```cpp
using VarArray = amrex::GpuArray<amrex::Real, NVar>;
```

**变量索引顺序** (`enum` 未显示，从代码推断):
```cpp
enum StateVariables {
    AlphaRho1 = 0,   // 相 1 密度 (α₁ρ₁)
    AlphaRho2,       // 相 2 密度 (α₂ρ₂)
    XMom,            // x 动量 (ρu)
    YMom,            // y 动量 (ρv)
    ZMom,            // z 动量 (ρw) [3D]
    Energy,          // 总能量 (ρE) 或 压力 (p)【注意上下文】
    Alpha1,          // 相 1 体积分数 (α₁)
    NVar             // 变量总数 (2D: 6, 3D: 7)
};
```

### 6.2 问题参数 (`ProbParm`)

```cpp
struct ProbParm {
    // 流体属性
    amrex::Real gamma1;      // 相 1 绝热指数
    amrex::Real gamma2;      // 相 2 绝热指数
    amrex::Real pInf1;       // 相 1 参考压力
    amrex::Real pInf2;       // 相 2 参考压力
    
    // 数值参数
    amrex::Real alphaEps;    // 体积分数阈值（用于单相区判断）
    
    // 设备端成员
    device_lambda_type cal_gamma;  // 计算混合 gamma 的 lambda 函数
};
```

**混合区 gamma 计算**:
```cpp
AMREX_GPU_HOST_DEVICE
amrex::GpuArray<amrex::Real, 2> cal_gamma(amrex::Real alpha1) const
{
    // 混合规则：1/γ_bar = α₁/γ₁ + α₂/γ₂
    amrex::Real gamma_bar = 1.0 / (alpha1/gamma1 + (1-alpha1)/gamma2);
    amrex::Real pInf_bar = alpha1*pInf1 + (1-alpha1)*pInf2;
    return {gamma_bar, pInf_bar};
}
```

### 6.3 通量求解参数 (`FluxSolveParams`)

已在 3.3 节给出，此处略。

---

## 7. 数值格式说明

### 7.1 空间重构格式对比

| 格式 | 阶数 | 适用场景 | 特点 |
|------|------|----------|------|
| **一阶** | 1st | 强间断 | 稳定但耗散大 |
| **MUSCL** | 2nd | 一般流动 | 平衡精度与稳定性 |
| **WENO5** | 5th | 湍流/复杂波系 | 高精度，低耗散 |
| **THINC** | - | 多相界面 | 锐利界面捕捉 |
| **THINC-WENO** | 5th | 多相流 | 兼顾界面与波系 |

### 7.2 Riemann 求解器对比

| 求解器 | 精度 | 稳定性 | 计算成本 | 适用场景 |
|--------|------|--------|----------|----------|
| **HLL** | 低 | 高 | 低 | 强激波 |
| **HLLC** | 高 | 中 | 中 | 一般流动 |
| **Roe** | 高 | 低 | 高 | 精确捕捉接触间断 |

当前代码使用:**HLLC** + 保正限制器

### 7.3 时间积分格式

**当前使用**: 三阶 Runge-Kutta (RK3)

**Butcher 表**:
```
0   |
1   |  1
1/2 | 1/4  1/4
----|-----------
    | 1/6  1/6  2/3
```

**CFL 条件**:
$$
\Delta t = \text{CFL} \cdot \min_i \left(\frac{\Delta x_i}{|u_i| + c_i}\right)
$$

典型 CFL 数：0.3~0.6

---

## 附录 A: 主要头文件依赖关系

```
AmrLevelCong.H
├── AMReX_AmrLevel.H
├── AMReX_FArrayBox.H
├── AMReX_FluxRegister.H
├── Prob_Parm.H
└── equation.H (间接)

equation.H
├── Prob_Parm.H
├── AMReX_Box.H
├── AMReX_FArrayBox.H
└── AMReX_REAL.H

fifth_order_thinc_recon.H
├── eigen_system.H
├── interpolation_schemes.H
├── linear_algebra.H
└── positive_preserving.H
```

---

## 附录 B: GPU/CPU 异构计算

**并行策略**:
- **CPU**: OpenMP 多线程 (`#pragma omp parallel`)
- **GPU**: CUDA/HIP 核函数 (`ParallelFor` + `AMREX_GPU_DEVICE`)

**宏定义**:
```cpp
#ifdef AMREX_USE_GPU
    // GPU 代码路径
    #define AMREX_GPU_DEVICE __device__
#else
    // CPU 代码路径
    #define AMREX_GPU_DEVICE
#endif
```

---

## 附录 C: 关键参考文献

1. **Allaire 五方程模型**:  
   Allaire, G., et al. "Numerical simulation of two phase flow by a multi-component solver." (2004).

2. **THINC 格式**:  
   Xiao, F., et al. "Revisit to the THINC scheme: A simple algebraic VOF algorithm." (2011).

3. **HLLC 求解器**:  
   Toro, E.F., et al. "Restoration of the contact surface in the HLL-Riemann solver." (1994).

4. **WENO 格式**:  
   Jiang, G.S., & Shu, C.W. "Efficient implementation of weighted ENO schemes." (1996).

5. **AMReX 框架**:  
   Zhang, W., et al. "AMReX: a framework for block-structured adaptive mesh refinement." (2019).

---

**文档版本**: v1.0  
**最后更新**: 2026-03-31  
**维护者**: CFD_inWSL_MF_initial 开发团队