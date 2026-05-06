/**
 * @file
 * @brief AmrLevelCong 的时间推进、Runge-Kutta 子步和通量装配实现。
 * @ingroup amr_core
 */

#include "AmrLevelCong.H"
#include "deriv.H"
#include "equation.H"

// 【临时定位】重构方法
// 重构方法
// #include "1th_01_recon.H"      // 一阶重构方法
// #include "2th_01_muscl_recon.H"  // 二阶MUSCL重构方法
// #include "5th_01_recon.H"      // 五阶重构方法
#include "5th_02_thinc_recon.H" // 五阶THINC重构方法

// 物理约束
// #include "positive_preserving.H" // 保正性处理

// 通量计算
// #include "roe_flux.H"           // Roe通量计算

// 输出功能
#include "output_tmp.H"             // 临时输出功能

/**
 * 在此级别上推进网格。
 * 实现时间步进、通量计算和重构等核心功能
 */
using namespace amrex;

/**
 * 推进当前AMR级别的时间步
 * @param time 当前时间
 * @param dt 时间步长
 * @param iteration 迭代次数
 * @param ncycle 循环次数
 * @return 实际使用的时间步长
 */
Real AmrLevelCong::advance(Real time,
    Real dt,
    int iteration,
    int ncycle)
{
    // 打印根级别的变量最大值和最小值
    if (level == 0) {
        MultiFab& S_mm = get_new_data(State_Type);
        for (int ivar = 0; ivar < NVar; ++ivar) {
            Real maxval = S_mm.max(ivar);
            Real minval = S_mm.min(ivar);
            amrex::Print() << "Variable " << ivar
                           << ": max = " << maxval
                           << ", min = " << minval
                           << std::endl;
        }
    }
    
    // 为所有状态类型分配旧数据并交换时间水平
    for (int k = 0; k < NUM_STATE_TYPE; k++) {
        state[k].allocOldData();
        state[k].swapTimeLevels(dt);
    }

    // 获取新旧状态数据
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& S_old = get_old_data(State_Type);

    // 计算时间相关变量
    const Real prev_time = state[State_Type].prevTime();
    const Real cur_time = state[State_Type].curTime();
    const Real ctr_time = 0.5 * (prev_time + cur_time);

    // 获取计算域的左下角坐标
    GpuArray<Real, BL_SPACEDIM> prob_lo = geom.ProbLoArray();

    // 获取通量寄存器指针
    FluxRegister* fine = 0;
    FluxRegister* current = 0;

    int finest_level = parent->finestLevel();

    // 如果需要回流且当前级别不是最细级别，获取细网格的通量寄存器
    if (do_reflux && level < finest_level) {
        fine = &getFluxReg(level + 1);
        fine->setVal(0.0);
    }

    // 如果需要回流且当前级别不是根级别，获取当前级别的通量寄存器
    if (do_reflux && level > 0) {
        current = &getFluxReg(level);
    }

    // 创建带虚拟网格的状态变量和时间导数
    MultiFab Sborder(grids, dmap, NUM_STATE, NUM_GROW, MFInfo(), Factory());
    MultiFab dSdt(grids, dmap, NUM_STATE, 0, MFInfo(), Factory());

    // constexpr int RKOrder = 1;
    // FillPatch(*this, Sborder, NUM_GROW, time, State_Type, 0, NUM_STATE);
    // compute_dSdt<RKOrder>(1, Sborder, dSdt, dt, current, fine);
    // MultiFab::LinComb(S_new, 1.0, Sborder, 0, -dt, dSdt, 0, 0, NUM_STATE, 0);

    // RK2 第一阶段
    // 填充边界数据到 Sborder，用于计算状态变量的导数
    // constexpr int RKOrder = 2;
    // FillPatch(*this, Sborder, NUM_GROW, time, State_Type, 0, NUM_STATE);
    // // 计算状态变量的导数 dSdt，时间步长为 0.5 * dt
    // compute_dSdt<RKOrder>(1, Sborder, dSdt, 0.5 * dt, current, fine);
    // // 更新 S_new：S_new = Sborder - dt * dSdt
    // MultiFab::LinComb(S_new, 1.0, Sborder, 0, dt, dSdt, 0, 0, NUM_STATE, 0);

    // // RK2 第二阶段
    // // 填充边界数据到 Sborder，时间点为 time + dt
    // FillPatch(*this, Sborder, NUM_GROW, time + dt, State_Type, 0, NUM_STATE);
    // // 计算状态变量的导数 dSdt，时间步长为 0.5 * dt
    // compute_dSdt<RKOrder>(2, Sborder, dSdt, 0.5 * dt, current, fine);
    // // // 更新 S_new：S_new = 0.5 * Sborder + 0.5 * S_old
    // // MultiFab::LinComb(S_new, 0.5, Sborder, 0, 0.5, S_old, 0, 0, NUM_STATE, 0);
    // // // 进一步更新 S_new：S_new = S_new - 0.5 * dt * dSdt
    // // MultiFab::Saxpy(S_new, -0.5 * dt, dSdt, 0, 0, NUM_STATE, 0);
    // LinComb3(S_new, 0.5, Sborder, 0, 0.5, S_old, 0, 0.5 * dt, dSdt, 0, 0, NUM_STATE, IntVect(AMREX_D_DECL(0, 0, 0)));

    // constexpr int RKOrder = 3;

    // // 预先声明临时MultiFab以保持作用域可见
    // MultiFab dSdt1(grids, dmap, NUM_STATE, 0, MFInfo(), Factory());
    // MultiFab dSdt2(grids, dmap, NUM_STATE, 0, MFInfo(), Factory());
    // MultiFab dSdt3(grids, dmap, NUM_STATE, 0, MFInfo(), Factory());

    // // 阶段1：计算u^(1) = u^n + dt * L(u^n)
    // FillPatch(*this, Sborder, NUM_GROW, time, State_Type, 0, NUM_STATE);

    // compute_dSdt<RKOrder>(1, Sborder, dSdt1, dt / 6.0, current, fine); // 正确传递dt

    // // 使用LinComb直接组合（两个项）
    // MultiFab::LinComb(S_new, 1.0, Sborder, 0, -dt, dSdt1, 0, 0, NUM_STATE, 0);

    // // 阶段2：计算u^(2) = u^n + 0.25*dt*(L(u^n)+L(u^(1)))
    // FillPatch(*this, Sborder, NUM_GROW, time + dt, State_Type, 0, NUM_STATE); // 时间修正为time

    // compute_dSdt<RKOrder>(2, Sborder, dSdt2, dt / 6.0, current, fine); // 正确传递dt

    // // 使用LinComb3组合三个项（u^n + 0.25*dt*(dSdt1+dSdt2)）
    // LinComb3(S_new, 1.0, S_old, 0, // 基础项u^n
    //     -0.25 * dt, dSdt1, 0, // L(u^n)
    //     -0.25 * dt, dSdt2, 0, // L(u^(1))
    //     0, NUM_STATE, IntVect(AMREX_D_DECL(0, 0, 0)));

    // // 阶段3：计算最终解u^{n+1} = u^n + dt/6*(dSdt1 + dSdt2 + 4*dSdt3)
    // FillPatch(*this, Sborder, NUM_GROW, time + 0.5 * dt, State_Type, 0, NUM_STATE); // 时间修正为time+dt/2

    // compute_dSdt<RKOrder>(3, Sborder, dSdt3, dt * (4.0 / 6.0), current, fine); // 正确传递dt

    // // // 使用LinComb初始化为u^n，然后分步Saxpy添加各贡献项
    // // MultiFab::LinComb(S_new, 1.0, S_old, 0, 0.0, S_old, 0, 0, NUM_STATE, 0);

    // // // 添加三个阶段的贡献（避免使用LinComb3）
    // // MultiFab::Saxpy(S_new, -dt / 6.0, dSdt1, 0, 0, NUM_STATE, 0); // 第一阶段
    // // MultiFab::Saxpy(S_new, -dt / 6.0, dSdt2, 0, 0, NUM_STATE, 0); // 第二阶段
    // // MultiFab::Saxpy(S_new, -(4.0 * dt) / 6.0, dSdt3, 0, 0, NUM_STATE, 0); // 第三阶段
    // LinComb4(S_new, 1.0, S_old, 0, -dt / 6.0, dSdt1, 0, -dt / 6.0, dSdt2, 0, -(4.0 * dt) / 6.0, dSdt3, 0, 0, NUM_STATE, IntVect(AMREX_D_DECL(0, 0, 0)));

    //     void
    // AmrLevel::FillRKPatch (int state_index, MultiFab& S, Real time,
    //                        int stage, int iteration, int ncycle)

    // constexpr int RKOrder = 3;

    // // FillPatch(*this, Sborder, NUM_GROW, time, State_Type, 0, NUM_STATE);
    // compute_dSdt<RKOrder>(1, Sborder, dSdt, dt / 6.0, current, fine);
    // MultiFab::LinComb(S_new, 1.0, S_old, 0, -dt, dSdt, 0, 0, NUM_STATE, 0);

    // FillPatch(*this, Sborder, NUM_GROW, time + dt, State_Type, 0, NUM_STATE);
    // compute_dSdt<RKOrder>(2, Sborder, dSdt, dt / 6.0, current, fine);

    // LinComb3(S_new, 0.75, S_old, 0,
    //     0.25, Sborder, 0,
    //     -0.25 * dt, dSdt, 0,
    //     0, NUM_STATE, IntVect(AMREX_D_DECL(0, 0, 0)));

    // FillPatch(*this, Sborder, NUM_GROW, time + dt, State_Type, 0, NUM_STATE);
    // compute_dSdt<RKOrder>(3, Sborder, dSdt, (2.0/ 3) * dt , current, fine); // 计算L(q^(2))

    //     // 使用LinComb3组合三个项
    // LinComb3(S_new, (2.0 / 3.0), Sborder, 0,
    //     (1.0 / 3.0), S_old, 0,
    //     -(4.0 * dt) / 6.0, dSdt, 0,
    //     0, NUM_STATE, IntVect(AMREX_D_DECL(0, 0, 0)));
    // ExportAllLevelsToFiles(*this->parent);

    // 使用3阶Runge-Kutta方法进行时间步进
    constexpr int RKOrder = 3;
    RK(RKOrder, State_Type, time, dt, iteration, ncycle,
        // 计算状态变量的时间导数
        [&](int stage, MultiFab& dSdt, MultiFab const& S, Real t, Real dtsub) {
            compute_dSdt<RKOrder>(stage, S, dSdt, dtsub, current, fine);
        },
        // 每个RK子步骤后的可选操作
        [&](int /*stage*/, MultiFab& S) {});

    return dt;
}

/**
 * 获取Runge-Kutta方法的时间系数
 * @tparam RKOrder Runge-Kutta方法的阶数
 * @return 对应阶数的时间系数数组
 */
template <int RKOrder>
constexpr auto get_RK_time_coefficients()
{
    static_assert(RKOrder >= 1 && RKOrder <= 4, "RKOrder must be 1, 2, 3, or 4.");
    if constexpr (RKOrder == 1) {
        return amrex::GpuArray<amrex::Real, 1> { 1.0 };
    } else if constexpr (RKOrder == 2) {
        return amrex::GpuArray<amrex::Real, 2> { 2.0, 2.0 };
    } else if constexpr (RKOrder == 3) {
        return amrex::GpuArray<amrex::Real, 3> { 6.0, 6.0, 3.0 / 2 };
    } else if constexpr (RKOrder == 4) {
        return amrex::GpuArray<amrex::Real, 4> { 0.0, 0.5, 0.5, 1.0 };
    }
}

// 通量计算所需的虚拟网格层数
constexpr long num_flux_ghost_cells = 2;

/**
 * 计算状态变量的时间导数
 * @tparam RKOrder Runge-Kutta方法的阶数
 * @param stage 当前Runge-Kutta阶段
 * @param Sborder 带虚拟网格的状态变量
 * @param dSdt 时间导数输出
 * @param dt 时间步长
 * @param current 当前级别的通量寄存器
 * @param fine 细网格的通量寄存器
 */
template <int RKOrder>
void AmrLevelCong::compute_dSdt(int stage, const MultiFab& Sborder, MultiFab& dSdt, Real dt,
    FluxRegister* current, FluxRegister* fine)
{
    // 获取旧数据，为通量保正限制器作准备
    const MultiFab& old_data = get_old_data(State_Type);

    // 创建原始变量数组
    MultiFab prim(grids, dmap, NUM_STATE, NUM_GROW, MFInfo(), Factory());
    for (MFIter mfi(Sborder, false); mfi.isValid(); ++mfi) {
        Box ghost_box = amrex::grow(mfi.validbox(), NUM_GROW);
        Array4<Real> const& prim_array = prim[mfi].array();
        Array4<Real const> const& statein_array = Sborder[mfi].const_array();
        ProbParm const& pp = *d_prob_parm;
        // 将守恒变量转换为原始变量
        cons_to_prim_for_array4(ghost_box, prim_array, statein_array, pp);
    }

    // 初始化时间导数为0
    dSdt.setVal(0.0);
    // 创建通量数组
    MultiFab fluxes[BL_SPACEDIM];
    // 获取网格间距
    GpuArray<Real, BL_SPACEDIM> dx = geom.CellSizeArray();

    // 如果需要回流，定义通量寄存器
    if (do_reflux) {
        for (int j = 0; j < BL_SPACEDIM; j++) {
            BoxArray ba = prim.boxArray();
            ba.surroundingNodes(j);
            fluxes[j].define(ba, dmap, NUM_STATE, 0);
        }
    }
#ifdef AMREX_USE_GPU
    // 在GPU上，定义非守恒通量数组
    MultiFab nonconservfluxes[BL_SPACEDIM];
    for (int j = 0; j < BL_SPACEDIM; j++) {
        BoxArray ba = prim.boxArray();
        ba.surroundingNodes(j);
        nonconservfluxes[j].define(ba, dmap, NUM_STATE * 2, 0);
    }
#endif
    /*------------------------------------------------------------------------------*/
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        // 守恒通量
        amrex::GpuArray<FArrayBox, AMREX_SPACEDIM> fluxfab;
        amrex::GpuArray<FArrayBox*, AMREX_SPACEDIM> flux;

        // 非守恒通量
        amrex::GpuArray<FArrayBox, AMREX_SPACEDIM> nonconserfluxfab;
        amrex::GpuArray<FArrayBox*, AMREX_SPACEDIM> nonconserflux;

        for (MFIter mfi(prim, false); mfi.isValid(); ++mfi) {
            // 设置 tileboxes 和 nodal tileboxes
            Box bx = mfi.tilebox();
            Box fluxBx = amrex::grow(bx, num_flux_ghost_cells);
            GpuArray<Box, BL_SPACEDIM> nbx;
            AMREX_D_TERM(nbx[0] = mfi.nodaltilebox(0);,
                nbx[1] = mfi.nodaltilebox(1);
                ,
                nbx[2] = mfi.nodaltilebox(2));

            const FArrayBox& primin = prim[mfi];
            const FArrayBox& consin = Sborder[mfi];
            FArrayBox& rhs = dSdt[mfi];

            // 遍历所有空间维度
            for (int idim = 0; idim < BL_SPACEDIM; idim++) {
#ifdef AMREX_USE_GPU
                // 在 GPU 上不使用分块，直接将通量和面速度的 fab 指针指向未分块的 fab
                flux[idim] = &(fluxes[idim][mfi]);
                nonconserflux[idim] = &(nonconservfluxes[idim][mfi]);
#else
                // 在 CPU 上，重新调整临时 fab 的大小以存储通量和面速度
                const Box& bxflux = amrex::surroundingNodes(bx, idim);
                fluxfab[idim].resize(bxflux, NUM_STATE);

                // 将通量和面速度的 fab 指针指向临时 fab
                flux[idim] = &(fluxfab[idim]);

                // 处理非守恒通量，只在mfi的局部使用
                nonconserfluxfab[idim].resize(bxflux, NUM_STATE * 2);
                nonconserflux[idim] = &(nonconserfluxfab[idim]);
#endif
            }
            // 以alpha*du/dx为例，这里的cellNcFlux指代的是将用于差分的u
            // 因此有多少项这种输运项就需要多少个组件
            FArrayBox cellFlux(fluxBx, NUM_STATE), cellNcFlux(fluxBx, AdvectionTerm::num_terms);
            // 获取当前Runge-Kutta阶段的时间系数
            Real RKTimeCoef = get_RK_time_coefficients<RKOrder>()[stage - 1];
            // 使用 AMREX_D_TERM 宏遍历所有空间维度，计算各维度的通量
            AMREX_D_TERM(
                flux_solve_in_dim<0>({ bx, flux, nonconserflux, cellFlux, cellNcFlux, primin, consin, rhs, nbx, dt, RKTimeCoef });,
                flux_solve_in_dim<1>({ bx, flux, nonconserflux, cellFlux, cellNcFlux, primin, consin, rhs, nbx, dt, RKTimeCoef });
                ,
                flux_solve_in_dim<2>({ bx, flux, nonconserflux, cellFlux, cellNcFlux, primin, consin, rhs, nbx, dt, RKTimeCoef });)

            // 通量保正限制器（当前被注释掉）
            // 获取当前 Runge - Kutta 阶段的时间系数
            // Real RKTimeCoef = get_RK_time_coefficients<RKOrder>()[stage - 1];

            // ProbParm const& pp = *d_prob_parm;

            // //     // 定义正限制函数对象
            // const PositiveQuantityFunctions pqf;
            // const FArrayBox& old_data_fab = old_data[mfi];

            // apply_positivive_preserving_flux_limiter(
            //     bx, // 计算区域
            //     Sborder[mfi], // 守恒量数组
            //     primin, // 原始量数组
            //     old_data_fab,
            //     NUM_STATE, // 守恒量的分量数
            //     rhs, // 右端项数组
            //     flux, // 通量数组
            //     nonconserflux, // 非守恒通量数组
            //     dt * RKTimeCoef, // RK凸组合内的时间步长系数
            //     stage, // Runge-Kutta 阶段
            //     pp, // 问题参数
            //     pqf, // 正限制函数
            //     dx // 网格间距
            // );

/*------------------------------------------------------------------------------*/
#ifndef AMREX_USE_GPU
            // 在CPU上，如果需要回流，将通量复制到通量寄存器
            if (do_reflux) {
                for (int i = 0; i < BL_SPACEDIM; i++)
                    fluxes[i][mfi].copy(*flux[i], mfi.nodaltilebox(i));
            }
#endif
        }
    }

    // 处理回流
    if (do_reflux) {
        // 计算每个维度的通量缩放因子
        GpuArray<Real, BL_SPACEDIM> fluxScalingFactors;

        for (int idim = 0; idim < BL_SPACEDIM; idim++) {
            Real scalingFactor = dt; // 临时变量存储时间步长
            for (int iidim = 0; iidim < BL_SPACEDIM; iidim++) {
                if (iidim != idim) 
                    scalingFactor *= dx[iidim]; // 在非当前方向上乘以网格间距
            }
            fluxScalingFactors[idim] = scalingFactor; // 将结果存储到 GpuArray 中
        }

        // 如果存在当前级别的通量寄存器，添加通量
        if (current) {
            for (int i = 0; i < BL_SPACEDIM; i++)
                current->FineAdd(fluxes[i], i, 0, 0, NUM_STATE, fluxScalingFactors[i]);
        }
        // 如果存在细网格的通量寄存器，初始化通量
        if (fine) {
            for (int i = 0; i < BL_SPACEDIM; i++)
                fine->CrseInit(fluxes[i], i, 0, 0, NUM_STATE, -fluxScalingFactors[i], FluxRegister::ADD);
        }
    }
}

/**
 * 在指定维度上计算通量
 * @tparam idim 空间维度索引
 * @param params 通量计算参数
 */
template <int idim>
void AmrLevelCong::flux_solve_in_dim(const FluxSolveParams& params)
{
    // 解包参数
    const auto& bx = params.bx;
    auto& flux = params.flux;
    auto& ncflux = params.ncflux;
    auto& cellFlux = params.cellFlux;
    auto& cellNcFlux = params.cellNcFlux;
    const auto& primin = params.primin;
    const auto& consin = params.consin;
    auto& rhs = params.rhs;
    const auto& nbx = params.nbx;
    const auto& dt = params.dt;
    const auto& RKTimeCoef = params.RKTimeCoef;

    // 断言通量计算区域包含计算盒
    AMREX_ASSERT(cellFlux.box().contains(bx));
    AMREX_ASSERT(cellNcFlux.box().contains(bx));

    // 维度偏移量（用于跨维度访问）
    Dim3 offset = amrex::IntVect::TheDimensionVector(idim).dim3();

    // 获取问题参数和几何信息
    ProbParm const& pp = *d_prob_parm;
    amrex::GpuArray<amrex::Real, 3> norm = get_direction<idim>();
    GpuArray<Real, BL_SPACEDIM> dx = geom.CellSizeArray();

    // 扩展通量计算区域（包含虚拟网格）
    Box fluxBx = amrex::grow(bx, num_flux_ghost_cells);

    // 初始化节点通量存储
    cellFlux.setVal<RunOn::Device>(0.0);

    // 初始化非守恒通量
    cellNcFlux.setVal<RunOn::Device>(0.0);
    ncflux[idim]->setVal<RunOn::Device>(0.0);

    // 获取数组视图
    Array4<Real const> const& primin_array = primin.const_array();
    Array4<Real const> const& consin_array = consin.const_array();
    Array4<Real> const& fluxRef = flux[idim]->array();
    Array4<Real> const& stateout_array = rhs.array();
    Array4<Real> const& nodeFlux_array = cellFlux.array();

    // 非守恒通量数组视图
    Array4<Real> const& nodencFlux_array = cellNcFlux.array();
    Array4<Real> const& ncfluxRef = ncflux[idim]->array();

    // ----------------------------
    // 阶段1：计算节点通量
    // ----------------------------
    ParallelFor(fluxBx, [=] AMREX_GPU_DEVICE(int xIdx, int yIdx, int zIdx) {
        // 提取原始状态变量
        VarArray primitiveState;
        for (int stateIdx = 0; stateIdx < NUM_STATE; stateIdx++) {
            primitiveState[stateIdx] = primin_array(xIdx, yIdx, zIdx, stateIdx);
        }

        // 计算守恒通量
        VarArray fluxState = prim_to_flux(primitiveState, pp, norm);
        for (int stateIdx = 0; stateIdx < NUM_STATE; stateIdx++) {
            nodeFlux_array(xIdx, yIdx, zIdx, stateIdx) = fluxState[stateIdx];
        }

        // 计算非守恒项通量
        AdvectionTerm advTerm;
        auto flux_var = advTerm.computeNodeFlux(primitiveState, pp, norm);
        for (int idx = 0; idx < AdvectionTerm::num_terms; idx++) {
            nodencFlux_array(xIdx, yIdx, zIdx, idx) += flux_var[idx];
        }
    });

    // ----------------------------
    // 阶段2：Riemann求解器计算通量
    // ----------------------------
    ParallelFor(nbx[idim], [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        // 重构左右状态
        amrex::GpuArray<VarArray, 2> stateLR = reconstruction<idim>(i, j, k, primin_array, NUM_STATE, pp, norm, dx);

        // 使用 Riemann 求解器计算通量
        VarArray iFlux = Riemann_solver_1<idim>(i, j, k, primin_array, stateLR[0], stateLR[1], norm, pp, dx);

        // 对计算出的通量进行修正，使用四阶精度的中心差分格式
        VarArray consL, consR, fluxCellL, fluxCellR;
        for (int ivar = 0; ivar < NUM_STATE; ivar++) {
            consL[ivar] = consin_array(i - offset.x, j - offset.y, k - offset.z, ivar); // 左侧单元的守恒变量
            consR[ivar] = consin_array(i, j, k, ivar); // 右侧单元的守恒变量
            fluxCellL[ivar] = nodeFlux_array(i - offset.x, j - offset.y, k - offset.z, ivar); // 左侧单元的通量
            fluxCellR[ivar] = nodeFlux_array(i, j, k, ivar); // 右侧单元的通量
        }
        
        // 通量修正
        for (int ivar = 0; ivar < NUM_STATE; ivar++) {
            iFlux[ivar] = 4.0 / 3.0 * iFlux[ivar]
                - 1.0 / 6.0 * (nodeFlux_array(i - offset.x, j - offset.y, k - offset.z, ivar) + nodeFlux_array(i, j, k, ivar));
        }

        // 非守恒通量修正
        AdvectionTerm advTerm; // 实例化非守恒项对象
        // 计算 alpha du/dx 中半节点处的 u_1/2
        auto flux_var = advTerm.computeFlux(stateLR[0], stateLR[1], pp, norm, dx);
        // 获取输运项对应的方程索引
        auto indices = advTerm.getTargetTermIndices();

        // 计算 alpha du/dx 节点处的 alpha_i 和 alpha_i+1
        VarArray primL, primR;
        for (int ivar = 0; ivar < NUM_STATE; ivar++) {
            primL[ivar] = primin_array(i - offset.x, j - offset.y, k - offset.z, ivar); // 左侧单元的原始变量
            primR[ivar] = primin_array(i, j, k, ivar); // 右侧单元的原始变量
        }
        // 使用 AdvectionTerm 计算单元中心处的系数 alpha
        auto coefL = advTerm.computeCoefficient(primL, pp, norm);
        auto coefR = advTerm.computeCoefficient(primR, pp, norm);

        // 初始化增量通量数组
        VarArray incFluxL, incFluxR;
        incFluxL.fill(0.0);
        incFluxR.fill(0.0);
        
        // 根据非守恒项计算增量通量
        for (int idx = 0; idx < AdvectionTerm::num_terms; idx++) {
            incFluxL[indices[idx]] += coefL[idx]
                * (4.0 / 3.0 * flux_var[idx]
                    - 1.0 / 6.0 * (nodencFlux_array(i - offset.x, j - offset.y, k - offset.z, idx) + nodencFlux_array(i, j, k, idx)));
            incFluxR[indices[idx]] += coefR[idx]
                * (4.0 / 3.0 * flux_var[idx]
                    - 1.0 / 6.0 * (nodencFlux_array(i - offset.x, j - offset.y, k - offset.z, idx) + nodencFlux_array(i, j, k, idx)));
        }

        // 计算时间步长相关的乘数因子
        Real multL = 2. * dt * RKTimeCoef / dx[idim] / computeAlpha(dx, primL, pp)[idim];
        Real multR = -2. * dt * RKTimeCoef / dx[idim] / computeAlpha(dx, primR, pp)[idim];

        // 计算中间状态的守恒变量
        VarArray consStarL, consStarR;
        PositiveQuantityFunctions pqf;
        for (int ivar = 0; ivar < NUM_STATE; ivar++) {
            consStarL[ivar] = consL[ivar] - multL * (iFlux[ivar] + incFluxL[ivar] - fluxCellL[ivar]);
            consStarR[ivar] = consR[ivar] - multR * (iFlux[ivar] + incFluxR[ivar] - fluxCellR[ivar]);
        }

        // 如果中间状态包含无效值，则应用凸组合限制器
        bool leftFlag = containsInvalidValues(consStarL, pp, pqf);
        bool rightFlag = containsInvalidValues(consStarR, pp, pqf);
        if (leftFlag || rightFlag) {
            VarArray foFlux = Riemann_solver_2(i, j, k, primin_array, primL, primR, norm, pp), foncFluxL, foncFluxR;
            foncFluxL.fill(0.0);
            foncFluxR.fill(0.0);
            auto foflux_var = advTerm.computeFlux(primL, primR, pp, norm, dx);
            for (int idx = 0; idx < AdvectionTerm::num_terms; idx++) {
                foncFluxL[indices[idx]] += coefL[idx] * foflux_var[idx];
                foncFluxR[indices[idx]] += coefR[idx] * foflux_var[idx];
            }
            // 计算限制后的守恒变量
            VarArray consFoL, consFoR;
            for (int ivar = 0; ivar < NUM_STATE; ivar++) {
                consFoL[ivar] = consL[ivar] - multL * (foFlux[ivar] + foncFluxL[ivar] - fluxCellL[ivar]);
                consFoR[ivar] = consR[ivar] - multR * (foFlux[ivar] + foncFluxR[ivar] - fluxCellR[ivar]);
            }
            // 使用 Lambda 表达式计算 thetaL 和 thetaR
            auto computeThetaWithCondition = [](bool flag, Real (*computeThetaFunc)(const VarArray&, const VarArray&, const ProbParm&, const PositiveQuantityFunctions&),
                                                 const VarArray& consStar, const VarArray& consFo, const ProbParm& pp, const PositiveQuantityFunctions& pqf) -> Real {
                return flag ? computeThetaFunc(consStar, consFo, pp, pqf) : 1.0;
            };
            Real thetaL = computeThetaWithCondition(leftFlag, compute_theta, consStarL, consFoL, pp, pqf);
            Real thetaR = computeThetaWithCondition(rightFlag, compute_theta, consStarR, consFoR, pp, pqf);
            Real theta = min(thetaL, thetaR);
            apply_convex_combination(foFlux, iFlux, theta, 0, NVar);
            apply_convex_combination(foncFluxL, incFluxL, theta, 0, NVar);
            apply_convex_combination(foncFluxR, incFluxR, theta, 0, NVar);
        }

        // 更新通量和非守恒通量数组
        for (int ivar = 0; ivar < NVar; ivar++) {
            fluxRef(i, j, k, ivar) = iFlux[ivar]; // 更新守恒通量
            ncfluxRef(i, j, k, ivar) = incFluxL[ivar]; // 更新左侧非守恒通量
            ncfluxRef(i, j, k, ivar + NVar) = incFluxR[ivar]; // 更新右侧非守恒通量
        }
    });

    // ----------------------------
    // 阶段3：更新右手边（RHS）
    // ----------------------------
    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        for (int ivar = 0; ivar < NUM_STATE; ivar++) {
            // 守恒项通量差分
            Real flux_diff = (fluxRef(i + offset.x, j + offset.y, k + offset.z, ivar) - fluxRef(i, j, k, ivar)) / dx[idim];
            stateout_array(i, j, k, ivar) -= flux_diff;

            // 非守恒项差分
            // 注意要用右面的左值减去左面的右值
            Real ncflux_diff = (ncfluxRef(i + offset.x, j + offset.y, k + offset.z, ivar) - ncfluxRef(i, j, k, ivar + NUM_STATE)) / dx[idim];
            stateout_array(i, j, k, ivar) -= ncflux_diff;
        }
    });
}