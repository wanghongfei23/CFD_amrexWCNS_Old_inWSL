/**
 * @file
 * @brief AmrLevelCong 的静态成员、参数读取与层级生命周期实现。
 * @ingroup amr_core
 */

// #include <AMReX_GpuMemory.H>
#include <AMReX_ParmParse.H>
// #include <AMReX_TagBox.H>
#include <AMReX_VisMF.H>

#include "AmrLevelCong.H"
#include "equation.H"

using namespace amrex;

// 类静态成员初始化
int AmrLevelCong::verbose = 0; //  verbose level
Real AmrLevelCong::cfl = 0.5; // CFL数
int AmrLevelCong::do_reflux = 0; // 是否进行通量回流

int AmrLevelCong::NUM_STATE = 3; // 状态变量数量
int AmrLevelCong::NUM_GROW = 3; // 鬼细胞数量

ProbParm* AmrLevelCong::h_prob_parm = nullptr; // 主机端问题参数
ProbParm* AmrLevelCong::d_prob_parm = nullptr; // 设备端问题参数

int AmrLevelCong::max_phierr_lev = -1; // 最大phi误差级别
int AmrLevelCong::max_phigrad_lev = -1; // 最大phi梯度级别

Vector<Real> AmrLevelCong::phierr; // phi误差数组
Vector<Real> AmrLevelCong::phigrad; // phi梯度数组

/**
 * 默认构造函数。构建无效对象。
 */
AmrLevelCong::AmrLevelCong()
{
    flux_reg = 0;
}

/**
 * 基本构造函数。
 */
AmrLevelCong::AmrLevelCong(Amr& papa,
    int lev,
    const Geometry& level_geom,
    const BoxArray& bl,
    const DistributionMapping& dm,
    Real time)
    : AmrLevel(papa, lev, level_geom, bl, dm, time)
{
    flux_reg = 0;
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids, dmap, crse_ratio, level, NUM_STATE);
}

/**
 * 析构函数。
 */
AmrLevelCong::~AmrLevelCong()
{
    delete flux_reg;
}

/**
 * 从检查点文件重启。
 */
void AmrLevelCong::restart(Amr& papa,
    std::istream& is,
    bool bReadSpecial)
{
    AmrLevel::restart(papa, is, bReadSpecial);

    BL_ASSERT(flux_reg == 0);
    if (level > 0 && do_reflux)
        flux_reg = new FluxRegister(grids, dmap, crse_ratio, level, NUM_STATE);
}

/**
 * 写入检查点文件。
 */
void AmrLevelCong::checkPoint(const std::string& dir,
    std::ostream& os,
    VisMF::How how,
    bool dump_old)
{
    AmrLevel::checkPoint(dir, os, how, dump_old);
}

/**
 * 写入绘图文件到指定目录。
 */
void AmrLevelCong::writePlotFile(const std::string& dir,
    std::ostream& os,
    VisMF::How how)
{

    AmrLevel::writePlotFile(dir, os, how);
}

/**
 * 在问题启动时初始化网格数据。
 */
void AmrLevelCong::initData()
{
    //
    // 遍历网格。
    //
    const Real* dx = geom.CellSize();
    const Real* prob_lo = geom.ProbLo();
    MultiFab& S_new = get_new_data(State_Type);
    Real cur_time = state[State_Type].curTime();

    const ProbParm& pp = *d_prob_parm;

    if (verbose) {
        amrex::Print() << "Initializing the data at level " << level << std::endl;
    }

#ifdef AMREX_USE_GPU
    // 在CPU上使用固定内存创建临时MultiFab
    MultiFab S_tmp(S_new.boxArray(),
        S_new.DistributionMap(),
        S_new.nComp(),
        S_new.nGrowVect(),
        MFInfo().SetArena(The_Pinned_Arena()));
#else
    // 使用MultiFab指针
    MultiFab& S_tmp = S_new;
#endif

    for (MFIter mfi(S_tmp); mfi.isValid(); ++mfi) {
        const Box& box = mfi.validbox();
        // auto dx = geom.CellSize();
        const int* lo = box.loVect();
        const int* hi = box.hiVect();

        // 使用Fortran子程序在CPU上初始化数据。
        // initdata(&level, &cur_time, AMREX_ARLIM_3D(lo), AMREX_ARLIM_3D(hi),
        //     BL_TO_FORTRAN_3D(S_tmp[mfi]), AMREX_ZFILL(dx),
        //     AMREX_ZFILL(prob_lo));
        //     void initdata(int level, amrex::Real time,
        // const amrex::Box& box, amrex::FArrayBox& state,
        // const amrex::Real* dx, const amrex::Real* problo)
        initdata(level, cur_time, box, S_tmp[mfi], dx, geom.ProbLo(), pp);
    }

#ifdef AMREX_USE_GPU
    // 显式将数据复制到GPU。
    amrex::htod_memcpy(S_new, S_tmp);
#endif

    if (verbose) {
        amrex::Print() << "Done initializing the level " << level
                       << " data " << std::endl;
    }
}

/**
 * 在重网格期间从另一个AmrLevelCong初始化此级别的数据。
 */
void AmrLevelCong::init(AmrLevel& old)
{
    AmrLevelCong* oldlev = (AmrLevelCong*)&old;

    //
    // 通过从旧网格填充来创建新网格数据。
    //
    Real dt_new = parent->dtLevel(level);
    Real cur_time = oldlev->state[State_Type].curTime();
    Real prev_time = oldlev->state[State_Type].prevTime();
    Real dt_old = cur_time - prev_time;
    setTimeLevel(cur_time, dt_old, dt_new);

    MultiFab& S_new = get_new_data(State_Type);

    FillPatch(old, S_new, 0, cur_time, State_Type, 0, NUM_STATE);
}

/**
 * 如果旧级别之前不存在，则在重网格后初始化此级别的数据
 */
void AmrLevelCong::init()
{
    Real dt = parent->dtLevel(level);
    Real cur_time = getLevel(level - 1).state[State_Type].curTime();
    Real prev_time = getLevel(level - 1).state[State_Type].prevTime();

    Real dt_old = (cur_time - prev_time) / (Real)parent->MaxRefRatio(level - 1);

    setTimeLevel(cur_time, dt_old, dt);
    MultiFab& S_new = get_new_data(State_Type);
    FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, NUM_STATE);
}

/**
 * 估计时间步长。
 */
Real AmrLevelCong::estTimeStep(Real)
{
    // 这只是一个初始的虚拟值
    ProbParm const& pp = *d_prob_parm;
    Real dt_est = 1.0e+20;

    GpuArray<Real, BL_SPACEDIM> dx = geom.CellSizeArray();
    GpuArray<Real, BL_SPACEDIM> prob_lo = geom.ProbLoArray();
    const Real cur_time = state[State_Type].curTime();
    const MultiFab& S_new = get_new_data(State_Type);
    Real pred_time = cur_time;
    if (cur_time > 0._rt) {
        pred_time += 0.5_rt * parent->dtLevel(level);
    }

    // 初始化 tau_xyz 为零
    for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
        tau_xyz[dir] = 0.0;
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(min : dt_est)
#endif
    {
        FArrayBox umaxFAB; // 用于存储局部最大特征速度的FArrayBox
        for (MFIter mfi(S_new, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox(); // 当前网格块的范围
            umaxFAB.resize(bx, BL_SPACEDIM); // 每个网格点存储BL_SPACEDIM个特征速度
            auto umaxArray = umaxFAB.array(); // 获取FArrayBox的数组

            // 获取密度、动量和总能量的数组视图
            Array4<Real const> const& state_array = S_new[mfi].const_array();

            // 并行计算每个网格点的局部最大特征速度
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                cal_u_max(i, j, k, state_array, umaxArray, pp, dx);
            });

            // 计算当前网格块的最大特征速度
            GpuArray<Real, BL_SPACEDIM> local_tau_xyz; // 局部最大特征速度
            for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
                local_tau_xyz[dir] = umaxFAB.norm<RunOn::Device>(0, dir, 1);
            }

            // 更新全局最大特征速度
            for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
                tau_xyz[dir] = amrex::max(tau_xyz[dir], local_tau_xyz[dir]);
            }
        }
    }

    // 归约操作：确保 tau_xyz 在所有进程中一致
    for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
        ParallelDescriptor::ReduceRealMax(tau_xyz[dir]);
    }

    // 计算时间步长估计值（放在循环外）
    Real norminator = 0.0;
    for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
        norminator += tau_xyz[dir];
    }
    if (norminator > 1.e-100) {
        dt_est = amrex::min(dt_est, 1.0 / norminator);
    }

    // 应用 CFL 条件
    dt_est *= cfl;

    if (verbose) {
        amrex::Print() << "AmrLevelCong::estTimeStep at level " << level
                       << ":  dt_est = " << dt_est << std::endl;
        amrex::Print() << "AmrLevelCong::estTimeStep: tau_xyz = ";
        for (int dir = 0; dir < BL_SPACEDIM; ++dir) {
            amrex::Print() << tau_xyz[dir] << " ";
        }
        amrex::Print() << std::endl;
    }

    return dt_est;
}

/**
 * 计算初始时间步长。
 */
Real AmrLevelCong::initialTimeStep()
{
    return estTimeStep(0.0);
}

/**
 * @brief 计算AMR层次结构中所有级别的初始时间步长。
 *
 * 此函数计算模拟中每个级别的初始时间步长（`dt`），
 * 确保所有级别之间的一致性，并在指定时尊重停止时间。
 *
 * @param finest_level AMR层次结构中的最细级别。
 * @param n_cycle 包含每个级别的子循环数的向量。
 * @param dt_level 用于存储每个级别的计算时间步长的向量。
 * @param stop_time 允许的最大模拟时间。调整`dt`以避免超过此值。
 */
void AmrLevelCong::computeInitialDt(int finest_level,
    int /*sub_cycle*/, // 未使用的参数
    Vector<int>& n_cycle,
    const Vector<IntVect>& /*ref_ratio*/, // 未使用的参数
    Vector<Real>& dt_level,
    Real stop_time)
{
    // 此函数应仅在最粗级别（level == 0）执行。
    // 如果当前级别大于0，立即退出函数。
    if (level > 0)
        return;

    // 初始化`dt_0`为一个大值，以找到所有级别中的最小时间步长。
    Real dt_0 = 1.0e+100;

    // `n_factor` 累积级别之间子循环的缩放因子。
    int n_factor = 1;

    // 从最粗级别到最细级别遍历所有级别。
    for (int i = 0; i <= finest_level; i++) {
        // 使用其特定逻辑计算当前级别的初始时间步长。
        dt_level[i] = getLevel(i).initialTimeStep();

        // 更新当前级别的子循环因子。
        n_factor *= n_cycle[i];

        // 找到所有级别中的最小缩放时间步长。
        dt_0 = amrex::min(dt_0, n_factor * dt_level[i]);
    }

    // 调整计算的`dt_0`以确保它不超过停止时间。
    const Real eps = 0.001 * dt_0; // 小公差以避免稍微超过stop_time。
    Real cur_time = state[State_Type].curTime(); // 当前模拟时间。
    if (stop_time >= 0.0) {
        // 如果应用`dt_0`后的总时间超过stop_time，调整`dt_0`。
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    // 将调整后的`dt_0`重新分配到所有级别。
    n_factor = 1; // 重置子循环因子。
    for (int i = 0; i <= finest_level; i++) {
        // 更新当前级别的子循环因子。
        n_factor *= n_cycle[i];

        // 为当前级别分配时间步长，按子循环因子缩放。
        dt_level[i] = dt_0 / n_factor;
    }
}

/**
 * 计算新的`dt'。
 */
void AmrLevelCong::computeNewDt(int finest_level,
    int /*sub_cycle*/,
    Vector<int>& n_cycle,
    const Vector<IntVect>& /*ref_ratio*/,
    Vector<Real>& dt_min,
    Vector<Real>& dt_level,
    Real stop_time,
    int post_regrid_flag)
{
    //
    // 我们在粗网格时间循环的末尾。
    // 计算下一次迭代的时间步长。
    //
    if (level > 0)
        return;

    for (int i = 0; i <= finest_level; i++) {
        AmrLevelCong& adv_level = getLevel(i);
        dt_min[i] = adv_level.estTimeStep(dt_level[i]);
    }

    if (post_regrid_flag == 1) {
        //
        // 按重网格前的dt限制dt
        //
        for (int i = 0; i <= finest_level; i++) {
            dt_min[i] = amrex::min(dt_min[i], dt_level[i]);
        }
    } else {
        //
        // 按change_max * 旧dt限制dt
        //
        static Real change_max = 10;
        for (int i = 0; i <= finest_level; i++) {
            dt_min[i] = amrex::min(dt_min[i], change_max * dt_level[i]);
        }
    }

    //
    // 找到所有级别的最小值
    //
    Real dt_0 = 1.0e+100;
    int n_factor = 1;
    for (int i = 0; i <= finest_level; i++) {
        n_factor *= n_cycle[i];
        dt_0 = amrex::min(dt_0, n_factor * dt_min[i]);
    }

    //
    // 按stop_time的值限制dt。
    //
    const Real eps = 0.001 * dt_0;
    Real cur_time = state[State_Type].curTime();
    if (stop_time >= 0.0) {
        if ((cur_time + dt_0) > (stop_time - eps))
            dt_0 = stop_time - cur_time;
    }

    n_factor = 1;
    for (int i = 0; i <= finest_level; i++) {
        n_factor *= n_cycle[i];
        dt_level[i] = dt_0 / n_factor;
    }
}

/**
 * 在timestep()之后工作。
 */
void AmrLevelCong::post_timestep(int iteration)
{
    //
    // 细级别网格的积分循环完成
    // 在这里做post_timestep的工作。
    //
    int finest_level = parent->finestLevel();

    if (do_reflux && level < finest_level)
        reflux();

    if (level < finest_level)
        avgDown();
}

/**
 * 在regrid()之后工作。
 */
void AmrLevelCong::post_regrid(int lbase, int /*new_finest*/)
{
    amrex::ignore_unused(lbase);
}

/**
 * 在restart()之后工作。
 */
void AmrLevelCong::post_restart()
{
}

/**
 * 在init()之后工作。
 */
void AmrLevelCong::post_init(Real /*stop_time*/)
{
    if (level > 0)
        return;
    //
    // 从更细的级别向下平均数据
    // 使守恒数据在级别之间一致。
    //
    int finest_level = parent->finestLevel();
    for (int k = finest_level - 1; k >= 0; k--)
        getLevel(k).avgDown();
}

/**
 * 从输入文件读取参数。
 */
void AmrLevelCong::read_params()
{
    static bool done = false;

    if (done)
        return;

    done = true;

    ParmParse pp("adv");

    pp.query("v", verbose);
    pp.query("cfl", cfl);
    pp.query("do_reflux", do_reflux);
    Geometry const* gg
        = AMReX::top()->getDefaultGeometry();

    // 本教程代码仅支持笛卡尔坐标。
    if (!gg->IsCartesian()) {
        amrex::Abort("Please set geom.coord_sys = 0");
    }

    // 本教程代码仅支持周期性边界。
    // if (!gg->isAllPeriodic()) {
    //     amrex::Abort("Please set geom.is_periodic = 1 1 1");
    // }

    // 从输入文件的tagging块读取标记参数。
    // 参见Src_nd/Tagging_params.cpp中的函数实现。
    get_tagging_params();
}

// 通量回流函数
void AmrLevelCong::reflux()
{
    BL_ASSERT(level < parent->finestLevel());

    const auto strt = amrex::second();

    getFluxReg(level + 1).Reflux(get_new_data(State_Type), 1.0, 0, 0, NUM_STATE, geom);

    if (verbose) {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();
        auto end = amrex::second() - strt;

        ParallelDescriptor::ReduceRealMax(end, IOProc);

        amrex::Print() << "AmrLevelCong::reflux() at level " << level
                       << " : time = " << end << std::endl;
    }
}

// 向下平均函数
void AmrLevelCong::avgDown()
{
    if (level == parent->finestLevel())
        return;
    avgDown(State_Type);
}

// 向下平均函数（带状态索引）
void AmrLevelCong::avgDown(int state_indx)
{
    if (level == parent->finestLevel())
        return;

    AmrLevelCong& fine_lev = getLevel(level + 1);
    MultiFab& S_fine = fine_lev.get_new_data(state_indx);
    MultiFab& S_crse = get_new_data(state_indx);

    amrex::average_down(S_fine, S_crse,
        fine_lev.geom, geom,
        0, S_fine.nComp(), parent->refRatio(level));
}
