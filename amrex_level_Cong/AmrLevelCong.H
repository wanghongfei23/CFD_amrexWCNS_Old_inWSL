/**
 * @file
 * @brief AMR 层级对象 AmrLevelCong 的接口定义。
 * @details 该头文件集中声明 AMR 生命周期、时间推进、回流、误差标记和边界辅助接口。
 * @ingroup amr_core
 */

#ifndef _AmrLevelAdv_H_
#define _AmrLevelAdv_H_

#include <AMReX_AmrLevel.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_PhysBCFunct.H>

#include <iostream>
#include <memory>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

// 问题相关结构体。在Exec目录下定义。
#include "Prob_Parm.H"

/**
 * @enum StateType
 * @brief AMReX 状态描述符中的数据槽编号。
 */
enum StateType { 
    State_Type = 0,  // 状态类型
    NUM_STATE_TYPE   // 状态类型数量
};

/**
 * @struct FluxSolveParams
 * @brief 单个方向通量离散所需的临时数据视图。
 */
struct FluxSolveParams {
    const amrex::Box& bx;                                  // 计算区域
    amrex::GpuArray<amrex::FArrayBox*, AMREX_SPACEDIM>& flux; // 通量数组
    amrex::GpuArray<amrex::FArrayBox*, AMREX_SPACEDIM>& ncflux; // 非守恒通量数组
    amrex::FArrayBox& cellFlux;                           // 单元格通量，大小=NUM_STATE
    amrex::FArrayBox& cellNcFlux;                         // 单元格非守恒通量，大小=AdvectionTerm::num_terms
    const amrex::FArrayBox& primin;                       // 原始变量数组
    const amrex::FArrayBox& consin;                       // 守恒变量数组
    amrex::FArrayBox& rhs;                                // 右端项数组
    const amrex::GpuArray<amrex::Box, BL_SPACEDIM>& nbx;   // 邻居边界框
    amrex::Real dt;                                       // 时间步长
    amrex::Real RKTimeCoef;                               // Runge-Kutta时间系数
};

/**
 * @brief 基于 amrex::AmrLevel 的双曲守恒律层级对象。
 * @details 该类负责管理状态变量、时间推进、回流、AMR 初始化与误差标记。
 */
class AmrLevelCong
    : public amrex::AmrLevel {
public:
    /**
     * 默认构造函数。构建无效对象。
     */
    AmrLevelCong();

    /**
     * 基本构造函数。
     */
    AmrLevelCong(amrex::Amr& papa,        // 父Amr对象
        int lev,                           // 级别
        const amrex::Geometry& level_geom,  // 级别几何信息
        const amrex::BoxArray& bl,         // 网格数组
        const amrex::DistributionMapping& dm, // 分布映射
        amrex::Real time);                 // 时间

    /**
     * 析构函数。
     */
    virtual ~AmrLevelCong() override;

    /**
     * 从检查点文件重启。
     */
    virtual void restart(amrex::Amr& papa,    // 父Amr对象
        std::istream& is,                     // 输入流
        bool bReadSpecial = false) override;  // 是否读取特殊数据

    /**
     * 写入检查点文件。
     */
    virtual void checkPoint(const std::string& dir,           // 目录
        std::ostream& os,                                     // 输出流
        amrex::VisMF::How how = amrex::VisMF::NFiles,          // 写入方式
        bool dump_old = true) override;                        // 是否转储旧数据

    /**
     * 写入绘图文件到指定目录。
     */
    virtual void writePlotFile(const std::string& dir,     // 目录
        std::ostream& os,                                 // 输出流
        amrex::VisMF::How how) override;                  // 写入方式

    /**
     * 定义数据描述符。
     */
    static void variableSetUp();

    /**
     * 在运行结束时清理数据描述符。
     */
    static void variableCleanUp();

    /**
     * 在问题启动时初始化网格数据。
     */
    virtual void initData() override;

    /**
     * 在重网格期间从另一个AmrLevelCong初始化此级别的数据。
     */
    virtual void init(amrex::AmrLevel& old) override;  // 旧AmrLevel对象

    /**
     * 如果旧级别之前不存在，则在重网格后初始化此级别的数据
     */
    virtual void init() override;

    /**
     * 在此级别上随时间推进网格。
     */
    virtual amrex::Real advance(amrex::Real time,    // 时间
        amrex::Real dt,                             // 时间步长
        int iteration,                              // 迭代次数
        int ncycle) override;                       // 循环次数

    /**
     * @brief 计算 Runge-Kutta 某一阶段的时间导数。
     * @tparam RKOrder Runge-Kutta 阶数。
     * @param stage 当前阶段编号。
     * @param S 当前阶段的守恒变量场。
     * @param dSdt 时间导数输出。
     * @param dt 当前子步长。
     * @param fr_as_crse 指向粗层通量寄存器的指针。
     * @param fr_as_fine 指向细层通量寄存器的指针。
     */
    template <int RKOrder>
    void compute_dSdt(int stage, const amrex::MultiFab& S, amrex::MultiFab& dSdt, amrex::Real dt,
        amrex::FluxRegister* fr_as_crse, amrex::FluxRegister* fr_as_fine);

    /**
     * 估计时间步长。
     */
    amrex::Real estTimeStep(amrex::Real dt_old);  // 旧时间步长

    /**
     * 计算初始时间步长。
     */
    amrex::Real initialTimeStep();

    /**
     * 计算初始`dt'。
     */
    virtual void computeInitialDt(int finest_level,              // 最细级别
        int sub_cycle,                                          // 子循环
        amrex::Vector<int>& n_cycle,                            // 循环次数向量
        const amrex::Vector<amrex::IntVect>& ref_ratio,          // 细化比率向量
        amrex::Vector<amrex::Real>& dt_level,                   // 级别时间步长向量
        amrex::Real stop_time) override;                        // 停止时间

    /**
     * 计算新的`dt'。
     */
    virtual void computeNewDt(int finest_level,               // 最细级别
        int sub_cycle,                                       // 子循环
        amrex::Vector<int>& n_cycle,                         // 循环次数向量
        const amrex::Vector<amrex::IntVect>& ref_ratio,       // 细化比率向量
        amrex::Vector<amrex::Real>& dt_min,                  // 最小时间步长向量
        amrex::Vector<amrex::Real>& dt_level,                // 级别时间步长向量
        amrex::Real stop_time,                               // 停止时间
        int post_regrid_flag) override;                      // 重网格后标志

    /**
     * 在timestep()之后工作。
     */
    virtual void post_timestep(int iteration) override;  // 迭代次数

    /**
     * 在regrid()之后工作。
     */
    virtual void post_regrid(int lbase, int new_finest) override;  // 基础级别和新的最细级别

    /**
     * 在restart()之后工作。
     */
    virtual void post_restart() override;

    /**
     * 在init()之后工作。
     */
    virtual void post_init(amrex::Real stop_time) override;  // 停止时间

    /**
     * 重网格的误差估计。
     */
    virtual void errorEst(amrex::TagBoxArray& tb,          // 标记框数组
        int clearval,                                     // 清除值
        int tagval,                                       // 标记值
        amrex::Real time,                                 // 时间
        int n_error_buf = 0, int ngrow = 0) override;     // 误差缓冲区和增长层数

    static int NUM_STATE;  // 状态变量数量
    static int NUM_GROW;   // 虚拟点数量

    // ProbParm类实例
    static ProbParm* h_prob_parm;  // 主机端问题参数
    static ProbParm* d_prob_parm;  // 设备端问题参数

    // 误差标记变量
    static int max_phierr_lev;     // 最大phi误差级别
    static int max_phigrad_lev;    // 最大phi梯度级别
    static amrex::Vector<amrex::Real> phierr;  // phi误差数组
    static amrex::Vector<amrex::Real> phigrad; // phi梯度数组

#if !defined(AMREX_USE_CUDA)
protected:
#endif

    /**
     * @brief 读取 adv 段中的求解控制参数。
     */
    static void read_params();

    /**
     * @brief 在指定方向上执行通量离散与右端项装配。
     * @tparam idim 空间方向编号。
     * @param params 当前方向的通量求解上下文。
     */
    template <int idim>
    void flux_solve_in_dim(const FluxSolveParams& params);  // 通量求解参数
    
    /**
     * 从输入文件读取标记参数的函数。
     * 参见Tagging_params.cpp中的实现。
     */
    static void get_tagging_params();

    // 下面实现的内联函数
    AmrLevelCong& getLevel(int lev);           // 获取指定级别的AmrLevelCong对象
    amrex::FluxRegister& getFluxReg();         // 获取通量寄存器
    amrex::FluxRegister& getFluxReg(int lev);  // 获取指定级别的通量寄存器

    void reflux();  // 通量回流

    void avgDown();  // 向下平均

    void avgDown(int state_indx);  // 向下平均（带状态索引）

    /*
     * 数据。
     */
    amrex::FluxRegister* flux_reg;  // 通量寄存器

    /*
     * 静态数据成员。
     */
    static int verbose;       //  verbose级别
    static amrex::Real cfl;   // CFL数
    static int do_reflux;     // 是否进行通量回流

    // 限制器使用的时间参数
    amrex::GpuArray<amrex::Real, BL_SPACEDIM> tau_xyz;  // 各方向的时间参数
};

/**
 * @brief 为需要用户填充 ghost cell 的边界类型提供占位实现。
 * @details 当前实现不写入任何值，主要用于保持接口完整性；真实用户边界应在该函数基础上扩展。
 */
void nullfill(amrex::Box const& bx,         // 边界框
    amrex::FArrayBox& data,                  // 数据数组
    const int dcomp,                         // 起始组件
    const int numcomp,                       // 组件数量
    amrex::Geometry const& geom,             // 几何信息
    const amrex::Real time,                  // 时间
    const amrex::Vector<amrex::BCRec>& bcr,  // 边界条件记录
    const int bcomp,                         // 边界组件
    const int scomp);                        // 起始组件

/**
 * @brief 计算每个面的 Godunov 速度。
 * @param time 当前时间。
 * @param vx x 方向面速度。
 * @param vy y 方向面速度。
 * @param dx 网格间距。
 * @param prob_lo 物理域下边界。
 */
AMREX_GPU_HOST
AMREX_FORCE_INLINE
void get_face_velocity(const amrex::Real time,                 // 时间
    amrex::FArrayBox& vx,                                     // x方向速度
    amrex::FArrayBox& vy,                                     // y方向速度
    amrex::GpuArray<amrex::Real, BL_SPACEDIM> dx,              // 网格间距
    amrex::GpuArray<amrex::Real, BL_SPACEDIM> prob_lo);        // 问题下界

/*
 * Inlines.
 */

inline AmrLevelCong&
AmrLevelCong::getLevel(int lev)
{
    return *(AmrLevelCong*)&parent->getLevel(lev);
}

inline amrex::FluxRegister&
AmrLevelCong::getFluxReg()
{
    BL_ASSERT(flux_reg);
    return *flux_reg;
}

inline amrex::FluxRegister&
AmrLevelCong::getFluxReg(int lev)
{
    return getLevel(lev).getFluxReg();
}

#endif /*_AmrLevelAdv_H_*/
