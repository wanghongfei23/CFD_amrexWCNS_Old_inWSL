/**
 * @file
 * @brief AMReX 层级构建器 LevelBldAdv 的实现。
 * @ingroup amr_core
 */

#include "AmrLevelCong.H"
#include <AMReX_LevelBld.H>

using namespace amrex;

// LevelBldAdv类，继承自LevelBld
// 用于构建和管理AMR级别
class LevelBldAdv
    : public LevelBld {
    // 变量设置函数
    virtual void variableSetUp() override;
    // 变量清理函数
    virtual void variableCleanUp() override;
    // 默认构造函数
    virtual AmrLevel* operator()() override;
    // 带参数的构造函数
    virtual AmrLevel* operator()(Amr& papa,
        int lev,
        const Geometry& level_geom,
        const BoxArray& ba,
        const DistributionMapping& dm,
        Real time) override;
};

// 全局LevelBldAdv实例
LevelBldAdv Adv_bld;

// 获取级别构建器函数
// 输出：LevelBld指针
LevelBld*
getLevelBld()
{
    return &Adv_bld;
}

// 变量设置函数实现
void LevelBldAdv::variableSetUp()
{
    AmrLevelCong::variableSetUp();
}

// 变量清理函数实现
void LevelBldAdv::variableCleanUp()
{
    AmrLevelCong::variableCleanUp();
}

// 默认构造函数实现
// 输出：AmrLevel指针
AmrLevel*
LevelBldAdv::operator()()
{
    return new AmrLevelCong;
}

// 带参数的构造函数实现
// 输入：Amr对象papa，级别lev，几何信息level_geom，网格数组ba，分布映射dm，时间time
// 输出：AmrLevel指针
AmrLevel*
LevelBldAdv::operator()(Amr& papa,
    int lev,
    const Geometry& level_geom,
    const BoxArray& ba,
    const DistributionMapping& dm,
    Real time)
{
    return new AmrLevelCong(papa, lev, level_geom, ba, dm, time);
}
