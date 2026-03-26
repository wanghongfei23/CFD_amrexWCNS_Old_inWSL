// 方程求解器实现文件（两相流）
// 包含两相流欧拉方程相关的计算和初始化函数

// #include "equation.H"

/**
 * 初始化数据函数，用于设置一维激波管问题的初始条件
 * @param level 网格级别
 * @param time 当前时间
 * @param box 网格边界框
 * @param state 状态数组
 * @param dx 网格间距
 * @param problo 计算域下界
 */
// AMREX_GPU_HOST_DEVICE
// inline void initdata(int level, amrex::Real time,
//     const amrex::Box& box, amrex::FArrayBox& state,
//     const amrex::Real* dx, const amrex::Real problo)
// {
//     using namespace amrex;
//     constexpr Real gamma = 1.4; // 气体绝热指数
//     constexpr Real split_pos = 0.0; // 激波管分割位置
//     AMREX_ASSERT(BL_SPACEDIM == 1); // 确保是一维问题
//     const int* lo = box.loVect(); // 网格下界
//     const int* hi = box.hiVect(); // 网格上界
//     Array4<Real> const& state_array = state.array(); // 状态数组

//     // 遍历所有网格点
// #ifdef AMREX_USE_OMP
// #pragma omp parallel for if (Gpu::notInLaunchRegion())
// #endif
//     for (int i = lo[0]; i <= hi[0]; ++i) {
//         // 计算当前网格中心的物理坐标
//         Real x = problo + (i + 0.5) * dx[0];

//         // 判断当前网格点是否在分割点左侧或右侧
//         Real rho, u, p;
//         if (x <= split_pos) {
//             // 左侧初始条件（高压区）
//             rho = 1.0; // 密度
//             u = 0.0; // 速度
//             p = 1.0; // 压力
//         } else {
//             // 右侧初始条件（低压区）
//             rho = 0.125; // 密度
//             u = 0.0; // 速度
//             p = 0.1; // 压力
//         }

//         // 根据欧拉方程计算守恒变量
//         state_array(i, 0, 0, 0) = rho; // 密度
//         state_array(i, 0, 0, 1) = rho * u; // 动量
//         state_array(i, 0, 0, 2) = p / (gamma - 1) + 0.5 * rho * u * u; // 总能量
//     }
// }