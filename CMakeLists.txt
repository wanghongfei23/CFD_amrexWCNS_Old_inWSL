# 设置 CMake 最小版本要求
cmake_minimum_required(VERSION 3.12)

# 项目名称
project(AMReXConvectionSolver LANGUAGES C CXX)

# 设置构建类型（如果未指定，默认为 Release）
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

# 设置编译器优化选项
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")

# 设置 C++ 标准
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 设置 AMReX 安装路径
# 【临时定位】维数控制
# set(AMReX_DIR /opt/amrex_1d/lib/cmake/AMReX)
set(AMReX_DIR /opt/amrex_2d/lib/cmake/AMReX)
# set(AMReX_DIR /opt/amrex_3d/lib/cmake/AMReX)

# 查找 AMReX 包
find_package(AMReX REQUIRED)

# 检查 AMReX 是否被正确找到
if(NOT AMReX_FOUND)
    message(FATAL_ERROR "AMReX not found! Please check your AMReX installation.")
else()
    message(STATUS "AMReX found at ${AMReX_DIR}")
endif()

# add_compile_options(-Wall -Wextra -Werror)
# add_compile_options(-fsanitize=undefined)

# 新增编译选项
# add_compile_options(-fno-builtin-pow)
# add_compile_options(-fno-builtin-tanh)
# add_compile_options(-fno-builtin-atanh)
# add_compile_options(-fno-tree-vectorize)

# add_compile_options(-fno-reorder-blocks)


# 新增两相流（2phase）相关逻辑
# 使用 file(GLOB) 收集两相流源文件
file(GLOB TWO_PHASE_SOURCES 
    "${CMAKE_SOURCE_DIR}/main.cpp"
    "${CMAKE_SOURCE_DIR}/2phase/*.cpp"
    "${CMAKE_SOURCE_DIR}/amrex_level_Cong/*.cpp"
    "${CMAKE_SOURCE_DIR}/reconstruction_scheme/*.cpp"
    "${CMAKE_SOURCE_DIR}/deriv/*.cpp"
)

# 创建可执行文件（两相流求解器）
add_executable(twoPhaseSolver ${TWO_PHASE_SOURCES})

# 链接 AMReX 库
target_link_libraries(twoPhaseSolver PRIVATE AMReX::amrex)

# 添加头文件路径
target_include_directories(twoPhaseSolver PUBLIC
    ${CMAKE_SOURCE_DIR}/reconstruction_scheme
    ${CMAKE_SOURCE_DIR}/deriv
    ${CMAKE_SOURCE_DIR}/2phase
    ${CMAKE_SOURCE_DIR}/amrex_level_Cong/
)

# 如果 AMReX 使用了其他依赖库（例如 MPI 或 HDF5），需要设置
find_package(MPI REQUIRED)
target_link_libraries(twoPhaseSolver PRIVATE MPI::MPI_CXX)


# # 使用 file(GLOB) 收集欧拉方程相关源文件
# file(GLOB SOURCES 
#     "${CMAKE_SOURCE_DIR}/main.cpp"
#     "${CMAKE_SOURCE_DIR}/amrex_level_Cong/*.cpp"
#     "${CMAKE_SOURCE_DIR}/reconstruction_scheme/*.cpp"
#     "${CMAKE_SOURCE_DIR}/deriv/*.cpp"
# )

# # 根据 AMReX_SPACEDIM 的值来选择是否添加特定目录的源文件（欧拉方程相关）
# if(AMReX_SPACEDIM EQUAL 1)
#     message(STATUS "Including sources from Euler_XD folder for 1D space")
#     file(GLOB SOURCES_1D "${CMAKE_SOURCE_DIR}/Euler1D/*.cpp")
#     list(APPEND SOURCES ${SOURCES_1D})
# elseif(AMReX_SPACEDIM EQUAL 2)
#     message(STATUS "Including sources from Euler_XD folder for 2D space")
#     file(GLOB SOURCES_2D "${CMAKE_SOURCE_DIR}/Euler2D/*.cpp")
#     list(APPEND SOURCES ${SOURCES_2D})
# elseif(AMReX_SPACEDIM EQUAL 3)
#     message(STATUS "Including sources from Euler_XD folder for 3D space")
#     file(GLOB SOURCES_3D "${CMAKE_SOURCE_DIR}/Euler3D/*.cpp")
#     list(APPEND SOURCES ${SOURCES_3D})
# endif()

# # 创建第一个可执行文件（欧拉方程求解器）
# add_executable(congSolver ${SOURCES})

# # 链接 AMReX 库
# target_link_libraries(congSolver PRIVATE AMReX::amrex)

# # 添加头文件路径
# target_include_directories(congSolver PUBLIC
#     ${CMAKE_SOURCE_DIR}/reconstruction_scheme
#     ${CMAKE_SOURCE_DIR}/deriv
#     ${CMAKE_SOURCE_DIR}/amrex_level_Cong/
# )

# if(AMReX_SPACEDIM EQUAL 1)
#     target_include_directories(congSolver PUBLIC ${CMAKE_SOURCE_DIR}/Euler1D)
# elseif(AMReX_SPACEDIM EQUAL 2)
#     target_include_directories(congSolver PUBLIC ${CMAKE_SOURCE_DIR}/Euler2D)
# elseif(AMReX_SPACEDIM EQUAL 3)
#     target_include_directories(congSolver PUBLIC ${CMAKE_SOURCE_DIR}/Euler3D)
# endif()

# # 如果 AMReX 使用了其他依赖库（例如 MPI 或 HDF5），需要设置
# find_package(MPI REQUIRED)
# target_link_libraries(congSolver PRIVATE MPI::MPI_CXX)