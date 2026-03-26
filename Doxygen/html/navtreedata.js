/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "CFD_inWSL_MF_initial", "index.html", [
    [ "代码结构详解", "d2/d86/md_project__analysis_2code__structure.html", [
      [ "根目录文件", "d2/d86/md_project__analysis_2code__structure.html#autotoc_md1", null ],
      [ "2phase 目录", "d2/d86/md_project__analysis_2code__structure.html#autotoc_md2", null ],
      [ "amrex_level_Cong 目录", "d2/d86/md_project__analysis_2code__structure.html#autotoc_md3", null ],
      [ "deriv 目录", "d2/d86/md_project__analysis_2code__structure.html#autotoc_md4", null ],
      [ "reconstruction_scheme 目录", "d2/d86/md_project__analysis_2code__structure.html#autotoc_md5", null ],
      [ "data 目录", "d2/d86/md_project__analysis_2code__structure.html#autotoc_md6", null ],
      [ "核心类和函数", "d2/d86/md_project__analysis_2code__structure.html#autotoc_md7", [
        [ "AmrLevelCong 类", "d2/d86/md_project__analysis_2code__structure.html#autotoc_md8", null ],
        [ "FirstOrderDerivativeTerm 模板", "d2/d86/md_project__analysis_2code__structure.html#autotoc_md9", null ],
        [ "HLLC 求解器", "d2/d86/md_project__analysis_2code__structure.html#autotoc_md10", null ],
        [ "状态方程 (EOS)", "d2/d86/md_project__analysis_2code__structure.html#autotoc_md11", null ]
      ] ],
      [ "代码执行流程", "d2/d86/md_project__analysis_2code__structure.html#autotoc_md12", null ],
      [ "依赖关系", "d2/d86/md_project__analysis_2code__structure.html#autotoc_md13", null ],
      [ "编译选项", "d2/d86/md_project__analysis_2code__structure.html#autotoc_md14", null ],
      [ "代码特点", "d2/d86/md_project__analysis_2code__structure.html#autotoc_md15", null ]
    ] ],
    [ "项目分析目录结构", "db/de5/md_project__analysis_2directory__structure.html", [
      [ "目录组织", "db/de5/md_project__analysis_2directory__structure.html#autotoc_md17", null ],
      [ "文件说明", "db/de5/md_project__analysis_2directory__structure.html#autotoc_md18", [
        [ "README.md", "db/de5/md_project__analysis_2directory__structure.html#autotoc_md19", null ],
        [ "code_structure.md", "db/de5/md_project__analysis_2directory__structure.html#autotoc_md20", null ],
        [ "visualization_guide.md", "db/de5/md_project__analysis_2directory__structure.html#autotoc_md21", null ],
        [ "quick_start.md", "db/de5/md_project__analysis_2directory__structure.html#autotoc_md22", null ]
      ] ],
      [ "如何使用这些文档", "db/de5/md_project__analysis_2directory__structure.html#autotoc_md23", null ],
      [ "相关资源", "db/de5/md_project__analysis_2directory__structure.html#autotoc_md24", null ],
      [ "联系信息", "db/de5/md_project__analysis_2directory__structure.html#autotoc_md25", null ]
    ] ],
    [ "快速上手指南", "d2/da0/md_project__analysis_2quick__start.html", [
      [ "环境准备", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md27", [
        [ "1. 系统要求", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md28", null ],
        [ "2. AMReX 安装", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md29", null ],
        [ "3. 项目配置", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md30", null ]
      ] ],
      [ "编译项目", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md31", null ],
      [ "运行示例", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md32", [
        [ "1. 2D 双相流算例", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md33", null ],
        [ "2. 1D 激波管问题", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md34", null ],
        [ "3. 高速流动算例", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md35", null ]
      ] ],
      [ "配置文件说明", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md36", [
        [ "inputs2phase2d 主要参数", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md37", null ]
      ] ],
      [ "结果查看", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md38", [
        [ "1. 输出文件", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md39", null ],
        [ "2. 使用 ParaView 查看", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md40", null ],
        [ "3. 数据导出", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md41", null ]
      ] ],
      [ "常见问题", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md42", [
        [ "1. 编译错误", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md43", null ],
        [ "2. 运行错误", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md44", null ],
        [ "3. 结果异常", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md45", null ]
      ] ],
      [ "调优建议", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md46", [
        [ "1. 性能优化", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md47", null ],
        [ "2. 数值参数调整", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md48", null ],
        [ "3. 内存管理", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md49", null ]
      ] ],
      [ "高级用法", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md50", [
        [ "1. 自定义初始条件", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md51", null ],
        [ "2. 自定义物理模型", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md52", null ],
        [ "3. 自定义数值格式", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md53", null ]
      ] ],
      [ "故障排除", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md54", [
        [ "1. 检查 AMReX 版本", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md55", null ],
        [ "2. 检查编译配置", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md56", null ],
        [ "3. 运行调试版本", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md57", null ],
        [ "4. 查看日志文件", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md58", null ]
      ] ],
      [ "下一步", "d2/da0/md_project__analysis_2quick__start.html#autotoc_md59", null ]
    ] ],
    [ "可视化指南", "de/d4b/md_project__analysis_2visualization__guide.html", [
      [ "推荐的可视化工具", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md83", null ],
      [ "建议的可视化图表", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md84", [
        [ "1. 密度分布图", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md85", null ],
        [ "2. 体积分数分布", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md86", null ],
        [ "3. 压力场分布", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md87", null ],
        [ "4. 速度矢量场", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md88", null ],
        [ "5. 网格加密图", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md89", null ],
        [ "6. 时间序列分析", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md90", null ]
      ] ],
      [ "数据导出方法", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md91", [
        [ "1. AMReX 内置输出", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md92", null ],
        [ "2. 自定义导出", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md93", null ],
        [ "3. 1D数据提取", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md94", null ]
      ] ],
      [ "可视化技巧", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md95", [
        [ "1. 多变量对比", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md96", null ],
        [ "2. 数值精度检查", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md97", null ],
        [ "3. 网格收敛性分析", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md98", null ]
      ] ],
      [ "常见问题及解决方案", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md99", [
        [ "1. 数据量过大", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md100", null ],
        [ "2. 可视化软件崩溃", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md101", null ],
        [ "3. 相界面模糊", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md102", null ],
        [ "4. 激波捕捉不准确", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md103", null ]
      ] ],
      [ "示例可视化命令", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md104", [
        [ "ParaView", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md105", null ],
        [ "VisIt", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md106", null ],
        [ "Python/Matplotlib", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md107", null ]
      ] ],
      [ "高级可视化", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md108", [
        [ "1. 三维可视化", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md109", null ],
        [ "2. 并行可视化", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md110", null ],
        [ "3. 交互式分析", "de/d4b/md_project__analysis_2visualization__guide.html#autotoc_md111", null ]
      ] ]
    ] ],
    [ "CFD求解器项目介绍", "d3/dcc/md__r_e_a_d_m_e.html", [
      [ "1. 项目概述", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md113", [
        [ "1.1 项目特点", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md114", null ]
      ] ],
      [ "2. 项目结构与文件说明", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md115", [
        [ "2.1 主要文件", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md116", null ],
        [ "2.2 目录结构", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md117", null ]
      ] ],
      [ "3. 核心功能与实现", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md118", [
        [ "3.1 网格系统", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md119", null ],
        [ "3.2 求解器", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md120", null ],
        [ "3.3 工具函数", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md121", null ]
      ] ],
      [ "4. CFD数值模拟具体逻辑分析", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md122", [
        [ "4.1 整体求解流程", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md123", null ],
        [ "4.2 各部分分块与数据传输", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md124", [
          [ "4.2.1 网格系统模块", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md125", null ],
          [ "4.2.2 求解器模块", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md126", null ],
          [ "4.2.3 工具函数模块", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md127", null ]
        ] ],
        [ "4.3 关键算法实现细节", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md128", [
          [ "4.3.1 有限体积法离散化", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md129", null ],
          [ "4.3.2 SIMPLE算法实现", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md130", null ],
          [ "4.3.3 边界条件处理", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md131", null ]
        ] ],
        [ "4.4 数据流动图", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md132", null ],
        [ "4.5 计算效率优化", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md133", null ]
      ] ],
      [ "5. C++项目构建知识", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md134", [
        [ "5.1 Makefile基础", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md135", null ],
        [ "5.2 编译与链接过程", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md136", null ],
        [ "5.3 常见编译选项", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md137", null ],
        [ "5.4 项目组织最佳实践", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md138", null ]
      ] ],
      [ "6. 构建与运行", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md139", [
        [ "6.1 构建项目", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md140", null ],
        [ "6.2 运行程序", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md141", null ],
        [ "6.3 清理构建文件", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md142", null ]
      ] ],
      [ "7. 扩展与改进", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md143", [
        [ "7.1 可能的扩展方向", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md144", null ],
        [ "7.2 代码优化建议", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md145", null ]
      ] ],
      [ "8. 总结", "d3/dcc/md__r_e_a_d_m_e.html#autotoc_md146", null ]
    ] ],
    [ "命名空间", "namespaces.html", [
      [ "命名空间列表", "namespaces.html", "namespaces_dup" ],
      [ "命名空间成员", "namespacemembers.html", [
        [ "全部", "namespacemembers.html", null ],
        [ "函数", "namespacemembers_func.html", null ]
      ] ]
    ] ],
    [ "类", "annotated.html", [
      [ "类列表", "annotated.html", "annotated_dup" ],
      [ "类索引", "classes.html", null ],
      [ "类继承关系", "hierarchy.html", "hierarchy" ],
      [ "类成员", "functions.html", [
        [ "全部", "functions.html", null ],
        [ "函数", "functions_func.html", null ],
        [ "变量", "functions_vars.html", null ]
      ] ]
    ] ],
    [ "文件", "files.html", [
      [ "文件列表", "files.html", "files_dup" ],
      [ "文件成员", "globals.html", [
        [ "全部", "globals.html", null ],
        [ "函数", "globals_func.html", null ],
        [ "变量", "globals_vars.html", null ],
        [ "类型定义", "globals_type.html", null ],
        [ "枚举", "globals_enum.html", null ],
        [ "枚举值", "globals_eval.html", null ],
        [ "宏定义", "globals_defs.html", null ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"annotated.html",
"d5/da8/fifth__order__thinc__recon_8_h.html#a6fc4f383a933ec63ca8493eb361d9738"
];

var SYNCONMSG = '点击 关闭 面板同步';
var SYNCOFFMSG = '点击 开启 面板同步';
var LISTOFALLMEMBERS = '所有成员列表';