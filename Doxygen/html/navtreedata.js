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
    [ "代码结构详解", "db/dcd/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_dae505ce405bbfdac7c1a7ca2dc6dd68.html", [
      [ "根目录文件", "db/dcd/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_dae505ce405bbfdac7c1a7ca2dc6dd68.html#autotoc_md1", null ],
      [ "2phase 目录", "db/dcd/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_dae505ce405bbfdac7c1a7ca2dc6dd68.html#autotoc_md2", null ],
      [ "amrex_level_Cong 目录", "db/dcd/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_dae505ce405bbfdac7c1a7ca2dc6dd68.html#autotoc_md3", null ],
      [ "deriv 目录", "db/dcd/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_dae505ce405bbfdac7c1a7ca2dc6dd68.html#autotoc_md4", null ],
      [ "reconstruction_scheme 目录", "db/dcd/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_dae505ce405bbfdac7c1a7ca2dc6dd68.html#autotoc_md5", null ],
      [ "data 目录", "db/dcd/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_dae505ce405bbfdac7c1a7ca2dc6dd68.html#autotoc_md6", null ],
      [ "核心类和函数", "db/dcd/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_dae505ce405bbfdac7c1a7ca2dc6dd68.html#autotoc_md7", [
        [ "AmrLevelCong 类", "db/dcd/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_dae505ce405bbfdac7c1a7ca2dc6dd68.html#autotoc_md8", null ],
        [ "FirstOrderDerivativeTerm 模板", "db/dcd/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_dae505ce405bbfdac7c1a7ca2dc6dd68.html#autotoc_md9", null ],
        [ "HLLC 求解器", "db/dcd/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_dae505ce405bbfdac7c1a7ca2dc6dd68.html#autotoc_md10", null ],
        [ "状态方程 (EOS)", "db/dcd/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_dae505ce405bbfdac7c1a7ca2dc6dd68.html#autotoc_md11", null ]
      ] ],
      [ "代码执行流程", "db/dcd/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_dae505ce405bbfdac7c1a7ca2dc6dd68.html#autotoc_md12", null ],
      [ "依赖关系", "db/dcd/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_dae505ce405bbfdac7c1a7ca2dc6dd68.html#autotoc_md13", null ],
      [ "编译选项", "db/dcd/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_dae505ce405bbfdac7c1a7ca2dc6dd68.html#autotoc_md14", null ],
      [ "代码特点", "db/dcd/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_dae505ce405bbfdac7c1a7ca2dc6dd68.html#autotoc_md15", null ]
    ] ],
    [ "项目分析目录结构", "d5/dba/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_085a407c629d12b8a812c5fe6461feef.html", [
      [ "目录组织", "d5/dba/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_085a407c629d12b8a812c5fe6461feef.html#autotoc_md17", null ],
      [ "文件说明", "d5/dba/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_085a407c629d12b8a812c5fe6461feef.html#autotoc_md18", [
        [ "README.md", "d5/dba/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_085a407c629d12b8a812c5fe6461feef.html#autotoc_md19", null ],
        [ "code_structure.md", "d5/dba/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_085a407c629d12b8a812c5fe6461feef.html#autotoc_md20", null ],
        [ "visualization_guide.md", "d5/dba/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_085a407c629d12b8a812c5fe6461feef.html#autotoc_md21", null ],
        [ "quick_start.md", "d5/dba/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_085a407c629d12b8a812c5fe6461feef.html#autotoc_md22", null ]
      ] ],
      [ "如何使用这些文档", "d5/dba/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_085a407c629d12b8a812c5fe6461feef.html#autotoc_md23", null ],
      [ "相关资源", "d5/dba/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_085a407c629d12b8a812c5fe6461feef.html#autotoc_md24", null ],
      [ "联系信息", "d5/dba/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_085a407c629d12b8a812c5fe6461feef.html#autotoc_md25", null ]
    ] ],
    [ "快速上手指南", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html", [
      [ "环境准备", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md27", [
        [ "1. 系统要求", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md28", null ],
        [ "2. AMReX 安装", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md29", null ],
        [ "3. 项目配置", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md30", null ]
      ] ],
      [ "编译项目", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md31", null ],
      [ "运行示例", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md32", [
        [ "1. 2D 双相流算例", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md33", null ],
        [ "2. 1D 激波管问题", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md34", null ],
        [ "3. 高速流动算例", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md35", null ]
      ] ],
      [ "配置文件说明", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md36", [
        [ "inputs2phase2d 主要参数", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md37", null ]
      ] ],
      [ "结果查看", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md38", [
        [ "1. 输出文件", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md39", null ],
        [ "2. 使用 ParaView 查看", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md40", null ],
        [ "3. 数据导出", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md41", null ]
      ] ],
      [ "常见问题", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md42", [
        [ "1. 编译错误", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md43", null ],
        [ "2. 运行错误", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md44", null ],
        [ "3. 结果异常", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md45", null ]
      ] ],
      [ "调优建议", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md46", [
        [ "1. 性能优化", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md47", null ],
        [ "2. 数值参数调整", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md48", null ],
        [ "3. 内存管理", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md49", null ]
      ] ],
      [ "高级用法", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md50", [
        [ "1. 自定义初始条件", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md51", null ],
        [ "2. 自定义物理模型", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md52", null ],
        [ "3. 自定义数值格式", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md53", null ]
      ] ],
      [ "故障排除", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md54", [
        [ "1. 检查 AMReX 版本", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md55", null ],
        [ "2. 检查编译配置", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md56", null ],
        [ "3. 运行调试版本", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md57", null ],
        [ "4. 查看日志文件", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md58", null ]
      ] ],
      [ "下一步", "d4/d7b/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_05eece26b5b0e259b3d61386a29a28b7.html#autotoc_md59", null ]
    ] ],
    [ "可视化指南", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html", [
      [ "推荐的可视化工具", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md83", null ],
      [ "建议的可视化图表", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md84", [
        [ "1. 密度分布图", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md85", null ],
        [ "2. 体积分数分布", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md86", null ],
        [ "3. 压力场分布", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md87", null ],
        [ "4. 速度矢量场", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md88", null ],
        [ "5. 网格加密图", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md89", null ],
        [ "6. 时间序列分析", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md90", null ]
      ] ],
      [ "数据导出方法", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md91", [
        [ "1. AMReX 内置输出", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md92", null ],
        [ "2. 自定义导出", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md93", null ],
        [ "3. 1D数据提取", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md94", null ]
      ] ],
      [ "可视化技巧", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md95", [
        [ "1. 多变量对比", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md96", null ],
        [ "2. 数值精度检查", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md97", null ],
        [ "3. 网格收敛性分析", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md98", null ]
      ] ],
      [ "常见问题及解决方案", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md99", [
        [ "1. 数据量过大", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md100", null ],
        [ "2. 可视化软件崩溃", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md101", null ],
        [ "3. 相界面模糊", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md102", null ],
        [ "4. 激波捕捉不准确", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md103", null ]
      ] ],
      [ "示例可视化命令", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md104", [
        [ "ParaView", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md105", null ],
        [ "VisIt", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md106", null ],
        [ "Python/Matplotlib", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md107", null ]
      ] ],
      [ "高级可视化", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md108", [
        [ "1. 三维可视化", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md109", null ],
        [ "2. 并行可视化", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md110", null ],
        [ "3. 交互式分析", "d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md111", null ]
      ] ]
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
"d6/d52/md__2_2wsl_8localhost_2arch__zyc_2home_2archwanghongfei_2_documents_2_git_hub_w_s_l_2_c_f_d__in_bc0712e76f60344e7fb65ff9331ce10b.html#autotoc_md98"
];

var SYNCONMSG = '点击 关闭 面板同步';
var SYNCOFFMSG = '点击 开启 面板同步';
var LISTOFALLMEMBERS = '所有成员列表';