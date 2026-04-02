/**
 * @file
 * @brief 两相状态方程与变量互转测试程序。
 * @ingroup analysis_docs
 */

#include "equation.H"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace amrex;

/**
 * @brief 打印状态数组。
 * @param name 标签名。
 * @param var 状态数组。
 * @param out 输出流。
 */
void printVarArray(const std::string& name, const VarArray& var, std::ostream& out)
{
    out << name << ": [ ";
    for (int i = 0; i < NVar; ++i) {
        out << std::setw(12) << var[i];
    }
    out << " ]" << std::endl;
}

/**
 * @brief 打印三维向量。
 * @param name 标签名。
 * @param arr 三维向量。
 * @param out 输出流。
 */
void printGpuArray3(const std::string& name, const amrex::GpuArray<amrex::Real, 3>& arr, std::ostream& out)
{
    out << name << ": [ " << arr[0] << ", " << arr[1] << ", " << arr[2] << " ]" << std::endl;
}

/**
 * @brief 运行单个状态方程测试用例。
 * @param testName 用例名称。
 * @param cons 守恒变量输入。
 * @param pp 物性参数。
 * @param outFile 输出流。
 * @param dimension 维度编号。
 */
void runTestCase(
    const std::string& testName,
    const VarArray& cons,
    const ProbParm& pp,
    std::ostream& outFile,
    int dimension = BL_SPACEDIM) // 问题的维度（1D, 2D, 3D）
{
    outFile << "=== Test Case: " << testName << " ===\n";
    std::cout << "=== Test Case: " << testName << " ===\n";

    // 打印守恒变量
    printVarArray("Conserved Variables", cons, outFile);

    // 转换为原始变量
    VarArray prim = cons_to_prim(cons, pp);
    printVarArray("Primitive Variables", prim, outFile);

    // 验证可逆性
    VarArray cons_reconstructed = prim_to_cons(prim, pp);
    printVarArray("Reconstructed Conserved Variables", cons_reconstructed, outFile);

    // 检查重构守恒变量是否与原始守恒变量一致（容差 1e-5）
    bool isConsEqual = true;
    for (int i = 0; i < cons.size(); ++i) {
        if (std::abs(cons[i] - cons_reconstructed[i]) > 1e-5) {
            isConsEqual = false;
            break;
        }
    }

    if (isConsEqual) {
        outFile << "Conservation check passed: cons_reconstructed == cons\n";
    } else {
        outFile << "Conservation check failed: cons_reconstructed != cons\n";
    }

    // 状态方程计算
    Real pressure = cons_to_eos<Pressure>(cons, pp);
    Real total_energy = prims_to_eos<TotalEnergy>(prim, pp);
    outFile << "Pressure from EOS: " << pressure << "\n";
    outFile << "Total Energy from EOS: " << total_energy << "\n";

    // 根据维度测试不同的法向量 norm
    std::vector<amrex::GpuArray<amrex::Real, 3>> norms;
    if (dimension == 1) {
        norms.push_back({ 1.0, 0.0, 0.0 }); // 1D: 只测试 norm = (1, 0, 0)
    } else if (dimension == 2) {
        norms.push_back({ 1.0, 0.0, 0.0 }); // 2D: 测试 norm = (1, 0, 0)
        norms.push_back({ 0.0, 1.0, 0.0 }); //      和 norm = (0, 1, 0)
    } else if (dimension == 3) {
        norms.push_back({ 1.0, 0.0, 0.0 }); // 3D: 测试 norm = (1, 0, 0)
        norms.push_back({ 0.0, 1.0, 0.0 }); //      和 norm = (0, 1, 0)
        norms.push_back({ 0.0, 0.0, 1.0 }); //      和 norm = (0, 0, 1)
    }

    for (const auto& norm : norms) {
        outFile << "Testing flux with norm = ("
                << norm[0] << ", " << norm[1] << ", " << norm[2] << ")\n";

        // 计算通量
        VarArray flux = prim_to_flux(prim, pp, norm);
        printVarArray("Flux", flux, outFile);
    }

    outFile << "\n";
}

/**
 * @brief 执行状态方程与通量测试程序。
 * @return 进程退出码。
 */
int main()
{
    // 打开输出文件
    std::ofstream outFile("test_results.txt");
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open output file." << std::endl;
        return 1;
    }

    // 定义问题参数
    ProbParm pp;

    // 测试用例 1: 纯相1区域
    {
        VarArray cons = {
            1.0, 0.0,
            AMREX_D_DECL(0.5, 0.5, 0.5), // XMom, YMom, ZMom (根据维度自动展开)
            10.0, 1.0
        }; // AlphaRho1, AlphaRho2, XMom, YMom, ZMom, Energy, Alpha1
        runTestCase("Pure Phase 1", cons, pp, outFile);
    }

    // 测试用例 2: 纯相2区域
    {
        VarArray cons = {
            0.0, 1.0,
            AMREX_D_DECL(0.5, 0.5, 0.5), // XMom, YMom, ZMom (根据维度自动展开)
            10.0, 0.0
        }; // AlphaRho1, AlphaRho2, XMom, YMom, ZMom, Energy, Alpha1
        runTestCase("Pure Phase 2", cons, pp, outFile);
    }

    // 测试用例 3: 混合区
    {
        VarArray cons = {
            0.5, 0.5,
            AMREX_D_DECL(0.5, 0.5, 0.5), // XMom, YMom, ZMom (根据维度自动展开)
            10.0, 0.5
        }; // AlphaRho1, AlphaRho2, XMom, YMom, ZMom, Energy, Alpha1
        runTestCase("Mixed Zone", cons, pp, outFile);
    }

    // 测试用例 4: 区域 0.25 <= x < 0.75
    {
        VarArray prim = {
            1000.0, 1.0e-8,
            AMREX_D_DECL(100.0, 0.0, 0.0), // XMom, YMom, ZMom (根据维度自动展开)
            101325.0, 1.0 - 1.0e-8
        }; // AlphaRho1, AlphaRho2, XMom, YMom, ZMom, Energy, Alpha1
        VarArray cons = prim_to_cons(prim, pp);
        runTestCase("Region 0.25 <= x < 0.75", cons, pp, outFile);
    }

    // 测试用例 5: 区域 otherwise
    {
        VarArray prim = {
            1.0e-8, 1.204,
            AMREX_D_DECL(100.0, 0.0, 0.0), // XMom, YMom, ZMom (根据维度自动展开)
            101325.0, 1.0e-8
        }; // AlphaRho1, AlphaRho2, XMom, YMom, ZMom, Energy, Alpha1
        VarArray cons = prim_to_cons(prim, pp);
        runTestCase("Region otherwise", cons, pp, outFile);
    }

    // 关闭文件
    outFile.close();
    std::cout << "Test results written to test_results.txt\n";
    return 0;
}