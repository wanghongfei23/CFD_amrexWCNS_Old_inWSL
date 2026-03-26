// 方程测试文件
// 包含方程求解器测试相关的函数

#include "equation.H"
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace amrex;

// 辅助函数：打印 VarArray
// 输入：名称name，变量数组var，输出流out
void printVarArray(const std::string& name, const VarArray& var, std::ostream& out)
{
    out << name << ": [ ";
    for (int i = 0; i < NVar; ++i) {
        out << std::setw(12) << var[i];
    }
    out << " ]" << std::endl;
}

// 辅助函数：打印 GpuArray<Real, 3>
// 输入：名称name，GPU数组arr，输出流out
void printGpuArray3(const std::string& name, const amrex::GpuArray<amrex::Real, 3>& arr, std::ostream& out)
{
    out << name << ": [ " << arr[0] << ", " << arr[1] << ", " << arr[2] << " ]" << std::endl;
}

// 辅助函数：运行单个测试用例
// 输入：测试名称testName，守恒变量cons，问题参数pp，输出文件outFile，维度dimension
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
    Real pressure = eos_cons<Pressure>(cons, pp);
    Real total_energy = eos_prim<TotalEnergy>(prim, pp);
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