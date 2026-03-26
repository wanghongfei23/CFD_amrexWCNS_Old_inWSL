// 特征系统测试文件
// 包含特征矩阵测试相关的函数

#include "eigen_system.H"
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>

// 辅助函数：打印矩阵（用于调试）
// 输入：矩阵mat，行数rows，列数cols
void printMatrix(const VarMatrix& mat, int rows, int cols)
{
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << mat[i * cols + j] << " ";
        }
        std::cout << "\n";
    }
}

// 辅助函数：检查矩阵是否为单位矩阵
// 输入：矩阵mat，大小size，容差tol
// 输出：是否为单位矩阵的布尔值
bool isIdentityMatrix(const VarMatrix& mat, int size, amrex::Real tol = 1e-6)
{
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            amrex::Real expected = (i == j) ? 1.0 : 0.0;
            if (std::abs(mat[i * size + j] - expected) > tol) {
                return false;
            }
        }
    }
    return true;
}

// 打印元素公式函数
// 输入：左特征矩阵leftEig，右特征矩阵rightEig，行row，列col，大小size
void printElementFormula(const VarMatrix& leftEig, const VarMatrix& rightEig, int row, int col, int size)
{
    std::cout << "Formula for product[" << row << "][" << col << "]:\n";
    std::cout << "product[" << row << "][" << col << "] = ";

    // 用于存储最终计算结果
    amrex::Real result = 0.0;

    for (int k = 0; k < size; ++k) {
        // 获取左特征矩阵和右特征矩阵的元素值
        amrex::Real leftValue = leftEig[row * size + k];
        amrex::Real rightValue = rightEig[k * size + col];

        // 打印公式和具体数值
        std::cout << "(" << leftValue << " * " << rightValue << ")";
        if (k < size - 1) {
            std::cout << " + ";
        }

        // 累加计算结果
        result += leftValue * rightValue;
    }

    // 打印最终结果
    std::cout << "\nResult: product[" << row << "][" << col << "] = " << result << "\n\n";
}

// 测试主函数
// 功能：测试特征系统的正确性
int main()
{
    // 初始化输入参数
    VarArray primL = {
        1.0, // AlphaRho1 (第一相密度)
        2.0, // AlphaRho2 (第二相密度)
        AMREX_D_DECL(3.0, 4.0, 5.0), // 动量分量 (XMom, YMom, ZMom)
        6.0, // 能量 (Energy)
        0.5 // Alpha1 (第一相体积分数)
    };

    VarArray primR = {
        1.5, // AlphaRho1 (第一相密度)
        2.5, // AlphaRho2 (第二相密度)
        AMREX_D_DECL(3.5, 4.5, 5.5), // 动量分量 (XMom, YMom, ZMom)
        6.5, // 能量 (Energy)
        0.6 // Alpha1 (第一相体积分数)
    };

    amrex::GpuArray<amrex::Real, 3> norm = {
        AMREX_D_DECL(1.0, 0.0, 0.0) // 法向量 (x, y, z)
        // AMREX_D_DECL(0.0, 1.0, 0.0) // 法向量 (x, y, z)
        // AMREX_D_DECL(0.0, 0.0, 1.0) // 法向量 (x, y, z)
    };

    ProbParm pp; // 假设 ProbParm 是一个空结构体

    // 创建 EigenSystem 对象
    EigenSystem eigenSystem(primL, primR, norm, pp);

    // 仅在调试模式下通过接口获取特征矩阵
    const VarMatrix& leftEig = eigenSystem.getLeftEig();
    const VarMatrix& rightEig = eigenSystem.getRightEig();

    // 矩阵维度
    constexpr int size = NVar; // 使用 NVar 定义矩阵大小

    printElementFormula(leftEig, rightEig, 1, 1, size); // 打印特定元素的公式

    // 计算左右特征矩阵的乘积
    VarMatrix product = { 0.0 }; // 初始化结果矩阵
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            for (int k = 0; k < size; ++k) {
                product[i * size + j] += leftEig[i * size + k] * rightEig[k * size + j];
            }
        }
    }

    // 打印乘积矩阵（可选，用于调试）
    std::cout << "Product of leftEig and rightEig:\n";
    printMatrix(product, size, size);

    // 检查乘积矩阵是否为单位矩阵
    bool isIdentity = isIdentityMatrix(product, size);
    if (isIdentity) {
        std::cout << "Test passed: The product of leftEig and rightEig is an identity matrix.\n";
    } else {
        std::cout << "Test failed: The product of leftEig and rightEig is NOT an identity matrix.\n";
    }

    return isIdentity ? 0 : 1;
}