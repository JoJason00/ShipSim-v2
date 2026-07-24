#include "SparsityAnalyzer.h"
#include <iostream>
#include <iomanip>
#include <stdexcept>

// 构造函数
MatrixSparsityAnalyzer::MatrixSparsityAnalyzer(double eps) : eps_(eps) {
    if (eps < 0) {
        throw std::invalid_argument("Epsilon must be non-negative!");
    }
}

// 统计单个矩阵的稀疏度
double MatrixSparsityAnalyzer::calculateSingleMatrixSparsity(const Eigen::MatrixXd& mat,
    bool print_detail) const {
    if (mat.size() == 0) {
        if (print_detail) {
            std::cout << "[SparsityAnalyzer] 矩阵为空，稀疏度为0" << std::endl;
        }
        return 0.0;
    }

    // 统计非零元素数量
    int non_zero_count = fastNonZeroCount(mat);
    double sparsity = static_cast<double>(non_zero_count) / mat.size();

    // 打印详细信息（可选）
    if (print_detail) {
        std::cout << "\n=== 单个矩阵稀疏度统计 ===" << std::endl;
        std::cout << "矩阵维度: " << mat.rows() << "x" << mat.cols() << std::endl;
        std::cout << "总元素数: " << mat.size() << std::endl;
        std::cout << "非零元素数: " << non_zero_count << std::endl;
        std::cout << "非零元素占比: " << std::fixed << std::setprecision(4) << sparsity * 100 << "%" << std::endl;
        std::cout << "零元素占比: " << std::fixed << std::setprecision(4) << (1 - sparsity) * 100 << "%" << std::endl;
    }

    return sparsity;
}

// 统计矩阵容器的平均稀疏度
double MatrixSparsityAnalyzer::calculateMatrixContainerSparsity(const std::vector<Eigen::MatrixXd>& mat_container,
    const std::string& container_name,
    int max_count) const {
    if (mat_container.empty()) {
        std::cout << "[SparsityAnalyzer] " << container_name << " 容器为空！" << std::endl;
        return 0.0;
    }

    double total_sparsity = 0.0;
    int valid_matrix_count = 0;
    int actual_max_count = (max_count <= 0) ? mat_container.size() : max_count;

    std::cout << "\n=== 开始统计容器 [" << container_name << "] 的稀疏度 ===" << std::endl;
    std::cout << "容器总矩阵数: " << mat_container.size() << std::endl;
    std::cout << "计划统计矩阵数: " << actual_max_count << std::endl;

    // 遍历容器中的矩阵
    for (int i = 0; i < mat_container.size() && i < actual_max_count; ++i) {
        const auto& mat = mat_container[i];
        if (mat.size() == 0) {
            std::cout << "[SparsityAnalyzer] 第" << i << "个矩阵为空，跳过" << std::endl;
            continue;
        }

        double sparsity = calculateSingleMatrixSparsity(mat, true); // 打印单个矩阵详情
        total_sparsity += sparsity;
        valid_matrix_count++;
    }

    // 计算平均稀疏度
    double avg_sparsity = valid_matrix_count > 0 ? (total_sparsity / valid_matrix_count) : 0.0;

    // 汇总输出
    std::cout << "\n=== 容器 [" << container_name << "] 稀疏度汇总 ===" << std::endl;
    std::cout << "有效统计矩阵数: " << valid_matrix_count << std::endl;
    std::cout << "平均非零元素占比: " << std::fixed << std::setprecision(4) << avg_sparsity * 100 << "%" << std::endl;
    std::cout << "平均零元素占比: " << std::fixed << std::setprecision(4) << (1 - avg_sparsity) * 100 << "%" << std::endl;

    // 给出优化建议
    std::cout << "\n=== 优化建议 ===" << std::endl;
    if (avg_sparsity < 0.3) {
        std::cout << "容器矩阵稀疏度极低（<30%），建议使用Eigen稀疏矩阵（SparseMatrix）优化内存！" << std::endl;
    }
    else if (avg_sparsity < 0.7) {
        std::cout << "容器矩阵为中等稀疏（30%~70%），可权衡使用稀疏矩阵或分块存储！" << std::endl;
    }
    else {
        std::cout << "容器矩阵接近稠密（>70%），无需使用稀疏矩阵，优先按需分配/分块存储！" << std::endl;
    }

    return avg_sparsity;
}

// 快速统计非零元素数量（无输出）
int MatrixSparsityAnalyzer::fastNonZeroCount(const Eigen::MatrixXd& mat) const {
    return (mat.array().abs() > eps_).count();
}

// 设置容差
void MatrixSparsityAnalyzer::setEpsilon(double eps) {
    if (eps < 0) {
        throw std::invalid_argument("Epsilon must be non-negative!");
    }
    eps_ = eps;
}

// 获取容差
double MatrixSparsityAnalyzer::getEpsilon() const {
    return eps_;
}