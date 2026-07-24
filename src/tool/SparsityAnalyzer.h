#pragma once

#include <Eigen/Dense>
#include <vector>
#include <string>

// 矩阵稀疏度分析工具类
class MatrixSparsityAnalyzer {
public:
    explicit MatrixSparsityAnalyzer(double eps = 1000);

    // 1. 统计单个Eigen稠密矩阵的稀疏度（返回非零元素占比，范围0~1）
    double calculateSingleMatrixSparsity(const Eigen::MatrixXd& mat,
        bool print_detail = true) const;

    // 2. 统计矩阵容器（如vector<MatrixXd>）的平均稀疏度
    //    max_count：最多统计前N个矩阵（避免输出过多，默认全部）
    double calculateMatrixContainerSparsity(const std::vector<Eigen::MatrixXd>& mat_container,
        const std::string& container_name = "MatrixContainer",
        int max_count = -1) const;

    // 3. 快速统计单个矩阵的非零元素数量（无输出，仅返回数值）
    int fastNonZeroCount(const Eigen::MatrixXd& mat) const;

    // 设置浮点容差（小于该值视为0）
    void setEpsilon(double eps);
    // 获取当前容差
    double getEpsilon() const;

private:
    double eps_; // 浮点容差，用于判断元素是否为0
};
