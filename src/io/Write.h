#pragma once

#include <string>
#include <Eigen/Dense>
#include "../config/CaseConfig.h"

namespace Write
{
	void start(const CaseConfig& cfg);

	void writefile(const std::string& filename, const Eigen::VectorXd& t, const Eigen::MatrixXd& matrix);
	void writefile(const std::string& filename, const Eigen::VectorXd& t, const Eigen::VectorXd& vector);
}