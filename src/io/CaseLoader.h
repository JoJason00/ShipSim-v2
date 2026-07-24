#pragma once

#include <string>
#include <Eigen/Dense>
#include "../json/json.h"
#include "../config/CaseConfig.h"

using ElementMatrix = Eigen::Matrix<double, 4, 3, Eigen::RowMajor>;

namespace CaseLoader 
{
    CaseConfig   loadcase(const std::string& json_file);

	//scale：单位，取1.0表示米，取0.01表示厘米，依此类推
    void UGtoElement(const std::string& datFile, const std::string& elementFile, const double scale = 1.0);

    std::unique_ptr<std::vector<ElementMatrix>> loadelement
        (std::string file, std::string ElementType, int NumElem);
};

