#pragma once

#include <array>

//有航速三维时域格林函数（级数展开、渐进展开求解）
//p是场点，q是源点，返回值分别为格林函数值、对各坐标的导数，及对时间的导数
extern std::array<double, 3> TDGF_ba(double belta, double miu);


