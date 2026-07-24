#include "Write.h"
#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "Write.h"

namespace
{
    void printSeparator(char c = '=', int len = 66) {
        std::cout << std::string(len, c) << std::endl;
    }

    std::string getCurrentTime() {
        std::time_t now = std::time(nullptr);
        std::tm tm_now{}; // 初始化tm结构体
        // 用localtime_s填充tm结构体，参数为：tm指针、time_t指针
        localtime_s(&tm_now, &now);
        std::ostringstream oss;
        oss << std::put_time(&tm_now, "%Y-%m-%d %H:%M:%S");
        return oss.str();
    }
}

void Write::start(const CaseConfig& cfg)
{
    printSeparator('*');
    std::cout << "*                                                                *" << std::endl;
    std::cout << "*         ShipSim - Ship Motion Simulation Program               *" << std::endl;
    //std::cout << "         Version: v1.0  | Start Time: " << getCurrentTime()<<std::endl;
    std::cout << "*                                                                *" << std::endl;
    std::cout << "*         Version: v1.0                                          *" << std::endl;
    std::cout << "*                                                                *" << std::endl;
    printSeparator('*');
   
    std::cout << "\n[Simulation Mode]:" << std::endl;
    
    if (cfg.enable_coupling)
        std::cout << "  - Mode: Coupled Seakeeping + Maneuvering " << std::endl;
    else if (cfg.enable_seakeeping && !cfg.enable_maneuvering)
        std::cout << "  - Mode: Seakeeping " << std::endl;
    else if (!cfg.enable_seakeeping && cfg.enable_maneuvering)
        std::cout << "  - Mode: Maneuvering " << std::endl;
    else if (cfg.enable_seakeeping && cfg.enable_maneuvering)
        std::cout << "  - Mode: Seakeeping + Maneuvering " << std::endl;
    else
        throw std::runtime_error("no simulation mode was selected!\n");

    std::cout << "\n[Simulation Ship]:\t" << cfg.Ship.Name << std::endl;

    std::cout << "\n[Status] Program initialization completed. Starting simulation...\n" << std::endl;
    printSeparator('-');
}


void Write::writefile(const std::string& filename, const Eigen::VectorXd& t, const Eigen::MatrixXd& matrix)
{
    std::ofstream out(filename);
    for (int i = 0; i < std::min(t.rows(), matrix.rows()); ++i)
    {
        out << t(i) << ",";
        for (int j = 0; j < matrix.cols(); ++j)
            out << matrix(i, j) << ",";
        out << std::endl;
    }  
    out.close();
}

void Write::writefile(const std::string& filename, const Eigen::VectorXd& t, const Eigen::VectorXd& vector)
{
    std::ofstream out(filename);
    for (int i = 0; i < std::min(t.rows(), vector.rows()); ++i)
        out << t(i) << "," << vector(i) << std::endl;
    out.close();
}
