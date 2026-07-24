#include "series.h"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "seakeeping/Greenf.h"
#include "seakeeping/Element.h"

void Green_compute(std::string filepath)
{
    double belta = 0.0;
    double miu = 0.0;
    std::array<double, 3> Time_Green;

    Greenf green;

    std::ofstream outFile1(filepath+"Green30_miu_1.csv");
    std::ofstream outFile2(filepath + "Green30_miu_2.csv");

    //计时开始
    auto start = std::chrono::high_resolution_clock::now();
    /* auto end = std::chrono::high_resolution_clock::now();
     auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);*/

    for (int i = 0; i < 21; i++)
    {

        miu = i * 0.05;
        //miu = 0.0;
           // std::cout << miu << std::endl;
           // start = std::chrono::high_resolution_clock::now();
        for (int j = 0; j < 3001; j++)
        {
            belta = j * 0.01;

            Time_Green = TDGF_ba(belta, miu);
            green.GreenFunctionCal(belta, miu);


            outFile1 << std::setprecision(15) << Time_Green[1] << ",";
            outFile2 << std::setprecision(15) << green.Gbd << ",";
        }
        outFile1 << std::endl;
        outFile2 << std::endl;
        /*end = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        std::cout << "计算时间: " << std::setprecision(15) << duration.count() << " 秒" << std::endl;*/


    }
    outFile1.close();
    outFile2.close();

    //计时结束
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    std::cout << "计算时间: " << std::setprecision(15) << duration.count() << " 秒" << std::endl;
}


void Green_integrate(std::string filepath)
{
    std::ofstream outFile(filepath + "Green_integrate_test.csv");

    double u = 2.5;
    std::vector<ElementMatrix> element_test;
    std::vector<Vector3d> point_test;

    ElementMatrix a;

    a << -0.99,	-1.0,	0.08999875-1,
        0.99498756,	-1.00000000,	-0.10999875-1,
        0.99498756,	1.00000000,	-0.10999875-1,
        -0.99498756,	1.00000000,	0.08999875-1;

    element_test.push_back(a);

	Vector3d b(-5, -2, -0.01 );
	point_test.push_back(b);

	Gsinteg gs(element_test, point_test, u);

	double dt = 0.001;   

	outFile << "Time, sG, xdG, ydG, zdG, tdG," << std::endl;
    for (int i = 0; i < 300; ++i)
    {
        double tn = i * dt;
        auto results = gs.GreenCalPanelGauss_test(tn, gs.GF1, 2);

		outFile << std::setprecision(15) << tn << ",";
        for (const auto& data : results)
        {
            outFile << std::setprecision(15) << data.sG << ","
                    << std::setprecision(15) << data.xdG << ","
                    << std::setprecision(15) << data.ydG << ","
                    << std::setprecision(15) << data.zdG << ","
                << std::setprecision(15) << data.tdG << ",";
		}   
		outFile << std::endl;
    }
	outFile.close();

	std::cout << "test Green integrate done." << std::endl; 
}