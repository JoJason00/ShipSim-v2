#include <string>
#include <iostream>
#include "../src/config/CaseConfig.h"
#include "../src/io/CaseLoader.h"
#include "../src/wave/Wave.h"
//#include "../jsoncpp.cpp"
#include "../src/seakeeping/Seakeeping.h"
#include "../src/mmg/mmg.h"
#include "../src/tool/Timer.h"
#include "../src/tool/findRoot.h"
#include "../src/io/Write.h"
#include <filesystem>
#include <seakeeping/LinearCumminsTDGF.h>
#include "seakeeping/WigleyIWritePanels.h"
#include "GreenFunction/Green_compute.hpp"


int main(int argc, char* argv[])
{
	std::string run_case = (argc >= 2) ? argv[1] : "GreenFunction";
	auto root = findProjectRootFromExe(argv[0]);
	std::string casePath = (root / "cases" / run_case).string() + "/";

	//Green_compute(casePath);

	Green_integrate(casePath);

	return 0;
}