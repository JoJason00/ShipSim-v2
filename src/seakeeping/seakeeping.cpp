#include "Seakeeping.h"
#include "Gsinteg.h"
#include "../const/Const.h"
#include <math.h>
#include <fstream>
#include <iostream>
#include <functional>
#include <omp.h>
#include <iomanip>
#include <sstream>
#include <vector>
#include "../wave/RegularWave.h"
#include "../wave/IrregularWave.h"
#include "../wave/CrossWave.h"
#include "../tool/Fit.h"
#include "../tool/ParallelGuard.h"
#include "../io/Write.h"
#include "../io/FKImpulseKernelIO.h"
#include "SeakeepingDOF.h"
#include "RollDamping.h"
#include "WaveForceRegion.h"

#include <filesystem>

Seakeeping::Seakeeping(const ShipConfig& Ship, std::string casePath, const SeakeepingConfig& Seakeeping)
	:filePath(casePath), ShipCfg(Ship), SeakeepingCfg(Seakeeping),
	ak(Seakeeping.Panel.NE), rVn(Seakeeping.Panel.NE)
{
}

void Seakeeping::run()
{
	//if (SeakeepingCfg.FreeRollDecay.enabled)
	//{
	//	runFreeRollDecay();
	//	return;
	//}

	ProcessElement();
	hs = SeakeepingDOF::hydrostaticsFromWaterline(ShipCfg, SeakeepingCfg, *element);
	//initGreenTable();

	rollVisc_ = RollDampingBuilder::build(
		filePath,
		SeakeepingCfg.RollDamping,
		ShipCfg.Mass.Mass,
		ShipCfg.Mass.GM
	);

	std::cout << "Roll damping loaded:\n"
		<< "  B44_lin  = " << rollVisc_.B44_lin << "\n"
		<< "  B44_quad = " << rollVisc_.B44_quad << "\n"
		<< "  B44_cube = " << rollVisc_.B44_cube << "\n";

	solve();
}

void Seakeeping::ProcessElement()
{
	std::string file = filePath + ShipCfg.Name + ".element";

	std::unique_ptr<std::vector<ElementMatrix>> ElementData =
		CaseLoader::loadelement(file, SeakeepingCfg.Panel.NEType, SeakeepingCfg.Panel.NE);

	element = std::make_shared<Element>(SeakeepingCfg.Solver, SeakeepingCfg.Panel.NE, std::move(ElementData));

	Eigen::Vector3d cg;
	cg << ShipCfg.Mass.CG.at(0), ShipCfg.Mass.CG.at(1), ShipCfg.Mass.CG.at(2);

	element->Geometry(cg);
	element->RankineSource();

	int n_WL = element->n_WL;

	N0sq.resize(n_WL);
	PotL_idx.resize(n_WL);

	double	n0;
	int		pk;
	for (int k = 0; k < n_WL; ++k)
	{
		pk = element->PotL[k];
		PotL_idx(k) = pk;
		n0 = element->Nvec(pk, 0);
		N0sq(k) = static_cast<GScalar>(n0 * n0);
	}
}

void Seakeeping::solve()
{
	RAO4DTable tab = initRAO4D();
	int n_waves = SeakeepingCfg.waves.size();

	for (int iFn = 0; iFn < (int)SeakeepingCfg.Fn.size(); ++iFn)
	{
		double Fn = SeakeepingCfg.Fn[iFn];
		for (int i_case = 0; i_case < n_waves; ++i_case)
		{
			CaseContext ctx = buildCaseContext(i_case, Fn, iFn);

			allocCaseBuffers(ctx);

			computeGreenTables(ctx);
			computeExciting(ctx);
			integrateAndFitStore(ctx, tab);
		}
	}

	std::vector<double> non_we;
	tab.writeCSV(filePath + "RAO.csv", non_we, 12, true);
}


void Seakeeping::adaptTimeByCircle(double omega)
{
	double T = 2.0 * PI / omega;

	SeakeepingCfg.Time.dt = T / 50.0;

	std::cout << "omiga:\t" << omega << "\tT:\t" << T << "\tdt:\t" << SeakeepingCfg.Time.dt << std::endl;

	SeakeepingCfg.Time.PreStep = int(SeakeepingCfg.Time.PreCircle * T / SeakeepingCfg.Time.dt);
	int steadyStep = int(SeakeepingCfg.Time.TimeCircle * T / SeakeepingCfg.Time.dt);
	SeakeepingCfg.Time.TimeStep = SeakeepingCfg.Time.PreStep + steadyStep + 2;

	SeakeepingCfg.Time.GreenStep = int(SeakeepingCfg.Time.GreenCircle * T / SeakeepingCfg.Time.dt);

	bufCols = SeakeepingCfg.Time.GreenStep + 1;

	if (SeakeepingCfg.Time.GreenStep >= SeakeepingCfg.Time.TimeStep)
		SeakeepingCfg.Time.GreenStep = SeakeepingCfg.Time.TimeStep - 1;
}


RAO4DTable Seakeeping::initRAO4D()
{
	RAO4DTable tab;

	tab.meta.L = ShipCfg.Geometry.Length;
	tab.meta.dirUnit = RAO4DMeta::DirUnit::Rad;
	tab.meta.freqUnit = RAO4DMeta::FreqUnit::RadPerSec;
	tab.meta.freqType = RAO4DMeta::FreqType::IncidentOmega;
	tab.meta.normType = RAO4DMeta::NormType::TransA_RotkA;

	tab.modeIds = SeakeepingCfg.modes;
	tab.Fns = SeakeepingCfg.Fn;

	std::vector<double> betas, oms;
	for (const auto& w : SeakeepingCfg.waves) {
		if (auto rw = std::dynamic_pointer_cast<RegularWave>(w)) {
			betas.push_back(rw->direction());
			oms.push_back(rw->getFreq());
		}
	}

	std::sort(betas.begin(), betas.end());
	betas.erase(std::unique(betas.begin(), betas.end()), betas.end());
	std::sort(oms.begin(), oms.end());
	oms.erase(std::unique(oms.begin(), oms.end()), oms.end());

	tab.dir = betas;
	tab.w = oms;

	betaAxisDeg = betas;
	omegaAxis = oms;

	tab.finalizeAndAllocate(true);
	return tab;
}


void Seakeeping::initGreenTable()
{
	GreenTable::Params p = GreenTable::Params::AdaptiveDefault(/*bmax=*/300.0);

	p.bmin = 1e-6;
	p.bmax = 300.0;
	p.mmin = 0.0;
	p.mmax = 1.0;

	p.mSplit1 = 0.005;
	p.mSplit2 = 0.01;
	p.mSplit3 = 0.1;

	//        ġ    Nm          ߡ Nb1   ׼               У 
	// p.NmTotal = 660;              //    ϸ   800/1000
	// p.Nb1Ref = 30000; p.bRef=300; //      Nb1   ׼  ͬ ͸     
	// p.nbPow = 1.0;                // bmax 仯ʱ Nb1      ţ 1   ԣ 2ƽ  

	p.throwOnOutOfRange = true;

	p.autoTune(); // <--  ؼ        bmax/NmTotal/Nb1Ref  Զ      Nm1..4 / Nb1..4

	gGreenTable = GreenTable(p);
	gGreenTable.init(filePath + "green_table.bin", G, true);
}



//            
CaseContext Seakeeping::buildCaseContext(const int i_case, const double Fn, const int iFn)
{
	CaseContext ctx;

	ctx.i_case = i_case;
	ctx.DOF = SeakeepingCfg.DOF;
	ctx.iFn = iFn;

	ctx.wave = SeakeepingCfg.waves.at(i_case);
	ctx.reg = std::dynamic_pointer_cast<RegularWave>(ctx.wave);

	ctx.Amp = ctx.wave->getAmp();
	ctx.W = ctx.wave->getFreq();

	ctx.Fn = Fn;
	ctx.U = Fn * std::sqrt(G * ShipCfg.Geometry.Length);
	ctx.UsquareG = ctx.U * ctx.U / G;

	//        д   ʹ ó Ա     U / UsquareG  Ĵ  루 Ȳ   ģ 
	U = ctx.U;
	UsquareG = ctx.UsquareG;

	// --- encounter frequency (we) + time setting ---
	if (ctx.reg)
	{
		ctx.dirRad = ctx.reg->direction();
		//ctx.we = ctx.W - ctx.W * ctx.W * ctx.U * cos(ctx.dirRad) / G;
		auto enc = calcEncounterInfo(ctx.W, ctx.U, ctx.dirRad);
		ctx.we = enc.we;
		adaptTimeByCircle(ctx.we);                  //    ޸  config.*

		ctx.iDir = findExactIndex(betaAxisDeg, ctx.dirRad, 1e-6);
		ctx.iw = findExactIndex(omegaAxis, ctx.W, 1e-9);

	}
	else
	{
		//  ǹ  򲨣 Ҫô       Ҫô     ˻     
		ctx.dirRad = 0.0;
		ctx.we = ctx.W;
		adaptTimeByCircle(ctx.we);
	}

	std::ostringstream ss;
	ss << std::fixed << std::setprecision(2)
		<< "_Fn" << Fn;

	if (ctx.reg)
		ss << "_Dir" << std::setprecision(3) << std::round(ctx.dirRad)
		<< "_We" << std::setprecision(3) << ctx.we;
	else
		ss << "_W" << std::setprecision(3) << ctx.W;

	ctx.tag = ss.str();

	ctx.dt = SeakeepingCfg.Time.dt;
	ctx.NE = SeakeepingCfg.Panel.NE;
	ctx.TS = SeakeepingCfg.Time.TimeStep;
	ctx.tMot = SeakeepingCfg.Time.PreStep;
	ctx.TG = SeakeepingCfg.Time.GreenStep;
	ctx.n_WL = element->n_WL;

	std::cout << "NE:\t" << ctx.NE << "\tTG:\t" << ctx.TG << " \tTS:\t" << ctx.TS
		<< " \ttMot:\t" << ctx.tMot << " \tdt:\t" << ctx.dt << " \twe:\t"<<ctx.we<<std::endl;

	//ctx.a_scale[0] = rho * ShipCfg.Geometry.Displacement * G * ctx.Amp / ShipCfg.Geometry.Length;

	double L = ShipCfg.Geometry.Length;
	double displacement = ShipCfg.Geometry.Displacement;	
	double a = ctx.Amp;
	ctx.a_scale[0] = 1;
	ctx.a_scale[1] = 1;
	ctx.a_scale[2] = rho * G * displacement * a / L;
	ctx.a_scale[3] = 1;
	ctx.a_scale[4] = rho * G * a * displacement;
	ctx.a_scale[5] = 1;

	ctx.ExcitingForce.resize(ctx.TS - ctx.tMot, ctx.DOF);
	ctx.ExcitingForce.setZero();

	return ctx;
}


void Seakeeping::allocCaseBuffers(const CaseContext& ctx)
{
	//	Gz     = std::make_unique<std::vector<Eigen::MatrixXd>>(ctx.TG, Eigen::MatrixXd::Zero(ctx.NE, ctx.NE));
	//	dGz    = std::make_unique<std::vector<Eigen::MatrixXd>>(ctx.TG, Eigen::MatrixXd::Zero(ctx.NE, ctx.NE));
	//	Gz_wl  = std::make_unique<std::vector<Eigen::MatrixXd>>(ctx.TG, Eigen::MatrixXd::Zero(ctx.n_WL, ctx.NE));
	//	dGz_wl = std::make_unique<std::vector<Eigen::MatrixXd>>(ctx.TG, Eigen::MatrixXd::Zero(ctx.n_WL, ctx.NE));

	Gz = std::make_unique<std::vector<GMat>>(ctx.TG, GMat::Zero(ctx.NE, ctx.NE));
	dGz = std::make_unique<std::vector<GMat>>(ctx.TG, GMat::Zero(ctx.NE, ctx.NE));

	// ת ô棺NE x n_WL
	Gz_wlT = std::make_unique<std::vector<GMat>>(ctx.TG, GMat::Zero(ctx.NE, ctx.n_WL));
	dGz_wlT = std::make_unique<std::vector<GMat>>(ctx.TG, GMat::Zero(ctx.NE, ctx.n_WL));

	motions.resize(ctx.TS - ctx.tMot - 2, ctx.DOF);
	eforce.resize(ctx.NE, ctx.TS, ctx.TG);
	rforce.resize(ctx.NE, ctx.TS, ctx.TG);
}

void Seakeeping::computeGreenTables(const CaseContext& ctx)
{
#pragma omp parallel
	{
		shipsim::EigenSingleThreadGuard _g;
		Gsinteg green(element, U);  // ÿ ߳ һ
#pragma omp for schedule(static)
		for (int tN = 0; tN < ctx.TG; ++tN)
		{
			double tn = (tN + 1) * ctx.dt;
			GreenFunction(tN, tn, green);
		}
	}
	std::cout << "compute Green tables done." << std::endl;
}

void Seakeeping::computeExciting(CaseContext& ctx)
{
	fkpData		    fkpdata;
	FKphi		    fkphi(SeakeepingCfg.Panel.NE);
	Eigen::VectorXd t(ctx.TS);

	initialFK(fkpdata);
	ctx.wave->loadData(fkpdata);

	std::string   excitingFile = filePath + "/ExcitingForce/ExcitingForce" + ctx.tag + ".csv";
	std::ofstream efile(excitingFile);
	if (!efile.is_open())
		throw std::runtime_error("can't create exciting file!\n");

	int		tN;
	double  tn, tnM;

	double dimension_t = sqrt(G / ShipCfg.Geometry.Length);
	for (tN = 0, tn = ctx.dt, tnM = -ctx.tMot * ctx.dt; tN < ctx.TS; tN++, tn += ctx.dt, tnM += ctx.dt)
	{
		//std::cout << "Exctied force tN:\t" << tN << std::endl;

		ctx.wave->Exciting(tnM, fkphi);
		ExcitingCal(tN, tn, fkphi, ctx);

		t(tN) = tnM * dimension_t;
		efile << tnM * dimension_t << ",";
		for (int i = 0; i < 6; ++i)
		{
			efile << eforce.force(i) / ctx.a_scale[i] << ",";
		}
		efile << std::endl;

		if (tN >= ctx.tMot)
		{
			for (int i = 0; i < ctx.DOF; i++)
			{
				ctx.ExcitingForce(tN - ctx.tMot, i) = eforce.force(SeakeepingCfg.modes[i]);
			}
		}
	}
	std::cout << "compute Exciting force done." << std::endl;

	writeFKImplese(filePath,ctx);
}

void Seakeeping::integrateAndFitStore(const CaseContext& ctx, RAO4DTable& tab)
{
	auto* regWave = dynamic_cast<RegularWave*>(ctx.wave.get());
	//double we = ctx.W - ctx.W * ctx.W * ctx.U * cos(regWave->direction()) / G;
	double we = ctx.we;
	Eigen::MatrixXd Force(6, ctx.TS);

	int tN;
	double tn;
	int DOF = SeakeepingCfg.DOF;
	for (int i = 0; i < DOF; ++i)
	{
		int mode = SeakeepingCfg.modes[i];

		for (tN = 0, tn = ctx.dt; tN < ctx.TS; tN++, tn += ctx.dt)
		{
			addedBoundary(tn, mode, ctx.Amp, we);
			RadiationCal(tN, tn, rVn);
			Force.col(tN) = rforce.force;
		}
		switch (mode)
		{
		case 2:
			AddedMassAndDamping(mode, Force, 2, we, ctx.Amp, added.a33, added.b33);
			AddedMassAndDamping(mode, Force, 3, we, ctx.Amp, added.a43, added.b43);
			AddedMassAndDamping(mode, Force, 4, we, ctx.Amp, added.a53, added.b53);
			break;
		case 3:
			AddedMassAndDamping(mode, Force, 2, we, ctx.Amp, added.a34, added.b34);
			AddedMassAndDamping(mode, Force, 3, we, ctx.Amp, added.a44, added.b44);
			AddedMassAndDamping(mode, Force, 4, we, ctx.Amp, added.a54, added.b54);
			break;
		case 4:
			AddedMassAndDamping(mode, Force, 2, we, ctx.Amp, added.a35, added.b35);
			AddedMassAndDamping(mode, Force, 3, we, ctx.Amp, added.a45, added.b45);
			AddedMassAndDamping(mode, Force, 4, we, ctx.Amp, added.a55, added.b55);
			break;
		}
	}
	double L = ShipCfg.Geometry.Length;

	//std::string outK30 = filePath + "K507_Fn" + std::to_string(Fn) + ".csv";
	//double afk = rho * G * L * L * L;
	//Eigen::VectorXd K30 = (eforce.fkForce.col(4) + eforce.dForce.col(4)) / afk;
	//Write::writefile(outK30, t, K30);


	//Fit fitExcit;
	//std::vector<std::string> name = { "    ","  ҡ" };
	//ExcitingForce.col(0) /= a;
	//ExcitingForce.col(1) /= (a*L);
	//fitExcit.setdata(ExcitingForce, dt, we, Amp, name);
	////fitExcit.setdata(eforce.fkForce, dt, we, Amp, name);
	//fitExcit.run();
	//std::string outExcit = filePath + "Excit_Fn" + std::to_string(Fn) + ".csv";
	//fitExcit.writeFile(L, outExcit);

	//added.a33 = 0.4515 * rho * config.displacement;
	//added.a55 = 0.0215 * rho * config.displacement * L * L;
	//added.a35 = -0.00007 * rho * config.displacement * L;
	//added.a53 = -added.a35;
	//added.b33 = 1.285 * rho * config.displacement * sqrt(G / L);
	//added.b55 = 0.049 * rho * config.displacement * L * L * sqrt(G / L);
	//added.b35 = 0.06 * rho * config.displacement * L * sqrt(G / L);
	//added.b53 = -added.b35;

	/*3.29942, 0.543041, 1.35336
		3.29942, -2.9228e-05, 6.00343e-06
		3.29942, 0.00140687, -0.034425
		3.29942, -0.000214087, 0.0792156
		3.29942, 9.24798e-05, -0.000142444
		3.29942, 0.0246141, 0.0444241*/

	//double displacement = ShipCfg.Geometry.Displacement;
	//added.a33 = 0.543041 * rho * displacement;
	//added.a55 = 0.0246141 * rho * displacement * L * L;
	//added.a35 = 0.00140687 * rho * displacement * L;
	//added.a53 = -added.a35;
	//added.b33 = 1.35336 * rho * displacement * sqrt(G / L);
	//added.b55 = 0.0444241 * rho * displacement * L * L * sqrt(G / L);
	//added.b35 = -0.034425 * rho * displacement * L * sqrt(G / L);
	//added.b53 = -added.b35;


	Eigen::MatrixXd M(2, 2);
	Eigen::MatrixXd b(2, 2);
	Eigen::MatrixXd c(2, 2);
	M << ShipCfg.Mass.Mass + added.a33, added.a35,
		added.a53, ShipCfg.Mass.Iyy + added.a55;
	b << added.b33, added.b35,
		added.b53, added.b55;
	c << hs.c33, hs.c35,
		hs.c53, hs.c55;

	//Eigen::MatrixXd M(3, 3);
	//Eigen::MatrixXd b(3, 3);
	//Eigen::MatrixXd c(3, 3);
	//M << ShipCfg.Mass.Mass + added.a33, added.a34, added.a35,
	//	added.a43, ShipCfg.Mass.Ixx + added.a44, added.a45,
	//	added.a53, added.a54, ShipCfg.Mass.Iyy + added.a55;
	//b << added.b33, added.b34, added.b35,
	//	added.b43, added.b44, added.b45,
	//	added.b53, added.b54, added.b55;
	//c << hs.c33, hs.c34, hs.c35,
	//	hs.c43, hs.c44, hs.c45,
	//	hs.c53, hs.c54, hs.c55;

		//Eigen::MatrixXd M, b, c;
		//SeakeepingDOF::buildSystemMatrices(config, hs, M, b, c);
	Eigen::MatrixXd Minv = M.inverse();

	//int DOF = ctx.DOF;
	ODEFunction ship_motion = [this, &Minv, &b, &c, &DOF](double t, const std::vector<double>& y, double alpha, const Eigen::MatrixXd& totalForce)
	{
		std::vector<double> dy(2 * DOF, 0.0);
		Eigen::VectorXd f_interp(DOF);

		for (int i = 0; i < DOF; ++i)
			f_interp[i] = totalForce(0, i) * (1.0 - alpha) + totalForce(1, i) * alpha;

		Eigen::VectorXd vel(DOF), disp(DOF);
		for (int i = 0; i < DOF; ++i) {
			disp(i) = y.at(i);
			vel(i) = y.at(DOF + i);
		}
		Eigen::VectorXd acc = Minv * (f_interp - b * vel - c * disp);

		for (int i = 0; i < DOF; ++i) {
			dy.at(i) = vel(i);
			dy.at(DOF + i) = acc(i);
		}
		return dy;
	};

	std::vector<double> y0(2 * DOF, 0.0);
	double t0 = 0.0;
	double t_end = (ctx.TS - ctx.tMot - 2) * ctx.dt;

	rk4_solve(ship_motion, ctx.ExcitingForce, y0, t0, t_end, ctx.dt);

	const double rotScale = G / (ctx.Amp * ctx.W * ctx.W);
	for (int k = 0; k < ctx.DOF; ++k) {
		int modeId = SeakeepingCfg.modes[k];
		if (modeId <= 2) motions.col(k) /= ctx.Amp;      // x,y,z
		else             motions.col(k) *= rotScale;     // roll,pitch,yaw
	}

	std::ofstream out(filePath + "motion/motion_Fn" + ctx.tag + ".csv");
	out << "t_nd";
	for (int k = 0; k < ctx.DOF; ++k) out << ",mode" << SeakeepingCfg.modes[k];
	out << "\n";

	for (int i = 0; i < motions.rows(); ++i) {
		out << (i + 1) * ctx.dt * sqrt(G / L);
		for (int k = 0; k < ctx.DOF; ++k) out << "," << motions(i, k);
		out << "\n";
	}

	Fit fitRAO;
	auto names = SeakeepingDOF::modeNames(SeakeepingCfg);
	fitRAO.setdata(motions, ctx.dt, ctx.we, 1.0, names);
	fitRAO.run();

	std::vector<std::complex<double>> Hk;
	Hk.reserve(ctx.DOF);
	for (auto& r : fitRAO.results()) {
		Hk.emplace_back(r.amplitude * std::cos(r.phase),
			r.amplitude * std::sin(r.phase));
	}
	tab.setAt(ctx.iFn, ctx.iDir, ctx.iw, Hk);
}



int Seakeeping::findExactIndex(const std::vector<double>& xs, double x, double tol)
{
	auto it = std::lower_bound(xs.begin(), xs.end(), x);
	if (it != xs.end() && std::abs(*it - x) <= tol) return (int)std::distance(xs.begin(), it);
	if (it != xs.begin()) {
		auto it2 = it - 1;
		if (std::abs(*it2 - x) <= tol) return (int)std::distance(xs.begin(), it2);
	}
	throw std::runtime_error("value not on axis (grid incomplete or unit mismatch)");
}



void Seakeeping::GreenFunction(int tN, double tn, Gsinteg& green)
{
	double tmp = SeakeepingCfg.Time.dt;
	int	   NE = SeakeepingCfg.Panel.NE;
	int    n_wl = element->n_WL;

	GreenData Gdata;
	double nV;
	double zg, zd;

	int p;

	for (int j = 0; j < NE; j++) {                    //j ǳ   Ԫ
		for (int i = 0; i < NE; i++) {			      //i  Դ  Ԫ

			//      άʱ    ֺ   
			Gdata = green.GreenCal(j, i, tn, green.GF1, gGreenTable);

			if (SeakeepingCfg.Solver == "Potential")
			{
				//  άʱ    ֺ     i  Ԫ      
				nV = -(element->A31[i] * Gdata.xdG + element->A32[i] * Gdata.ydG + element->A33[i] * Gdata.zdG);
			}
			else
			{
				//  άʱ    ֺ     j  Ԫ      
				nV = -(element->A31[j] * Gdata.xdG + element->A32[j] * Gdata.ydG + element->A33[j] * Gdata.zdG);
			}

			zg = Gdata.sG * tmp;

			//nV     ֺ       Ԫ      
			zd = nV * tmp;

			Gz->at(tN)(j, i) = static_cast<GScalar>(zg);
			dGz->at(tN)(j, i) = static_cast<GScalar>(zd);
		}

		//  ˮ  Ԫ ϵĸ  ֺ       	
		for (int i = 0; i < n_wl; i++)
		{
			Gdata = green.GreenCal_WL(j, i, tn, green.GF1, gGreenTable);
			if (SeakeepingCfg.Solver == "Potential")
			{
				//  άʱ       i  Ԫ      
				// nV=(-A31[i]*xdG-A32[i]*ydG+A33[i]*zdG);	
				//nV = (-(element->ypl(i,0) - element->ypl(i, 1)) * Gdata.xdG - 
				//	(element->xpl(i, 0) - element->xpl(i, 1)) * Gdata.ydG) / 
				//	sqrt(pow(element->ypl(i, 0) - element->ypl(i, 1), 2) + 
				//		element->ypl(i, 0) - element->ypl(i, 1));
				nV = -(element->A31[j] * Gdata.xdG + element->A32[j] * Gdata.ydG + element->A33[j] * Gdata.zdG);
			}
			else
			{
				//p = PotL_idx(j);
				//  άʱ       j  Ԫ      
				nV = -(element->A31[j] * Gdata.xdG + element->A32[j] * Gdata.ydG + element->A33[j] * Gdata.zdG);
			}

			zg = Gdata.sG * tmp;
			zd = nV * tmp;

			Gz_wlT->at(tN)(j, i) = static_cast<GScalar>(zg);
			dGz_wlT->at(tN)(j, i) = static_cast<GScalar>(zd);
		}
	}
}

void Seakeeping::initialFK(fkpData& fkpdata)
{
	fkpdata.NE = SeakeepingCfg.Panel.NE;
	fkpdata.U = U;
	fkpdata.A31 = element->A31;
	fkpdata.A32 = element->A32;
	fkpdata.A33 = element->A33;
	fkpdata.xcr = element->xcr;
	fkpdata.ycr = element->ycr;
	fkpdata.zcr = element->zcr;
}


void Seakeeping::ExcitingCal(int tN, double tn, FKphi& fkphi, const CaseContext& ctx)
{
	int j, tn1;
	double ft1, temp = 0.0;
	auto &wave = ctx.wave;
	int tMot = ctx.tMot;

	if (!(element->method))
	{
		//PotentialConvolution(tN, Kpz, veSg);
		//nForceCalculation(tN, veSg[tN], eFt, eForce);
		//temp = 1.0;

		SourceConvolution(tN, fkphi.df, eforce.veSg_d, eforce.veSg_g, eforce.sPot);
		ForceCal(tN, eforce.sPot, eforce.eFdt, eforce.force);
		temp = -1.0;                             //method==1    Ԫ      ָ       ⣬      
	}
	else
	{
		//         
		SourceConvolution(tN, fkphi.df, eforce.veSg_d, eforce.veSg_g, eforce.sPot);
		ForceCal(tN, eforce.sPot, eforce.eFdt, eforce.force);
		temp = -1.0;                             //method==1    Ԫ      ָ       ⣬      
	}
	//        δ ˲  ߣ 
	eforce.dForce.row(tN) = eforce.force;

	double dt = SeakeepingCfg.Time.dt;
	for (j = 0; j < 6; ++j) {

		// FK  û г˲  ߵĲ   
		eforce.fkForce(tN, j) = temp * element->ArInt.col(j).dot(fkphi.fk) * rho;

		ft1 = 0.0;

		// if(tN>=tMot-1)
		//int tn1_min = std::max(0, tN - 2 * ctx.tMot);
		for (tn1 = 0; tn1 < tN; tn1++)
			//   Բ  ߣ     ʱ       ֣ Eta ǲ  ߣ ά ȣ Ƶ    *ʱ 䲽  
			ft1 += (eforce.dForce(tN - tn1, j) + eforce.fkForce(tN - tn1, j)) * wave->Eta(tn1 * dt) * dt;

		eforce.force(j) = ft1;
	}
}

void Seakeeping::writeFKImplese(const std::string filePath, const CaseContext& ctx)
{
	const WaveForceRegion region =
		wave_force_region::classify(ctx.U, ctx.dirRad, ctx.W);

	const std::filesystem::path dir = std::filesystem::path(filePath) / "fkImpulse";
	std::filesystem::create_directories(dir);

	std::ostringstream fname;
	fname << "fkImpulse_nd_";
	if (SeakeepingCfg.FKImpulseKernel.IncludeShipName && !ShipCfg.Name.empty())
		fname << ShipCfg.Name << "_";
	fname << "Fn"   << FKImpulseKernelIO::keyDouble(ctx.Fn)
	      << "_U"   << FKImpulseKernelIO::keyDouble(ctx.U)
	      << "_Dir" << FKImpulseKernelIO::keyDouble(ctx.dirRad)
	      << "_Bkt" << wave_force_region::impulseBucketTag(region)
	      << ".csv";

	const std::string outFile = (dir / fname.str()).string();
	std::ofstream out(outFile);
	if (!out.is_open())
	{
		std::cerr << "writeFKImplese: cannot open " << outFile << "\n";
		return;
	}

	const double L = ShipCfg.Geometry.Length;
	const double displacement = ShipCfg.Geometry.Displacement;

	const double t_scale = std::sqrt(G / L);

	// K30/K37 and K50/K57: project scaling; K30_plus_K37 / K50_plus_K57: literature sum scale.
	const double scaleK3 = rho * G * displacement / L * std::sqrt(G / L); // K30, K37
	const double scaleK5 = rho * G * displacement * std::sqrt(G / L);   // K50, K57

	const double scaleK3sum = rho * G * L * L;
	const double scaleK5sum = rho * G * L * L * L;

	out << std::setprecision(17);
	out << "t_nd,K30,K37,K50,K57,K30_plus_K37,K50_plus_K57\n";

	for (int i = 0; i < eforce.fkForce.rows(); ++i)
	{
		const double tnd = (i - ctx.tMot) * SeakeepingCfg.Time.dt * t_scale;

		const double K30 = eforce.fkForce(i, 2) / scaleK3;
		const double K37 = eforce.dForce(i, 2) / scaleK3;

		const double K50 = eforce.fkForce(i, 4) / scaleK5;
		const double K57 = eforce.dForce(i, 4) / scaleK5;

		const double K3sum = (eforce.fkForce(i, 2) + eforce.dForce(i, 2)) / scaleK3sum;
		const double K5sum = (eforce.fkForce(i, 4) + eforce.dForce(i, 4)) / scaleK5sum;

		out << tnd << ","
			<< K30 << "," << K37 << ","
			<< K50 << "," << K57 << ","
			<< K3sum << "," << K5sum << "\n";
	}
	std::cout << "write FK impulse (nondim) done: " << outFile << std::endl;
}

void Seakeeping::addedBoundary(double tn, int mode, double Amp, double we)
{
	double vn = -Amp * we * sin(we * tn);
	switch (mode) {
	case 1: case 2: case 3:
		rVn = vn * element->Nvec.col(mode);
		break;

	case 4:
		rVn = vn * element->Nvec.col(mode) + U * Amp * cos(we * tn) * element->Nvec.col(2);
		break;
	}
}

void Seakeeping::RadiationCal(int tN, double tn, const Eigen::VectorXd rVn)
{
	if (SeakeepingCfg.Solver == "Potential") {
		////     ٶ   
		//PotentialConvolution(tN, rVn, vSg);
		////      	
		//nForceCalculation(tN, vSg[tN], rFt, rForce);
	}
	else {
		// ֲ Դ       ٶ  ƣ ԴǿvSg(TS,NE), ٶ  ƣ sPot(NE)
		SourceConvolution(tN, rVn, rforce.vSg_d, rforce.vSg_g, rforce.sPot);
		//     ٶ  ƵĻ   rFt   䵼       ַ     rForce
		ForceCal(tN, rforce.sPot, rforce.rFdt, rforce.force);
	}
}

void Seakeeping::AddedMassAndDamping(int mode, Eigen::MatrixXd& Force, int n, double we, double Amp, double& addedMass, double& dampingCoefficient)
{
	//   㸽       ͸       ϵ  
	addedMass = 0;
	dampingCoefficient = 0;

	double nond_addMass = 0;
	double nond_dampingCoefficient = 0;

	//tΪ   ֵĿ ʼʱ  
	 //double t = 50*dt;
	//int i0 = 200;
	int i0 = 100;
	double dt = SeakeepingCfg.Time.dt;
	double t = i0 * dt;
	double tn = t;
	double T = 2 * PI / we;


	while (tn < SeakeepingCfg.Time.TimeStep * dt)
	{
		tn += T;
	}
	tn -= T;

	// double t=0;
	int i = i0;
	while (t < tn)
	{
		addedMass += Force(n, i) * cos(we * (t + dt / 2)) * dt;
		dampingCoefficient += Force(n, i) * sin(we * (t + dt / 2)) * dt;
		t += dt;
		i++;
	}

	addedMass *= (2 / (tn - i0 * dt));
	dampingCoefficient *= (2 / (tn - i0 * dt));

	addedMass /= (we * we * Amp);
	dampingCoefficient /= (we * Amp);

	double L = ShipCfg.Geometry.Length;
	double displacement = ShipCfg.Geometry.Displacement;
	if (mode == 2)
	{
		switch (n)
		{
		case 2: nond_addMass = addedMass / (displacement * rho);
			nond_dampingCoefficient = dampingCoefficient / (displacement * rho * sqrt(G / L));
			break;
		case 3: nond_addMass = addedMass;
			nond_dampingCoefficient = dampingCoefficient;
			break;
		case 4: nond_addMass = addedMass / (displacement * rho * L);
			nond_dampingCoefficient = dampingCoefficient / (displacement * rho * sqrt(G / L) * L);
			break;
		}
	}
	else if (mode == 4)
	{
		switch (n)
		{
		case 2: nond_addMass = addedMass / (displacement * rho * L);
			nond_dampingCoefficient = dampingCoefficient / (displacement * rho * sqrt(G / L) * L);
			break;
		case 3: nond_addMass = addedMass;
			nond_dampingCoefficient = dampingCoefficient;
			break;
		case 4: nond_addMass = addedMass / (displacement * rho * L * L);
			nond_dampingCoefficient = dampingCoefficient / (displacement * rho * sqrt(G / L) * L * L);
			break;
		}
	}

	else if (mode == 3)
	{
		nond_addMass = addedMass;
		nond_dampingCoefficient = dampingCoefficient;
	}

	std::string   added = filePath + "added.txt";
	std::ofstream addedfile(added, std::ios::app);
	addedfile << we * sqrt(L / G) << "," << nond_addMass << "," << nond_dampingCoefficient << std::endl;
	addedfile.close();
}

void Seakeeping::SourceConvolution(int& tN, const Eigen::VectorXd& Vn,
	Eigen::MatrixXd& Sg_d, SgMatG& Sg_g, Eigen::VectorXd& sPot)
{
	const int NE = SeakeepingCfg.Panel.NE;
	const int nWL = element->n_WL;
	const int t0 = (tN > SeakeepingCfg.Time.GreenStep) ? (tN - SeakeepingCfg.Time.GreenStep) : 0;

	ak = 4 * PI * Vn;

	const int nth = omp_get_max_threads();
	std::vector<Eigen::VectorXd, Eigen::aligned_allocator<Eigen::VectorXd>> partial(nth);
	for (auto& v : partial) v = Eigen::VectorXd::Zero(NE);

#pragma omp parallel
	{
		shipsim::EigenSingleThreadGuard _g;
		const int tid = omp_get_thread_num();
		Eigen::VectorXd& acc = partial[tid];

		//  ߳ ˽    ʱ    һ η  䣬ѭ     ã
		GVec sg_g(NE);
		GVec wl_weight_g(nWL);
		GVec tmpNE1(NE), tmpNE2(NE); // float  м         ÿ   cast       ʱ

#pragma omp for schedule(static)
		for (int tn = t0; tn < tN; ++tn)
		{
			const int t_idx = tN - tn - 1;

			//dGz_wl_matά ȣ NE    n_WL   dGz_matά ȣ NE    NE
			const auto& dGz_wlT_m = dGz_wlT->at(t_idx);
			const auto& dGz_mat = dGz->at(t_idx);

			sg_g.noalias() = Sg_g.col(col(tn));

			wl_weight_g = sg_g(PotL_idx).cwiseProduct(N0sq);

			tmpNE1.noalias() = (dGz_wlT_m * wl_weight_g);
			tmpNE2.noalias() = (dGz_mat * sg_g);

			//origin
			//acc.noalias() += UsquareG * tmpNE1.cast<double>();
			//acc.noalias() -= tmpNE2.cast<double>();

			acc.noalias() += UsquareG * tmpNE1.cast<double>();
			acc.noalias() -= tmpNE2.cast<double>();
		}
	}

	for (auto& v : partial) ak += v;

	// LU.solve    double
	Eigen::VectorXd x = element->lu.solve(ak);

	// д       ring buffer
	Sg_d.col(col(tN)) = x;
	Sg_g.col(col(tN)) = x.cast<GScalar>();

	// calculate potential
	SourcePotential(t0, tN, sPot, Sg_d, Sg_g);
}


void Seakeeping::ForceCal(int& tN, Eigen::VectorXd& Pt, Eigen::RowVectorXd& Fdt, Eigen::RowVectorXd& F)
{
	int k;
	double force, derivativeofforce;

	for (k = 0; k < 6; ++k) {

		//ArInt  Ԫ             6*NE  Pt: ٶ   
		force = -element->ArInt.col(k).dot(Pt);

		// ٶ  ƶ ʱ 䵼    rFdt: 洢  һ  ʱ 䲽   ٶ  ƻ   
		derivativeofforce = (force - Fdt[k]) / SeakeepingCfg.Time.dt;
		Fdt[k] = force;

		switch (k) {
		case 0: case 1: case 2: case 3:
			F[k] = derivativeofforce * rho;
			break;
		case 4:
			F[k] = (derivativeofforce + U * element->ArInt.col(2).dot(Pt)) * rho;
			break;
		case 5:
			F[k] = (derivativeofforce - U * element->ArInt.col(1).dot(Pt)) * rho;
			break;
		}
	}
}

void Seakeeping::SourcePotential(int t0, int tN, Eigen::VectorXd& Pt,
	const Eigen::MatrixXd& Sg_d, const SgMatG& Sg_g)
{
	int		tn;
	int		t_idx;

	Pt = element->Rz * Sg_d.col(col(tN));

	GVec sg_g(SeakeepingCfg.Panel.NE);
	GVec wl_weight_g(element->n_WL);

	for (tn = t0; tn < tN; ++tn)
	{
		t_idx = tN - tn - 1;

		const auto& Gz_mat = Gz->at(t_idx);
		const auto& Gz_wlT_m = Gz_wlT->at(t_idx);

		sg_g.noalias() = Sg_g.col(col(tn));
		wl_weight_g = sg_g(PotL_idx).cwiseProduct(N0sq);

		Pt.noalias() += (Gz_mat * sg_g).cast<double>();
		Pt.noalias() -= UsquareG * (Gz_wlT_m * wl_weight_g).cast<double>();
	}

	//      ٶ   ʵ  ȡ ˸  ţ   Ϊ          ʱ       ķ     
	Pt /= (4 * PI);
	//double coef = Usquare / (4 * PI * G);
	//Pt *= coef;
}

void Seakeeping::rk4_solve(
	const ODEFunction& f, const Eigen::MatrixXd& ExcitingForce,
	const std::vector<double>& y0, double t0, double t_end, double h)
{
	//  Զ               Ŀ  ʱ  С ڳ ʼʱ 䣩
	if ((t_end - t0) * h < 0.0) {
		h = -h;
	}

	const double eps = 1e-5;
	if (std::abs(t0 - t_end) < eps) {
		return;
	}

	double t = t0;
	std::vector<double> y = y0;

	//   ʼ    ʱ              һ Σ 
	std::vector<double> k1, k2, k3, k4, y_temp;

	int tN = 0;

	int DOF = ExcitingForce.cols();
	Eigen::MatrixXd totalForce = Eigen::MatrixXd::Zero(2, DOF);
	//totalForce(0, 0) = ExcitingForce(0, 0);
	//totalForce(0, 1) = ExcitingForce(0, 1);

	const RollViscDamping& rvd = rollVisc_;
	while (true) {
		//    㵱ǰ         һ      С  h     ⳬ  Ŀ  ʱ 䣩
		double current_h = h;
		if ((h > 0.0f && t + h > t_end) || (h < 0.0f && t + h < t_end))
			current_h = t_end - t;

		y = rk4_step(f, ExcitingForce, tN, t, y, current_h, k1, k2, k3, k4, y_temp, totalForce, rvd);

		for (int i = 0; i < SeakeepingCfg.DOF; ++i)
			motions(tN, i) = y[i];

		t += current_h;
		tN++;

		if (std::abs(t - t_end) < eps)
			break;
	}
}

std::vector<double> Seakeeping::rk4_step(
	const ODEFunction& f, const Eigen::MatrixXd& ExcitingForce, int tN, double t0,
	const std::vector<double>& y0, double h,
	std::vector<double>& k1, std::vector<double>& k2,
	std::vector<double>& k3, std::vector<double>& k4,
	std::vector<double>& y_temp, Eigen::MatrixXd& totalForce, const RollViscDamping& vd)
{
	const int n = y0.size();

	//rVn = y0[2] * element->Nvec.col(2) + y0[3] * element->Nvec.col(4) + U * y0[1] * element->Nvec.col(2);

	//SeakeepingDOF::buildRadiationVn(config, element, y0, U, rVn);

	//RadiationCal(tN, t0, rVn);

	for (int i = 0; i < SeakeepingCfg.DOF; ++i) {

		//double rad = rforce.force(config.modes[i]);
		double rad = 0.0;
		totalForce(0, i) = totalForce(1, i);
		totalForce(1, i) = ExcitingForce(tN, i) + rad;
	}

	int iR = SeakeepingDOF::findModeIndex(SeakeepingCfg, MODE_ROLL);
	if (iR >= 0)
	{
		const double phidot = y0[SeakeepingCfg.DOF + iR];
		//totalForce(1, iR) += -vd.B44_lin * phidot - vd.B44_quad * phidot * abs(phidot) + vd.B44_cube * phidot * phidot * phidot;
		totalForce(1, iR) += vd.moment(phidot);
	}


	//             С  ƥ  ״̬ά  
	auto resize_buffer = [n](std::vector<double>& v) {
		if (v.size() != n) {
			v.resize(n);
		}
	};
	resize_buffer(k1);
	resize_buffer(k2);
	resize_buffer(k3);
	resize_buffer(k4);
	resize_buffer(y_temp);

	//     k1 = f(t0, y0)
	k1 = f(t0, y0, 0, totalForce);
	if (k1.size() != n) {
		throw std::runtime_error("k1ά    ״̬  ƥ  ");
	}

	//     k2 = f(t0 + h/2, y0 + h/2 * k1)
	for (size_t i = 0; i < n; ++i) {
		y_temp[i] = y0[i] + h * k1[i] / 2.0f;  //         ʹ  double
	}
	k2 = f(t0 + h / 2.0f, y_temp, 0.5, totalForce);
	if (k2.size() != n) {
		throw std::runtime_error("k2ά    ״̬  ƥ  ");
	}

	//     k3 = f(t0 + h/2, y0 + h/2 * k2)
	for (size_t i = 0; i < n; ++i) {
		y_temp[i] = y0[i] + h * k2[i] / 2.0f;
	}
	k3 = f(t0 + h / 2.0f, y_temp, 0.5, totalForce);
	if (k3.size() != n) {
		throw std::runtime_error("k3ά    ״̬  ƥ  ");
	}

	//     k4 = f(t0 + h, y0 + h * k3)
	for (size_t i = 0; i < n; ++i) {
		y_temp[i] = y0[i] + h * k3[i];
	}
	k4 = f(t0 + h, y_temp, 1, totalForce);
	if (k4.size() != n) {
		throw std::runtime_error("k4ά    ״̬  ƥ  ");
	}

	//       һ  ״̬  y0 + h/6*(k1 + 2*k2 + 2*k3 + k4)
	std::vector<double> y_next(n);
	for (size_t i = 0; i < n; ++i) {
		y_next[i] = y0[i] + h * (k1[i] + 2.0f * k2[i] + 2.0f * k3[i] + k4[i]) / 6.0f;
	}

	return y_next;
}

//GreenfData Seakeeping::moveGreenf()
//{
//	return GreenfData{
//			std::move(Gz),
//			std::move(dGz),
//			std::move(Gz_wl),
//			std::move(dGz_wl)
//	};
//}



//void Seakeeping::runFreeRollDecay()
//{
//	std::cout << "==== Free roll decay mode ====\n";
//
//	// 1) 先建立横摇阻尼
//	double Ieff = ShipCfg.Mass.Ixx; // 默认直接用船体横摇转动惯量
//
//	if (SeakeepingCfg.RollDamping.enabled &&
//		SeakeepingCfg.RollDamping.mode == RollDampingMode::FromDecayCsv)
//	{
//		auto samples = RollDampingBuilder::readDecayCsv(
//			filePath + SeakeepingCfg.RollDamping.decay.csvPath,
//			SeakeepingCfg.RollDamping.decay.angleInDeg
//		);
//
//		auto result = RollDampingBuilder::identifyFromDecay(
//			samples,
//			ShipCfg.Mass.Mass,
//			ShipCfg.Mass.GM,
//			SeakeepingCfg.RollDamping.decay.polyOrder,
//			SeakeepingCfg.RollDamping.decay.minPeakGap,
//			SeakeepingCfg.RollDamping.decay.IeffOverride
//		);
//
//		rollVisc_ = result.damping;
//		Ieff = result.Ieff;
//
//		std::cout << "Roll damping identified from CSV:\n"
//			<< "  omega_n   = " << result.poly.omega_n << "\n"
//			<< "  T_n       = " << result.poly.T_n << "\n"
//			<< "  Ieff      = " << result.Ieff << "\n"
//			<< "  B44_lin   = " << rollVisc_.B44_lin << "\n"
//			<< "  B44_quad  = " << rollVisc_.B44_quad << "\n"
//			<< "  B44_cube  = " << rollVisc_.B44_cube << "\n";
//	}
//	else
//	{
//		rollVisc_ = RollDampingBuilder::build(
//
//			SeakeepingCfg.RollDamping,
//			ShipCfg.Mass.Mass,
//			ShipCfg.Mass.GM
//		);
//
//		std::cout << "Roll damping loaded directly:\n"
//			<< "  Ieff      = " << Ieff << "\n"
//			<< "  B44_lin   = " << rollVisc_.B44_lin << "\n"
//			<< "  B44_quad  = " << rollVisc_.B44_quad << "\n"
//			<< "  B44_cube  = " << rollVisc_.B44_cube << "\n";
//	}
//
//	// 2) 线性恢复矩系数
//	// C44 = m * g * GM
//	const double C44 = ShipCfg.Mass.Mass * G * ShipCfg.Mass.GM;
//
//	// 3) 初始条件
//	const double dt = SeakeepingCfg.FreeRollDecay.dt;
//	const double duration = SeakeepingCfg.FreeRollDecay.duration;
//
//	double phi = SeakeepingCfg.FreeRollDecay.phi0_deg * PI / 180.0;
//	double phidot = SeakeepingCfg.FreeRollDecay.phidot0_deg_s * PI / 180.0;
//
//	const int nStep = static_cast<int>(duration / dt) + 1;
//
//	const double Tn = 2.0 * PI * std::sqrt(Ieff / C44);
//	const double omega_n = std::sqrt(C44 / Ieff);
//
//	std::cout << "Free roll decay setup:\n"
//		<< "  phi0 [deg]   = " << SeakeepingCfg.FreeRollDecay.phi0_deg << "\n"
//		<< "  phidot0 [deg/s] = " << SeakeepingCfg.FreeRollDecay.phidot0_deg_s << "\n"
//		<< "  dt [s]       = " << dt << "\n"
//		<< "  duration [s] = " << duration << "\n"
//		<< "  C44          = " << C44 << "\n"
//		<< "  Ieff         = " << Ieff << "\n"
//		<< "  omega_n      = " << omega_n << "\n"
//		<< "  Tn [s]       = " << Tn << "\n";
//
//	// 4) 输出文件
//	std::string outFile = filePath + "free_roll_decay.csv";
//	std::ofstream out(outFile);
//	if (!out.is_open())
//		throw std::runtime_error("cannot create free_roll_decay.csv");
//
//	out << "t,phi_rad,phi_deg,phidot_rad_s,phidot_deg_s,M_visc\n";
//
//	// 5) 定义 1DOF 横摇方程
//	auto rhs = [&](double phi_now, double phidot_now)
//	{
//		const double Mvisc = rollVisc_.moment(phidot_now);
//		const double phiddot = (Mvisc - C44 * phi_now) / Ieff;
//		return std::pair<double, double>{ phidot_now, phiddot };
//	};
//
//	// 6) RK4 时间积分
//	double t = 0.0;
//	for (int i = 0; i < nStep; ++i)
//	{
//		const double Mvisc = rollVisc_.moment(phidot);
//
//		out << t << ","
//			<< phi << ","
//			<< phi * 180.0 / PI << ","
//			<< phidot << ","
//			<< phidot * 180.0 / PI << ","
//			<< Mvisc << "\n";
//
//		auto [k1_phi, k1_phidot] = rhs(phi, phidot);
//		auto [k2_phi, k2_phidot] = rhs(
//			phi + 0.5 * dt * k1_phi,
//			phidot + 0.5 * dt * k1_phidot
//		);
//		auto [k3_phi, k3_phidot] = rhs(
//			phi + 0.5 * dt * k2_phi,
//			phidot + 0.5 * dt * k2_phidot
//		);
//		auto [k4_phi, k4_phidot] = rhs(
//			phi + dt * k3_phi,
//			phidot + dt * k3_phidot
//		);
//
//		phi += dt * (k1_phi + 2.0 * k2_phi + 2.0 * k3_phi + k4_phi) / 6.0;
//		phidot += dt * (k1_phidot + 2.0 * k2_phidot + 2.0 * k3_phidot + k4_phidot) / 6.0;
//
//		t += dt;
//	}
//
//	out.close();
//	std::cout << "Free roll decay finished. Output: " << outFile << std::endl;
//}
