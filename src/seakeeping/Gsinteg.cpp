#include "Gsinteg.h"
#include "../const/Const.h"
#include <iostream>
#include "GreenTable.h"

namespace
{
	void getGaussLegendre1D(
		int order,
		std::vector<double>& xi,
		std::vector<double>& w)
	{
		xi.clear();
		w.clear();

		switch (order)
		{
		case 1:
			xi = { 0.0 };
			w = { 2.0 };
			break;

		case 2:
			xi = {
				-1.0 / std::sqrt(3.0),
				 1.0 / std::sqrt(3.0)
			};
			w = { 1.0, 1.0 };
			break;

		case 3:
			xi = {
				-std::sqrt(3.0 / 5.0),
				 0.0,
				 std::sqrt(3.0 / 5.0)
			};
			w = {
				5.0 / 9.0,
				8.0 / 9.0,
				5.0 / 9.0
			};
			break;

		case 4:
		{
			const double a = std::sqrt((3.0 - 2.0 * std::sqrt(6.0 / 5.0)) / 7.0);
			const double b = std::sqrt((3.0 + 2.0 * std::sqrt(6.0 / 5.0)) / 7.0);
			const double w1 = (18.0 + std::sqrt(30.0)) / 36.0;
			const double w2 = (18.0 - std::sqrt(30.0)) / 36.0;

			xi = { -b, -a, a, b };
			w = { w2, w1, w1, w2 };
			break;
		}

		case 16:
		{
			xi = {
				-0.989400934991649932596154173450,
				-0.944575023073232576077988415535,
				-0.865631202387831743880467897712,
				-0.755404408355003033895101194847,
				-0.617876244402643748446671764049,
				-0.458016777657227386342419442984,
				-0.281603550779258913230460501460,
				-0.095012509837637440185319335425,
				 0.095012509837637440185319335425,
				 0.281603550779258913230460501460,
				 0.458016777657227386342419442984,
				 0.617876244402643748446671764049,
				 0.755404408355003033895101194847,
				 0.865631202387831743880467897712,
				 0.944575023073232576077988415535,
				 0.989400934991649932596154173450
			};

			w = {
				0.027152459411754094851780572456,
				0.062253523938647892862843836994,
				0.095158511682492784809925107602,
				0.124628971255533872052476282192,
				0.149595988816576732081501730547,
				0.169156519395002538189312079030,
				0.182603415044923588866763667969,
				0.189450610455068496285396723208,
				0.189450610455068496285396723208,
				0.182603415044923588866763667969,
				0.169156519395002538189312079030,
				0.149595988816576732081501730547,
				0.124628971255533872052476282192,
				0.095158511682492784809925107602,
				0.062253523938647892862843836994,
				0.027152459411754094851780572456
			};

			break;
		}

		default:
			throw std::runtime_error("GreenCalPanelGauss: unsupported Gauss order.");
		}
	}

	inline void evalKernelAtPoint(
		double xField, double yField, double zField,
		double xSrc, double ySrc, double zSrc,
		double tn, double U,
		Greenf& GF,
		double& sG, double& xdG, double& ydG, double& zdG, double& tdG)
	{
		const double tx = U * tn + xField - xSrc;
		const double ty = yField - ySrc;
		const double tz = zSrc + zField;

		const double r0 = std::max(1.0e-14, std::sqrt(tx * tx + ty * ty + tz * tz));
		const double r02 = r0 * r0;

		const double Gp = std::sqrt(G / (r0 * r0 * r0));
		const double b = std::sqrt(G / r0) * tn;
		const double m = -tz / r0;

		GF.GreenFunctionCal(b, m);

		double dbx, dby, dbz;
		double dmx, dmy, dmz;
		double dGpx, dGpy, dGpz;

		{
			const double mlt = -0.5 * b / r02;
			dbx = mlt * tx;
			dby = mlt * ty;
			dbz = mlt * tz;
		}

		{
			const double mlt = -m / r02;
			dmx = mlt * tx;
			dmy = mlt * ty;
			dmz = (m * m - 1.0) / r0;
		}

		{
			const double mlt = -1.5 * Gp / r02;
			dGpx = mlt * tx;
			dGpy = mlt * ty;
			dGpz = mlt * tz;
		}

		xdG = Gp * (GF.Gmd * dmx + GF.Gbd * dbx) + GF.Gf * dGpx;
		ydG = Gp * (GF.Gmd * dmy + GF.Gbd * dby) + GF.Gf * dGpy;
		zdG = Gp * (GF.Gmd * dmz + GF.Gbd * dbz) + GF.Gf * dGpz;
		sG = GF.Gf * Gp;

		// 这里按“点值”给出，再由外层乘 Jacobian 和权重积分
		tdG = U * xdG + G * GF.Gbd / (r0 * r0);
	}

	inline void evalKernelAtPointTDGF(
		double xField, double yField, double zField,
		double xSrc, double ySrc, double zSrc,
		double tn, double U,
		const TDGFProvider& tdgf,
		double& sG, double& xdG, double& ydG, double& zdG, double& tdG)
	{
		const double tx = U * tn + xField - xSrc;
		const double ty = yField - ySrc;
		const double tz = zSrc + zField;

		const double r0 = std::max(1.0e-14, std::sqrt(tx * tx + ty * ty + tz * tz));
		const double r02 = r0 * r0;

		const double Gp = std::sqrt(G / (r0 * r0 * r0));
		const double b = std::sqrt(G / r0) * tn;
		const double m = -tz / r0;

		TDGFValue gf = tdgf.eval(b, m);

		double dbx, dby, dbz;
		double dmx, dmy, dmz;
		double dGpx, dGpy, dGpz;

		{
			const double mlt = -0.5 * b / r02;
			dbx = mlt * tx;
			dby = mlt * ty;
			dbz = mlt * tz;
		}

		{
			const double mlt = -m / r02;
			dmx = mlt * tx;
			dmy = mlt * ty;
			dmz = (m * m - 1.0) / r0;
		}

		{
			const double mlt = -1.5 * Gp / r02;
			dGpx = mlt * tx;
			dGpy = mlt * ty;
			dGpz = mlt * tz;
		}

		xdG = Gp * (gf.Fm * dmx + gf.Ft * dbx) + gf.F * dGpx;
		ydG = Gp * (gf.Fm * dmy + gf.Ft * dby) + gf.F * dGpy;
		zdG = Gp * (gf.Fm * dmz + gf.Ft * dbz) + gf.F * dGpz;

		sG = gf.F * Gp;
		tdG = U * xdG + G * gf.Ft / r02;
	}

	constexpr int GAUSS12_N = 12;

	// 12-point Gauss-Legendre nodes on [-1, 1]
	constexpr double GAUSS12_X[GAUSS12_N] =
	{
		-0.981560634246719250690549090149,
		-0.904117256370474856678465866119,
		-0.769902674194304687036893833213,
		-0.587317954286617447296702418941,
		-0.367831498998180193752691536644,
		-0.125233408511468915472441369464,
		 0.125233408511468915472441369464,
		 0.367831498998180193752691536644,
		 0.587317954286617447296702418941,
		 0.769902674194304687036893833213,
		 0.904117256370474856678465866119,
		 0.981560634246719250690549090149
	};

	// 12-point Gauss-Legendre weights on [-1, 1]
	constexpr double GAUSS12_W[GAUSS12_N] =
	{
		0.047175336386511827194615961485,
		0.106939325995318430960254718194,
		0.160078328543346226334652529543,
		0.203167426723065921749064455810,
		0.233492536538354808760849898925,
		0.249147045813402785000562436043,
		0.249147045813402785000562436043,
		0.233492536538354808760849898925,
		0.203167426723065921749064455810,
		0.160078328543346226334652529543,
		0.106939325995318430960254718194,
		0.047175336386511827194615961485
	};


	constexpr int GAUSS16_N = 16;

	// 16-point Gauss-Legendre nodes on [-1, 1]
	constexpr double GAUSS16_X[GAUSS16_N] =
	{
		-0.989400934991649932596154173450,
		-0.944575023073232576077988415535,
		-0.865631202387831743880467897712,
		-0.755404408355003033895101194847,
		-0.617876244402643748446671764049,
		-0.458016777657227386342419442984,
		-0.281603550779258913230460501460,
		-0.095012509837637440185319335425,
		 0.095012509837637440185319335425,
		 0.281603550779258913230460501460,
		 0.458016777657227386342419442984,
		 0.617876244402643748446671764049,
		 0.755404408355003033895101194847,
		 0.865631202387831743880467897712,
		 0.944575023073232576077988415535,
		 0.989400934991649932596154173450
	};

	// 16-point Gauss-Legendre weights on [-1, 1]
	constexpr double GAUSS16_W[GAUSS16_N] =
	{
		0.027152459411754094851780572456,
		0.062253523938647892862843836994,
		0.095158511682492784809925107602,
		0.124628971255533872052476282192,
		0.149595988816576732081501730547,
		0.169156519395002538189312079030,
		0.182603415044923588866763667969,
		0.189450610455068496285396723208,
		0.189450610455068496285396723208,
		0.182603415044923588866763667969,
		0.169156519395002538189312079030,
		0.149595988816576732081501730547,
		0.124628971255533872052476282192,
		0.095158511682492784809925107602,
		0.062253523938647892862843836994,
		0.027152459411754094851780572456
	};
}

inline void Gsinteg::GpMultiply(double& k, double& m, double& p)
{
	m = -1.5;
	m *= p;
	m /= k;
}

//void Gsinteg::Gderivative(double& dGpx, double& dbx, double& dmx, double& Gp, double& Gxr)
//{
//	double tx, ty;
//
//	//tx=dG[0];
//	tx = GF.Gf;
//	tx *= dGpx;
//	//ty=dG[1];
//	ty = GF.Gbd;
//	ty *= dbx;
//	Gxr = dmx;
//	Gxr *= GF.Gmd;//dG[2];
//	Gxr += ty;
//	Gxr *= Gp;
//	Gxr += tx;
//
//}



Gsinteg::Gsinteg(const std::shared_ptr<Element> element_s, double u)
:element(element_s), U(u), sG(0), xdG(0), ydG(0), zdG(0), tdG(0){}

Gsinteg::Gsinteg(const std::vector<ElementMatrix>& element, const std::vector<Vector3d>& point_test, double u)
	: element_test(element), point_test(point_test), U(u), sG(0), xdG(0), ydG(0), zdG(0), tdG(0) {
}

GreenData Gsinteg::GreenCal(const int j, const int i, const double tn, Greenf& GF, GreenTable& gGreenTable)
{
	double tx = U * tn + element->xnr(j) - element->xcr(i);
	double ty = element->ynr(j) - element->ycr(i);
	double tz = element->zcr(i) + element->znr(j);
	double r0 = sqrt(tx * tx + ty * ty + tz * tz);	
	double r02 = r0 * r0;
	double mlt, dbx, dby, dbz, dmx, dmy, dmz, dGpx, dGpy, dGpz;

	double Gp = sqrt(G / (r0 * r0 * r0));
	double b = sqrt(G / r0) * tn;
	double m = -tz / r0;

	//std::cout << "b:\t" << b << "\tm:\t" << m << std::endl;

	GF.GreenFunctionCal(b, m);
	//double gf, gbd, gmd;
	//gGreenTable.eval(b, m, gf, gbd, gmd);

	////   С Ķ   д   GF
	//GF.Gf = gf;
	//GF.Gbd = gbd;
	//GF.Gmd = gmd;
	////GF.GreenFunctionCal(b, m);


	//beta  x,y,z  ƫ  
	mlt = -0.5 * b / r02;
	dbx = mlt * tx;
	dby = mlt * ty;
	dbz = mlt * tz;

	//mu  x,y,z  ƫ  
	mlt = -m / r02;
	dmx = mlt * tx;
	dmy = mlt * ty;
	dmz = (m * m - 1.0) / r0;

	mlt = -1.5 * Gp / r02;
	dGpx = mlt * tx;
	dGpy = mlt * ty;
	dGpz = mlt * tz;

	xdG = Gp * (GF.Gmd * dmx + GF.Gbd * dbx) + GF.Gf * dGpx;
	ydG = Gp * (GF.Gmd * dmy + GF.Gbd * dby) + GF.Gf * dGpy;
	zdG = Gp * (GF.Gmd * dmz + GF.Gbd * dbz) + GF.Gf * dGpz;

	//sG:    ֺ      Ե i    Ԫ              
	double ar = element->Area[i];						
	xdG *= ar;		
	ydG *= ar;			
	zdG *= ar;
	sG = GF.Gf * Gp * ar;

	//   ֺ     ʱ   ƫ      tdG，
	tdG = U * xdG + G * GF.Gbd / (r0 * r0) * ar;

	return GreenData{ sG, xdG, ydG, zdG, tdG };
}


GreenData Gsinteg::GreenCalPanelGauss(
	const int j, const int i, const double tn,
	Greenf& GF, GreenTable& gGreenTable,
	int order)
{
	(void)gGreenTable; 

	std::vector<double> gp, gw;
	getGaussLegendre1D(order, gp, gw);

	double sG_sum = 0.0;
	double xdG_sum = 0.0;
	double ydG_sum = 0.0;
	double zdG_sum = 0.0;
	double tdG_sum = 0.0;

	// 场点/控制点：保持与旧 GreenCal 一致
	const double xField = element->xnr(j);
	const double yField = element->ynr(j);
	const double zField = element->znr(j);

	// 面元四角点：ElementData 是 4x3
	const ElementMatrix& P = element->ElementData->at(i);

	for (std::size_t a = 0; a < gp.size(); ++a)
	{
		const double xi = gp[a];
		const double wa = gw[a];

		for (std::size_t b = 0; b < gp.size(); ++b)
		{
			const double eta = gp[b];
			const double wb = gw[b];

			// 双线性形函数
			const double N1 = 0.25 * (1.0 - xi) * (1.0 - eta);
			const double N2 = 0.25 * (1.0 + xi) * (1.0 - eta);
			const double N3 = 0.25 * (1.0 + xi) * (1.0 + eta);
			const double N4 = 0.25 * (1.0 - xi) * (1.0 + eta);

			Eigen::Vector3d x =
				N1 * P.row(0).transpose() +
				N2 * P.row(1).transpose() +
				N3 * P.row(2).transpose() +
				N4 * P.row(3).transpose();

			const double dN1_dxi = -0.25 * (1.0 - eta);
			const double dN2_dxi = 0.25 * (1.0 - eta);
			const double dN3_dxi = 0.25 * (1.0 + eta);
			const double dN4_dxi = -0.25 * (1.0 + eta);

			const double dN1_deta = -0.25 * (1.0 - xi);
			const double dN2_deta = -0.25 * (1.0 + xi);
			const double dN3_deta = 0.25 * (1.0 + xi);
			const double dN4_deta = 0.25 * (1.0 - xi);

			Eigen::Vector3d dx_dxi =
				dN1_dxi * P.row(0).transpose() +
				dN2_dxi * P.row(1).transpose() +
				dN3_dxi * P.row(2).transpose() +
				dN4_dxi * P.row(3).transpose();

			Eigen::Vector3d dx_deta =
				dN1_deta * P.row(0).transpose() +
				dN2_deta * P.row(1).transpose() +
				dN3_deta * P.row(2).transpose() +
				dN4_deta * P.row(3).transpose();

			const double J = dx_dxi.cross(dx_deta).norm();

			double sG_pt, xdG_pt, ydG_pt, zdG_pt, tdG_pt;
			evalKernelAtPoint(
				xField, yField, zField,
				x(0), x(1), x(2),
				tn, U, GF,
				sG_pt, xdG_pt, ydG_pt, zdG_pt, tdG_pt);

			const double w2 = wa * wb * J;

			sG_sum += sG_pt * w2;
			xdG_sum += xdG_pt * w2;
			ydG_sum += ydG_pt * w2;
			zdG_sum += zdG_pt * w2;
			tdG_sum += tdG_pt * w2;
		}
	}

	return GreenData{ sG_sum, xdG_sum, ydG_sum, zdG_sum, tdG_sum };
}


std::vector<GreenData> Gsinteg::GreenCalPanelGauss_test(
	const double tn, Greenf& GF, int order)
{
	std::vector<double> gp, gw;
	getGaussLegendre1D(order, gp, gw);

	ElementMatrix P;
	Vector3d single_point;

	auto center = [](ElementMatrix P)
	{
		return 0.25 * (P.row(0) + P.row(1) + P.row(2) + P.row(3));
	};

	int n1 = element_test.size();
	int n2 = point_test.size();

	std::vector<GreenData> results;

	for (int i = 0; i < n1; ++i)
	{
		P = element_test[i];
		for(int j= 0; j < n2; ++j)
		{
			single_point = point_test[j];

			double sG_sum = 0.0;
			double xdG_sum = 0.0;
			double ydG_sum = 0.0;
			double zdG_sum = 0.0;
			double tdG_sum = 0.0;

			auto center_p = center(P);
			const double xField = center_p(0);
			const double yField = center_p(1);
			const double zField = center_p(2);

			// 场点/控制点：保持与旧 GreenCal 一致
			//const double xField = element->xnr(j);
			//const double yField = element->ynr(j);
			//const double zField = element->znr(j);

			// 面元四角点：ElementData 是 4x3
			//const ElementMatrix& P = element->ElementData->at(i);

			for (std::size_t a = 0; a < gp.size(); ++a)
			{
				const double xi = gp[a];
				const double wa = gw[a];

				for (std::size_t b = 0; b < gp.size(); ++b)
				{
					const double eta = gp[b];
					const double wb = gw[b];

					// 双线性形函数
					const double N1 = 0.25 * (1.0 - xi) * (1.0 - eta);
					const double N2 = 0.25 * (1.0 + xi) * (1.0 - eta);
					const double N3 = 0.25 * (1.0 + xi) * (1.0 + eta);
					const double N4 = 0.25 * (1.0 - xi) * (1.0 + eta);

					Eigen::Vector3d x =
						N1 * P.row(0).transpose() +
						N2 * P.row(1).transpose() +
						N3 * P.row(2).transpose() +
						N4 * P.row(3).transpose();

					const double dN1_dxi = -0.25 * (1.0 - eta);
					const double dN2_dxi = 0.25 * (1.0 - eta);
					const double dN3_dxi = 0.25 * (1.0 + eta);
					const double dN4_dxi = -0.25 * (1.0 + eta);

					const double dN1_deta = -0.25 * (1.0 - xi);
					const double dN2_deta = -0.25 * (1.0 + xi);
					const double dN3_deta = 0.25 * (1.0 + xi);
					const double dN4_deta = 0.25 * (1.0 - xi);

					Eigen::Vector3d dx_dxi =
						dN1_dxi * P.row(0).transpose() +
						dN2_dxi * P.row(1).transpose() +
						dN3_dxi * P.row(2).transpose() +
						dN4_dxi * P.row(3).transpose();

					Eigen::Vector3d dx_deta =
						dN1_deta * P.row(0).transpose() +
						dN2_deta * P.row(1).transpose() +
						dN3_deta * P.row(2).transpose() +
						dN4_deta * P.row(3).transpose();

					const double J = dx_dxi.cross(dx_deta).norm();

					double sG_pt, xdG_pt, ydG_pt, zdG_pt, tdG_pt;
					evalKernelAtPoint(
						xField, yField, zField,
						x(0), x(1), x(2),
						tn, U, GF,
						sG_pt, xdG_pt, ydG_pt, zdG_pt, tdG_pt);

					const double w2 = wa * wb * J;

					sG_sum += sG_pt * w2;
					xdG_sum += xdG_pt * w2;
					ydG_sum += ydG_pt * w2;
					zdG_sum += zdG_pt * w2;
					tdG_sum += tdG_pt * w2;
				}
			}
			results.push_back(GreenData{ sG_sum, xdG_sum, ydG_sum, zdG_sum, tdG_sum });
		}
	}
	return results;
}


GreenData Gsinteg::GreenCal_WL_dl(int& j, int& i, double& tn, Greenf& GF, GreenTable& gGreenTable)
{
	double Gp, tx, ty, tz;
	double mlt, b, m, r0, xps, xng, yps, yng;
	double r02, dx1, dy1, dz1;
	double dbx, dby, dbz, dmx, dmy, dmz, dGpx, dGpy, dGpz;

	xdG = ydG = zdG = sG = tdG = 0.0;

	const Eigen::MatrixXd& xpl = element->xpl;
	const Eigen::MatrixXd& ypl = element->ypl;
	const Eigen::VectorXd& xnr = element->xnr;
	const Eigen::VectorXd& ynr = element->ynr;
	const Eigen::VectorXd& znr = element->znr;

	tx = xpl(i, 0) - xpl(i, 1);
	ty = ypl(i, 0) - ypl(i, 1);

	double JJ = sqrt(pow(tx, 2) + pow(ty, 2)) / 2;
	//double JJ = 1.0;

	double xc[3] = { (xpl(i, 0) + xpl(i, 1)) / 2 - sqrt(0.6) * tx / 2,(xpl(i, 0) + xpl(i, 1)) / 2,
		(xpl(i, 0) + xpl(i, 1)) / 2 + sqrt(0.6) * tx / 2 };
	double yc[3] = { (ypl(i, 0) + ypl(i, 1)) / 2 - sqrt(0.6) * ty / 2,(ypl(i, 0) + ypl(i, 1)) / 2,
		(ypl(i, 0) + ypl(i, 1)) / 2 + sqrt(0.6) * ty / 2 };

	double w[3] = { 5.0 / 9.0,8.0 / 9.0,5.0 / 9.0 };

	for (int k = 0; k < 3; k++) {
		tx = U * tn + xnr[j] - xc[k];
		ty = ynr[j] - yc[k];
		tz = znr[j];
		r0 = sqrt(tx * tx + ty * ty + tz * tz);     
		Gp = sqrt(G / (r0 * r0 * r0));
		b  = sqrt(G / r0) * tn;
		m  = -tz / r0;

		GF.GreenFunctionCal(b, m);
		//double gf, gbd, gmd;
		//gGreenTable.eval(b, m, gf, gbd, gmd);

		////   С Ķ   д   GF
		//GF.Gf = gf;
		//GF.Gbd = gbd;
		//GF.Gmd = gmd;

		r02 = r0 * r0;
		mlt = -0.5 * b / r02;
		dbx = mlt * tx;
		dby = mlt * ty;
		dbz = mlt * tz;

		mlt = -m / r02;
		dmx = mlt * tx;
		dmy = mlt * ty;
		dmz = (m * m - 1.0) / r0;

		mlt = -1.5 * Gp / r02;
		dGpx = mlt * tx;
		dGpy = mlt * ty;
		dGpz = mlt * tz;

		dx1 = Gp * (GF.Gmd * dmx + GF.Gbd * dbx) + GF.Gf * dGpx;
		dy1 = Gp * (GF.Gmd * dmy + GF.Gbd * dby) + GF.Gf * dGpy;
		dz1 = Gp * (GF.Gmd * dmz + GF.Gbd * dbz) + GF.Gf * dGpz;

		mlt = JJ * w[k];
		xdG += dx1 * mlt;
		ydG += dy1 * mlt;
		zdG += dz1 * mlt;
		sG += GF.Gf * Gp * mlt;

		tdG += (U * dx1 + G * GF.Gbd / (r0 * r0)) * mlt;
	}
	//tdG = U * xdG + G * GF.Gbd / (r0 * r0);

	return GreenData{ sG, xdG, ydG, zdG, tdG };
}



GreenData Gsinteg::GreenCal_WL(int& j, int& i, double& tn, Greenf& GF, GreenTable& gGreenTable)
{
	double Gp, tx, ty, tz;
	double b, m, r0, r02;
	double dbx, dby, dbz, dmx, dmy, dmz, dGpx, dGpy, dGpz;

	xdG = ydG = zdG = sG = tdG = 0.0;

	const Eigen::MatrixXd& xpl = element->xpl;
	const Eigen::MatrixXd& ypl = element->ypl;
	const Eigen::VectorXd& xnr = element->xnr;
	const Eigen::VectorXd& ynr = element->ynr;
	const Eigen::VectorXd& znr = element->znr;

	// 线元代表点：中点
	const double xc = 0.5 * (xpl(i, 0) + xpl(i, 1));
	const double yc = 0.5 * (ypl(i, 0) + ypl(i, 1));

	// 场点 j 到线元中点的相对坐标
	tx = U * tn + xnr[j] - xc;
	ty = ynr[j] - yc;
	tz = znr[j];

	r0 = std::sqrt(tx * tx + ty * ty + tz * tz);
	r02 = r0 * r0;

	Gp = std::sqrt(G / (r0 * r0 * r0));
	b = std::sqrt(G / r0) * tn;
	m = -tz / r0;

	GF.GreenFunctionCal(b, m);
	// 如果你想用查表版，就改成下面三句：
	// double gf, gbd, gmd;
	// gGreenTable.eval(b, m, gf, gbd, gmd);
	// GF.Gf = gf; GF.Gbd = gbd; GF.Gmd = gmd;

	// beta 对 x,y,z 的偏导
	{
		const double mlt = -0.5 * b / r02;
		dbx = mlt * tx;
		dby = mlt * ty;
		dbz = mlt * tz;
	}

	// mu 对 x,y,z 的偏导
	{
		const double mlt = -m / r02;
		dmx = mlt * tx;
		dmy = mlt * ty;
		dmz = (m * m - 1.0) / r0;
	}

	// Gp 对 x,y,z 的偏导
	{
		const double mlt = -1.5 * Gp / r02;
		dGpx = mlt * tx;
		dGpy = mlt * ty;
		dGpz = mlt * tz;
	}

	// 点值导数（不再乘线长）
	xdG = Gp * (GF.Gmd * dmx + GF.Gbd * dbx) + GF.Gf * dGpx;
	ydG = Gp * (GF.Gmd * dmy + GF.Gbd * dby) + GF.Gf * dGpy;
	zdG = Gp * (GF.Gmd * dmz + GF.Gbd * dbz) + GF.Gf * dGpz;

	// 点值格林函数（不再乘线长）
	sG = GF.Gf * Gp;

	// 时间导数点值
	tdG = U * xdG + G * GF.Gbd / (r0 * r0);

	return GreenData{ sG, xdG, ydG, zdG, tdG };
}



GreenData Gsinteg::GreenCal_WL_Gauss12(int& j, int& i, double& tn, Greenf& GF, GreenTable& gGreenTable)
{
	double Gp, tx, ty, tz;
	double b, m, r0, r02;
	double dbx, dby, dbz, dmx, dmy, dmz, dGpx, dGpy, dGpz;

	xdG = ydG = zdG = sG = tdG = 0.0;

	const Eigen::MatrixXd& xpl = element->xpl;
	const Eigen::MatrixXd& ypl = element->ypl;
	const Eigen::VectorXd& xnr = element->xnr;
	const Eigen::VectorXd& ynr = element->ynr;
	const Eigen::VectorXd& znr = element->znr;

	// 当前水线线元两个端点
	const double x1 = xpl(i, 0);
	const double y1 = ypl(i, 0);
	const double x2 = xpl(i, 1);
	const double y2 = ypl(i, 1);

	const double dx = x2 - x1;
	const double dy = y2 - y1;

	// 线段长度的一半：dl = J dξ
	const double J = 0.5 * std::sqrt(dx * dx + dy * dy);

	// 防止退化线元
	if (J <= 1.0e-14)
	{
		return GreenData{ 0.0, 0.0, 0.0, 0.0, 0.0 };
	}

	for (int k = 0; k < GAUSS12_N; ++k)
	{
		const double xi = GAUSS12_X[k];
		const double wk = GAUSS12_W[k];

		// 从标准区间 [-1,1] 映射到实际线段
		// xc = 0.5*(1-xi)*x1 + 0.5*(1+xi)*x2
		// yc = 0.5*(1-xi)*y1 + 0.5*(1+xi)*y2
		const double xc = 0.5 * ((1.0 - xi) * x1 + (1.0 + xi) * x2);
		const double yc = 0.5 * ((1.0 - xi) * y1 + (1.0 + xi) * y2);

		// 场点 j 到当前高斯点的相对坐标
		tx = U * tn + xnr[j] - xc;
		ty = ynr[j] - yc;
		tz = znr[j];

		r02 = tx * tx + ty * ty + tz * tz;
		if (r02 < 1.0e-30) r02 = 1.0e-30;   // 防止除零
		r0 = std::sqrt(r02);

		Gp = std::sqrt(G / (r0 * r0 * r0));
		b = std::sqrt(G / r0) * tn;
		m = -tz / r0;

		// 直接函数计算版
		GF.GreenFunctionCal(b, m);

		// 如果你要切换成查表版，就把上面那句改成下面三句
		// double gf, gbd, gmd;
		// gGreenTable.eval(b, m, gf, gbd, gmd);
		// GF.Gf = gf; GF.Gbd = gbd; GF.Gmd = gmd;

		// beta 对 x,y,z 的偏导
		{
			const double mlt = -0.5 * b / r02;
			dbx = mlt * tx;
			dby = mlt * ty;
			dbz = mlt * tz;
		}

		// mu 对 x,y,z 的偏导
		{
			const double mlt = -m / r02;
			dmx = mlt * tx;
			dmy = mlt * ty;
			dmz = (m * m - 1.0) / r0;
		}

		// Gp 对 x,y,z 的偏导
		{
			const double mlt = -1.5 * Gp / r02;
			dGpx = mlt * tx;
			dGpy = mlt * ty;
			dGpz = mlt * tz;
		}

		// 当前高斯点处的点值
		const double dx1 = Gp * (GF.Gmd * dmx + GF.Gbd * dbx) + GF.Gf * dGpx;
		const double dy1 = Gp * (GF.Gmd * dmy + GF.Gbd * dby) + GF.Gf * dGpy;
		const double dz1 = Gp * (GF.Gmd * dmz + GF.Gbd * dbz) + GF.Gf * dGpz;
		const double sg1 = GF.Gf * Gp;

		// tdG 也要在每个高斯点上形成 integrand，再一起积分
		const double td1 = U * dx1 + G * GF.Gbd / r02;

		const double fac = J * wk;

		xdG += dx1 * fac;
		ydG += dy1 * fac;
		zdG += dz1 * fac;
		sG += sg1 * fac;
		tdG += td1 * fac;
	}

	return GreenData{ sG, xdG, ydG, zdG, tdG };
}

GreenData Gsinteg::GreenCal(
	const int j,
	const int i,
	const double tn,
	const TDGFProvider& tdgf)
{
	const double tx = U * tn + element->xnr(j) - element->xcr(i);
	const double ty = element->ynr(j) - element->ycr(i);
	const double tz = element->zcr(i) + element->znr(j);

	const double r0 = std::max(1.0e-14, std::sqrt(tx * tx + ty * ty + tz * tz));
	const double r02 = r0 * r0;

	const double Gp = std::sqrt(G / (r0 * r0 * r0));
	const double b = std::sqrt(G / r0) * tn;
	const double m = -tz / r0;

	TDGFValue gf = tdgf.eval(b, m);

	double dbx, dby, dbz;
	double dmx, dmy, dmz;
	double dGpx, dGpy, dGpz;

	{
		const double mlt = -0.5 * b / r02;
		dbx = mlt * tx;
		dby = mlt * ty;
		dbz = mlt * tz;
	}

	{
		const double mlt = -m / r02;
		dmx = mlt * tx;
		dmy = mlt * ty;
		dmz = (m * m - 1.0) / r0;
	}

	{
		const double mlt = -1.5 * Gp / r02;
		dGpx = mlt * tx;
		dGpy = mlt * ty;
		dGpz = mlt * tz;
	}

	xdG = Gp * (gf.Fm * dmx + gf.Ft * dbx) + gf.F * dGpx;
	ydG = Gp * (gf.Fm * dmy + gf.Ft * dby) + gf.F * dGpy;
	zdG = Gp * (gf.Fm * dmz + gf.Ft * dbz) + gf.F * dGpz;

	const double ar = element->Area[i];

	xdG *= ar;
	ydG *= ar;
	zdG *= ar;

	sG = gf.F * Gp * ar;
	tdG = U * xdG + G * gf.Ft / r02 * ar;

	return GreenData{ sG, xdG, ydG, zdG, tdG };
}


GreenData Gsinteg::GreenCalPanelGauss(
	const int j,
	const int i,
	const double tn,
	const TDGFProvider& tdgf,
	int order)
{
	std::vector<double> gp, gw;
	getGaussLegendre1D(order, gp, gw);

	double sG_sum = 0.0;
	double xdG_sum = 0.0;
	double ydG_sum = 0.0;
	double zdG_sum = 0.0;
	double tdG_sum = 0.0;

	const double xField = element->xnr(j);
	const double yField = element->ynr(j);
	const double zField = element->znr(j);

	const ElementMatrix& P = element->ElementData->at(i);

	for (std::size_t a = 0; a < gp.size(); ++a)
	{
		const double xi = gp[a];
		const double wa = gw[a];

		for (std::size_t b = 0; b < gp.size(); ++b)
		{
			const double eta = gp[b];
			const double wb = gw[b];

			const double N1 = 0.25 * (1.0 - xi) * (1.0 - eta);
			const double N2 = 0.25 * (1.0 + xi) * (1.0 - eta);
			const double N3 = 0.25 * (1.0 + xi) * (1.0 + eta);
			const double N4 = 0.25 * (1.0 - xi) * (1.0 + eta);

			Eigen::Vector3d x =
				N1 * P.row(0).transpose() +
				N2 * P.row(1).transpose() +
				N3 * P.row(2).transpose() +
				N4 * P.row(3).transpose();

			const double dN1_dxi = -0.25 * (1.0 - eta);
			const double dN2_dxi = 0.25 * (1.0 - eta);
			const double dN3_dxi = 0.25 * (1.0 + eta);
			const double dN4_dxi = -0.25 * (1.0 + eta);

			const double dN1_deta = -0.25 * (1.0 - xi);
			const double dN2_deta = -0.25 * (1.0 + xi);
			const double dN3_deta = 0.25 * (1.0 + xi);
			const double dN4_deta = 0.25 * (1.0 - xi);

			Eigen::Vector3d dx_dxi =
				dN1_dxi * P.row(0).transpose() +
				dN2_dxi * P.row(1).transpose() +
				dN3_dxi * P.row(2).transpose() +
				dN4_dxi * P.row(3).transpose();

			Eigen::Vector3d dx_deta =
				dN1_deta * P.row(0).transpose() +
				dN2_deta * P.row(1).transpose() +
				dN3_deta * P.row(2).transpose() +
				dN4_deta * P.row(3).transpose();

			const double J = dx_dxi.cross(dx_deta).norm();

			double sG_pt, xdG_pt, ydG_pt, zdG_pt, tdG_pt;

			evalKernelAtPointTDGF(
				xField, yField, zField,
				x(0), x(1), x(2),
				tn, U, tdgf,
				sG_pt, xdG_pt, ydG_pt, zdG_pt, tdG_pt);

			const double w2 = wa * wb * J;

			sG_sum += sG_pt * w2;
			xdG_sum += xdG_pt * w2;
			ydG_sum += ydG_pt * w2;
			zdG_sum += zdG_pt * w2;
			tdG_sum += tdG_pt * w2;
		}
	}

	return GreenData{ sG_sum, xdG_sum, ydG_sum, zdG_sum, tdG_sum };
}



GreenData Gsinteg::GreenCal_WL_Gauss12(
	int& j,
	int& i,
	double& tn,
	const TDGFProvider& tdgf)
{
	xdG = ydG = zdG = sG = tdG = 0.0;

	const Eigen::MatrixXd& xpl = element->xpl;
	const Eigen::MatrixXd& ypl = element->ypl;
	const Eigen::VectorXd& xnr = element->xnr;
	const Eigen::VectorXd& ynr = element->ynr;
	const Eigen::VectorXd& znr = element->znr;

	const double x1 = xpl(i, 0);
	const double y1 = ypl(i, 0);
	const double x2 = xpl(i, 1);
	const double y2 = ypl(i, 1);

	const double dx = x2 - x1;
	const double dy = y2 - y1;

	const double J = 0.5 * std::sqrt(dx * dx + dy * dy);

	if (J <= 1.0e-14)
	{
		return GreenData{ 0.0, 0.0, 0.0, 0.0, 0.0 };
	}

	for (int k = 0; k < GAUSS12_N; ++k)
	{
		const double xi = GAUSS12_X[k];
		const double wk = GAUSS12_W[k];

		const double xc = 0.5 * ((1.0 - xi) * x1 + (1.0 + xi) * x2);
		const double yc = 0.5 * ((1.0 - xi) * y1 + (1.0 + xi) * y2);

		const double tx = U * tn + xnr[j] - xc;
		const double ty = ynr[j] - yc;
		const double tz = znr[j];

		double r02 = tx * tx + ty * ty + tz * tz;
		if (r02 < 1.0e-30)
			r02 = 1.0e-30;

		const double r0 = std::sqrt(r02);

		const double Gp = std::sqrt(G / (r0 * r0 * r0));
		const double b = std::sqrt(G / r0) * tn;
		const double m = -tz / r0;

		TDGFValue gf = tdgf.eval(b, m);

		double dbx, dby, dbz;
		double dmx, dmy, dmz;
		double dGpx, dGpy, dGpz;

		{
			const double mlt = -0.5 * b / r02;
			dbx = mlt * tx;
			dby = mlt * ty;
			dbz = mlt * tz;
		}

		{
			const double mlt = -m / r02;
			dmx = mlt * tx;
			dmy = mlt * ty;
			dmz = (m * m - 1.0) / r0;
		}

		{
			const double mlt = -1.5 * Gp / r02;
			dGpx = mlt * tx;
			dGpy = mlt * ty;
			dGpz = mlt * tz;
		}

		const double dx1 = Gp * (gf.Fm * dmx + gf.Ft * dbx) + gf.F * dGpx;
		const double dy1 = Gp * (gf.Fm * dmy + gf.Ft * dby) + gf.F * dGpy;
		const double dz1 = Gp * (gf.Fm * dmz + gf.Ft * dbz) + gf.F * dGpz;
		const double sg1 = gf.F * Gp;
		const double td1 = U * dx1 + G * gf.Ft / r02;

		const double fac = J * wk;

		xdG += dx1 * fac;
		ydG += dy1 * fac;
		zdG += dz1 * fac;
		sG += sg1 * fac;
		tdG += td1 * fac;
	}

	return GreenData{ sG, xdG, ydG, zdG, tdG };
}


GreenData Gsinteg::GreenCal_WL_Gauss16(
	int& j,
	int& i,
	double& tn,
	const TDGFProvider& tdgf)
{
	xdG = ydG = zdG = sG = tdG = 0.0;

	const Eigen::MatrixXd& xpl = element->xpl;
	const Eigen::MatrixXd& ypl = element->ypl;
	const Eigen::VectorXd& xnr = element->xnr;
	const Eigen::VectorXd& ynr = element->ynr;
	const Eigen::VectorXd& znr = element->znr;

	const double x1 = xpl(i, 0);
	const double y1 = ypl(i, 0);
	const double x2 = xpl(i, 1);
	const double y2 = ypl(i, 1);

	const double dx = x2 - x1;
	const double dy = y2 - y1;

	const double J = 0.5 * std::sqrt(dx * dx + dy * dy);

	if (J <= 1.0e-14)
	{
		return GreenData{ 0.0, 0.0, 0.0, 0.0, 0.0 };
	}

	for (int k = 0; k < GAUSS16_N; ++k)
	{
		const double xi = GAUSS16_X[k];
		const double wk = GAUSS16_W[k];

		const double xc = 0.5 * ((1.0 - xi) * x1 + (1.0 + xi) * x2);
		const double yc = 0.5 * ((1.0 - xi) * y1 + (1.0 + xi) * y2);

		const double tx = U * tn + xnr[j] - xc;
		const double ty = ynr[j] - yc;
		const double tz = znr[j];

		double r02 = tx * tx + ty * ty + tz * tz;
		if (r02 < 1.0e-30)
			r02 = 1.0e-30;

		const double r0 = std::sqrt(r02);
		const double Gp = std::sqrt(G / (r0 * r02));
		const double b = std::sqrt(G / r0) * tn;
		const double m = -tz / r0;

		TDGFValue gf = tdgf.eval(b, m);

		double dbx, dby, dbz;
		double dmx, dmy, dmz;
		double dGpx, dGpy, dGpz;

		{
			const double mlt = -0.5 * b / r02;
			dbx = mlt * tx;
			dby = mlt * ty;
			dbz = mlt * tz;
		}

		{
			const double mlt = -m / r02;
			dmx = mlt * tx;
			dmy = mlt * ty;
			dmz = (m * m - 1.0) / r0;
		}

		{
			const double mlt = -1.5 * Gp / r02;
			dGpx = mlt * tx;
			dGpy = mlt * ty;
			dGpz = mlt * tz;
		}

		const double xdG_pt =
			Gp * (gf.Fm * dmx + gf.Ft * dbx) + gf.F * dGpx;

		const double ydG_pt =
			Gp * (gf.Fm * dmy + gf.Ft * dby) + gf.F * dGpy;

		const double zdG_pt =
			Gp * (gf.Fm * dmz + gf.Ft * dbz) + gf.F * dGpz;

		const double sG_pt = gf.F * Gp;

		const double tdG_pt =
			U * xdG_pt + G * gf.Ft / r02;

		const double w1 = wk * J;

		sG += sG_pt * w1;
		xdG += xdG_pt * w1;
		ydG += ydG_pt * w1;
		zdG += zdG_pt * w1;
		tdG += tdG_pt * w1;
	}

	return GreenData{ sG, xdG, ydG, zdG, tdG };
}