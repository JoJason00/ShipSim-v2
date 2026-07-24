#include "WaveBase.h"
#include "RegularWave.h"
#include "../const/Const.h"
#include <iostream>
#include <complex>
#include <cmath>

using ComplexD = std::complex<double>;

RegularWave::RegularWave(const RegularWaveConfig& regularwave)
:config(regularwave)
{
	wt[0][0] = 0.342901327223705;     wt[0][1] = 6.108626337353e-1;
	wt[1][0] = 1.036610829789514;     wt[1][1] = 2.401386110823e-1;
	wt[2][0] = 1.756683649299882;     wt[2][1] = 3.387439445548e-2;
	wt[3][0] = 2.532731674232790;     wt[3][1] = 1.343645746781e-3;
	wt[4][0] = 3.436159118837738;     wt[4][1] = 7.640432855233e-6;
	wt[5][0] = -0.342901327223705;    wt[5][1] = 6.108626337353e-1;
	wt[6][0] = -1.036610829789514;    wt[6][1] = 2.401386110823e-1;
	wt[7][0] = -1.756683649299882;    wt[7][1] = 3.387439445548e-2;
	wt[8][0] = -2.532731674232790;    wt[8][1] = 1.343645746781e-3;
	wt[9][0] = -3.436159118837738;    wt[9][1] = 7.640432855233e-6;
}

void RegularWave::loadData(const fkpData& Data)
{
	//if (!data.has_value())
	//{
		data.emplace(Data);
		w = config.w;
		//we = w - w * w * data->U * cos(config.direction) / G;
		encounter_ = calcEncounterInfo(w, data->U, config.direction);
		we = encounter_.we;
	//}
	//else
	//{
	//	throw std::runtime_error("data already has value! (RegularWave.cpp)\n");
	//}
};

double RegularWave::Eta(double t) const
{
	return config.H * cos(we * t + config.phase0);
}

double RegularWave::initialPhase()
{
	return config.phase0;
}

void RegularWave::Exciting(double tn, FKphi& phi)
{
	double Vn, kz, fz;

	double eVx = 0;
	double eVy = 0;
	double eVz = 0;
	double ePr = 0;

	Eigen::VectorXd& Kpz = phi.df;
	Eigen::VectorXd& FKz = phi.fk;

	//0    б  
	int sw = (config.direction >= PI / 2 && config.direction <= (3.0 * PI / 2.0)) ? 0 : 1;

	for (int i = 0; i < data->NE; i++) {
		//       䡢           Ӧ    K=  eVx,eVy,eVz    p=ePr
		if (sw == 0)
			Headseas(tn, i, eVx, eVy, eVz, ePr);
		else     
			FollowingSeas(tn, i, eVx, eVy, eVz, ePr);

		//      Ԫ        A31  A32  A33  
		Vn = (data->A31[i] * eVx + data->A32[i] * eVy + data->A33[i] * eVz);
		kz = Vn;
		fz = ePr;
     
		//Kpz:fai7,FKz:p(fai0)    			
		Kpz[i] = kz;
		FKz[i] = fz;


		//   岨 ٶ  ƶ ʱ  ĵ        ֵ FK    
		//FKz[tN][i] = G*exp(W*W/G*(zcr[i]))*(cos(xcr[i]*cos(HeadingAngle)+ycr[i]*sin(HeadingAngle))*cos(W*tn)
	//+sin(xcr[i]*cos(HeadingAngle)+ycr[i]*sin(HeadingAngle))*sin(W*tn));

	}
}

double RegularWave::direction()
{
	return config.direction;
}

double RegularWave::getAmp()
{
	return config.H;
}

double RegularWave::getFreq()
{
	return config.w;
}

void RegularWave::Headseas(double tn, int i, double& eVx, double& eVy, double& eVz, double& ePr)
{
	double x, y, csb, snb;

	double beta = config.direction;
	std::array<double, 2> Riw;

	x = -data->zcr[i] / G;
	csb = cos(beta);
	snb = sin(beta);
	y = (data->xcr[i] * csb + data->ycr[i] * snb + data->U * tn * csb) / G;

	ComplexD a(x, y);  			//x+iy
	ComplexD ap = std::sqrt(a);		//a  ƽ    
	ComplexD b(0, -tn / 2);
	ComplexD tp(0, 1);			//      λ
	ComplexD ba = tp * b / ap;			//    i*b/sqrt(a)
	x = std::real(ba);
	y = std::imag(ba);

	//Riw:      
	Cerrorfun(x, y, Riw);

	ComplexD wz(Riw[0], Riw[1]);
	ComplexD I0 = 0.5 * sqrt(PI) * wz / ap;
	ComplexD I1 = (-1.0 * b * I0 + 0.5) / a;
	ComplexD I2 = (-1.0 * b * I1 + 0.5 * I0) / a;

	//      Ӧ    K  p
	eVx = csb * (std::real(I1) - 2.0 * data->U * csb * std::real(I2) / G) / PI;
	eVy = snb * (std::real(I1) - 2.0 * data->U * csb * std::real(I2) / G) / PI;
	eVz = -(std::imag(I1) - 2.0 * data->U * csb * std::imag(I2) / G) / PI;
	ePr = G * (std::real(I0) - 2.0 * data->U * csb * std::real(I1) / G) / PI;
}


void RegularWave::FollowingSeas(double tn, int i, double& eVx, double& eVy, double& eVz, double& ePr)
{
	const double beta = config.direction;
	const double csb = std::cos(beta);
	const double snb = std::sin(beta);
	const double speedProjection = data->U * csb;

	// The closed-form following-seas expressions degenerate to the head-seas form as U cos(beta) -> 0.
	if (speedProjection <= 1e-12) {
		Headseas(tn, i, eVx, eVy, eVz, ePr);
		return;
	}

	const ComplexD a(
		-data->zcr[i] / G,
		(data->xcr[i] * csb + data->ycr[i] * snb + speedProjection * tn) / G);
	const ComplexD b(0.0, -tn / 2.0);
	const ComplexD imagUnit(0.0, 1.0);

	auto complexError = [this](const ComplexD& z) {
		double xr = std::real(z);
		double yi = std::imag(z);
		std::array<double, 2> Riw{};
		Cerrorfun(xr, yi, Riw);
		return ComplexD(Riw[0], Riw[1]);
	};

	auto intervalIntegrals = [&](double lower, double upper, bool upperInfinite) {
		const ComplexD rootA = std::sqrt(a);
		const auto tailIntegral = [&](double bound) {
			const ComplexD exponent = -(a * bound * bound + 2.0 * b * bound);
			const ComplexD z = imagUnit * (b + a * bound) / rootA;
			return 0.5 * std::sqrt(PI) * std::exp(exponent) * complexError(z) / rootA;
		};

		const ComplexD I0Lower = tailIntegral(lower);
		const ComplexD I0Upper = upperInfinite ? ComplexD(0.0, 0.0) : tailIntegral(upper);
		const ComplexD I0 = I0Lower - I0Upper;

		const ComplexD expLower = std::exp(-(a * lower * lower + 2.0 * b * lower));
		const ComplexD expUpper = upperInfinite
			? ComplexD(0.0, 0.0)
			: std::exp(-(a * upper * upper + 2.0 * b * upper));

		const ComplexD I1 = (-b * I0 - 0.5 * (expUpper - expLower)) / a;

		const ComplexD upperBoundary = upperInfinite ? ComplexD(0.0, 0.0) : upper * expUpper;
		const ComplexD lowerBoundary = lower * expLower;
		const ComplexD I2 = (-b * I1 + 0.5 * I0 - 0.5 * (upperBoundary - lowerBoundary)) / a;

		return std::array<ComplexD, 3>{ I0, I1, I2 };
	};

	const double split1 = encounter_.split1;
	const double split2 = encounter_.split2;

	std::array<ComplexD, 3> integrals;
	if (config.w <= split1) {
		integrals = intervalIntegrals(0.0, split1, false);
	}
	else if (config.w < split2) {
		integrals = intervalIntegrals(split1, split2, false);
	}
	else {
		integrals = intervalIntegrals(split2, 0.0, true);
	}

	const ComplexD& I0 = integrals[0];
	const ComplexD& I1 = integrals[1];
	const ComplexD& I2 = integrals[2];
	const double speedFactor = 2.0 * speedProjection / G;

	eVx = csb * std::real(I1 - speedFactor * I2) / PI;
	eVy = snb * std::real(I1 - speedFactor * I2) / PI;
	eVz = -std::imag(I1 - speedFactor * I2) / PI;
	ePr = G * std::real(I0 - speedFactor * I1) / PI;
}

void RegularWave::Cerrorfun(double& x, double& y, std::array<double, 2>& Riw)
{
	double sgnx, sgny;

	sgnx = x;
	sgny = y;
	x = fabs(x);
	y = fabs(y);
	if (x < 5.33 && y < 4.29)  Taylor(Riw, x, y);
	else					   Gaussh(Riw, x, y);
	// Լ        д   
	Quadrant(Riw, x, y, sgnx, sgny);
}

void RegularWave::Taylor(std::array<double, 2>& Riw, double& x, double& y)
{
	int N, vv, n;
	double sz, hh, rr, ri, sr, si, t1, t2, c1;

	t1 = x; t1 /= 5.33;
	t1 = pow(t1, 2.0); t1 *= -1.0;
	t1 += 1.0; t1 = sqrt(t1);
	sz = -1.0; sz *= y; sz /= 4.29;
	sz += 1.0; sz *= t1;
	//sz=(1.0-y/4.29)*sqrt(1.0-pow((x/5.33),2));
	hh = 1.6; hh *= sz;
	t1 = 23.0; t1 *= sz; t1 += 6.0;
	N = (int)t1;
	//n=(int)(6.0+23.0*sz);
	t1 = 21.0; t1 *= sz; t1 += 9.0;
	vv = (int)t1;
	//vv=(int)(9.0+21.0*sz);
	rr = ri = sr = si = 0.0;
	for (n = vv; n >= 0; n--) {
		t1 = rr; t1 *= (n + 1); t1 += hh; t1 += y;
		//t1=y+hh+(n+1)*rr;	
		t2 = -1.0; t2 *= ri; t2 *= (n + 1); t2 += x;
		//t2=x-(n+1)*ri;
		sz = pow(t1, 2.0);
		sz += pow(t2, 2.0);
		c1 = 0.5; c1 /= sz;
		//c1=0.5/(t1*t1+t2*t2);
		rr = c1; rr *= t1;
		ri = c1; ri *= t2;
		if (n <= N) {
			t1 = 2.0; t1 *= hh; t1 = pow(t1, n); t1 += sr;
			//t1=pow((2.0*hh),n)+sr;
			sr = rr; sr *= t1;
			t2 = ri; t2 *= si;
			sr -= t2;
			//sr=rr*t1-ri*si;
			si = si;
			si *= rr;
			t2 = ri;
			t2 *= t1;
			si += t2;
			//si=ri*t1+rr*si;
		}
	}
	t1 = 2.0;
	t1 /= sqrt(PI);
	Riw[0] = t1 * sr;
	//Riw[0]=2.0/sqrt(PI)*sr;
	Riw[1] = t1 * si;
	//Riw[1]=2.0/sqrt(PI)*si;	
}

void RegularWave::Gaussh(std::array<double, 2>& Riw, double& x, double& y)
{
	int s;
	double Rew, Imw, rn1, rn2, rn3, Re, Im;

	Rew = 0.0; Imw = 0.0;
	for (s = 0; s < 10; s++) {
		rn1 = y;
		rn2 = x;
		rn2 -= wt[s][0];
		rn2 = pow(rn2, 2.0);
		rn2 += pow(y, 2.0);
		//rn2=pow((x-wt[s][0]),2)+y*y;
		rn3 = rn1;
		rn3 /= PI;
		rn3 /= rn2;
		//rn3=rn1/(PI*rn2);				
		Re = wt[s][1];
		Re *= rn3;
		//Re=wt[s][1]*rn3;
		Rew += Re;
		//Rew=Re+Rew;
		rn1 = x;
		rn1 -= wt[s][0];
		//rn1=(x-wt[s][0]);
		rn3 = rn1;
		rn3 /= PI;
		rn3 /= rn2;
		//rn3=rn1/(PI*rn2);
		Im = wt[s][1];
		Im *= rn3;
		//Im=wt[s][1]*rn3;
		Imw += Im;
		//Imw=Im+Imw;
	}
	Riw[0] = Rew;
	Riw[1] = Imw;
}

void RegularWave::Quadrant(std::array<double, 2>& Riw, double& x, double& y, double& sgnx, double& sgny)
{
	double Rew, Imw, Rew1, Imw1;

	Rew = Riw[0];
	Imw = Riw[1];
	if (sgnx < 0 && sgny>0) {
		Rew = Rew;
		Imw = -Imw;
	}
	else if (sgnx < 0 && sgny < 0) {
		//Rew1=2.0*exp(-x*x+y*y)*cos(2.0*x*y);
		//Imw1=-2.0*exp(-x*x+y*y)*sin(2.0*x*y);		
		//Rew=Rew1-Rew;
		//Imw=Imw1-Imw;
		Quadrant1(-1.0, x, y, Rew1, Imw1, Rew, Imw);

	}
	else if (sgnx > 0 && sgny < 0) {
		//Rew1=2.0*exp(-x*x+y*y)*cos(2.0*x*y);
		//Imw1=2.0*exp(-x*x+y*y)*sin(2.0*x*y);		
		//Rew=Rew1-Rew;
		//Imw=Imw1+Imw;
		Quadrant1(1.0, x, y, Rew1, Imw1, Rew, Imw);
	}
	Riw[0] = Rew;
	Riw[1] = Imw;
}

void RegularWave::Quadrant1(double k, double& x, double& y, double& Rew1, double& Imw1, double& Rew, double& Imw)
{
	double m, n, t;
	m = 2.0;
	m *= x;
	m *= y;
	n = m;
	m = cos(m);
	n = sin(n);
	t = pow(y, 2.0);
	t -= pow(x, 2.0);
	t = exp(t);
	t *= 2.0;
	Rew1 = t;
	Rew1 *= m;
	//Rew1=2.0*exp(-x*x+y*y)*cos(2.0*x*y);
	Imw1 = t;
	Imw1 *= n;
	//Imw1=2.0*exp(-x*x+y*y)*sin(2.0*x*y);		
	t = Rew1;
	t -= Rew;
	Rew = t;
	//Rew=Rew1-Rew;
	t = Imw;
	t += Imw1;
	t *= k;
	Imw = t;
	//Imw=k*(Imw1+Imw);
}


