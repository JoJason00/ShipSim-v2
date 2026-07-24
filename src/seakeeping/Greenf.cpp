#include <math.h>
#include "greenf.h"
#include "../const/Const.h"
#include <iostream>
#include "GreenFunction/series.h"

const double EPS = 1.0e-5;

Greenf::Greenf()
	:Gf(0),Gbd(0),Gmd(0)//,cp1(4), cp2(6),
{
	//GreenfHummerPoints();
}

const std::array<double, 47> Greenf::factor = [] {
	std::array<double, 47> f{};
	f[0] = 1.0; f[1] = 1.0;
	for (int r = 2; r < 47; ++r) f[r] = f[r - 1] * r;
	return f;
	}();

const std::array<double, 23> Greenf::anh = [] {
	std::array<double, 23> a{};
	for (int n = 0; n < 23; ++n) {
		double an1 = 0.0;
		for (int r = 0; r <= n; ++r) {
			double an2 = Greenf::factor[n + r];
			an2 *= pow(5.0, -2.0 * r - 1.0);
			an2 *= (pow(-1.0, r + n) - exp(-25.0));
			an2 /= (Greenf::factor[n - r] * Greenf::factor[r]);
			an1 += an2;
		}
		a[n] = an1;
	}
	return a;
	}();

const std::array<double, 6> Greenf::an  = { 1.0, -2.2499997, 1.2656208, -0.3163866, 0.0444479, 0.00021 };
const std::vector<double> Greenf::cp1   = { 0.51242420, 0.2752551, 0.05176536, 2.7247450 };
const std::vector<double> Greenf::cp2   = { 0.461313500, 0.190163500, 0.099992160, 1.784492700, 0.002883894, 5.525343700 };


//void
//Greenf::GreenfHummerPoints()
//{
//	int n, r;
//	double an1, an2;
//
//	//    ׳ˣ factor[i]=i!
//	factor[1] = factor[0] = 1.0;
//	for (r = 2; r < 47; r++) factor[r] = factor[r - 1] * r;
//
//	for (n = 0; n < 23; n++) {
//		an1 = 0.0;
//		for (r = 0; r <= n; r++) {
//			an2 = factor[n + r];
//			an2 *= pow(5.0, (-2 * r - 1));
//			an2 *= (pow(-1, (r + n)) - exp(-25.0));
//			an2 /= (factor[n - r] * factor[r]);
//			an1 += an2;
//		}
//		anh[n] = an1;
//	}
//
//	an[0] = 1.0000000; an[1] = -2.2499997;
//	an[2] = 1.2656208; an[3] = -0.3163866;
//	an[4] = 0.0444479; an[5] = 0.0002100;
//
//	cp1[0] = 0.51242420; cp1[1] = 0.2752551;
//	cp1[2] = 0.05176536; cp1[3] = 2.7247450;
//
//	cp2[0] = 0.461313500; cp2[1] = 0.190163500;
//	cp2[2] = 0.099992160; cp2[3] = 1.784492700;
//	cp2[4] = 0.002883894; cp2[5] = 5.525343700;
//}

void Greenf::Series1(const double b, const double m)
{
	const int smax = 135;

	double Pn0  = 1.0,  Pn1 = m;
	double Pnd0 = 0.0,  Pnd1 = 1.0;

	double Gf1  = 2.0 * b * Pn1;
	double Gbd1 = 2.0 * Pn1;
	double Gmd1 = 2.0 * b * Pnd1;

	double fact    = 1.0;
	double sign    = 1.0;       // (-1)^(s-1)
	double btk     = b;          // b^(2*s-1)
	double bsquare = b * b;

	int tk, tp;
	double Pn, Pnd;

	for (int v = 0, r = 1, s = 2; s <= smax; ++v, ++r, ++s)
	{
		fact = fact * (s + v) * (s + r) / s;

		tk = 2 * s - 1;
		tp = s - 1;

		Pn  = (Pn1 * m * tk - tp * Pn0) / s;
		Pnd = (tk * m * Pnd1 + tk * Pn1 - tp * Pnd0) / s;

		//Pn=((2*s-1)*m*Pn1-(s-1)*Pn0)/s;
		//Pnd=((2*s-1)*m*Pnd1+(2*s-1)*Pn1-(s-1)*Pnd0)/s;

		Pn0  = Pn1;  Pn1 = Pn;
		Pnd0 = Pnd1; Pnd1 = Pnd;

		sign = -sign;
		btk *= bsquare;
	
		Gf  = 2.0 * sign / fact * btk;
		Gmd = Gf * Pnd;
		Gf *= Pn;
		Gbd = Gf * tk / b;					// Gbd = Gf * tk * b^(tk-1)

		//Gf=2.0*pow(b,(2*s-1))*pow(-1,(s-1))*Pn/fact;                  		
		//Gbd=2.0*(2.0*s-1)*pow(b,(2*s-2))*pow(-1,(s-1))*Pn/fact;				  
		//Gmd=2.0*pow(b,(2*s-1))*pow(-1,(s-1))*Pnd/fact; 

		Gf  += Gf1;
		Gbd += Gbd1;
		Gmd += Gmd1;

		if (fabs(Gf - Gf1) < EPS)
			if (fabs(Gbd - Gbd1) < EPS)
				if (fabs(Gmd - Gmd1) < EPS)
					break;

		Gf1 = Gf; Gbd1 = Gbd; Gmd1 = Gmd;
	}
}


void Greenf::Series2(const double b, const double m)
{
	const int smax = 135;

	double Pn0 = 1.0, Pn1 = m;
	double Pnd0 = 0.0, Pnd1 = 1.0;

	double Gf1 = 2.0 * b * Pn1;
	double Gbd1 = 2.0 * Pn1;
	double Gmd1 = 2.0 * b * Pnd1;

	double fact = 1.0;
	double sign = 1.0;       // (-1)^(s-1)
	double btk = b;          // b^(2*s-1)
	double bsquare = b * b;

	int tk, tp;
	double Pn, Pnd;

	for (int v = 0, r = 1, s = 2; s <= smax; ++v, ++r, ++s)
	{
		fact = fact * (s + v) * (s + r) / s;

		tk = 2 * s - 1;
		tp = s - 1;

		Pn = (Pn1 * m * tk - tp * Pn0) / s;
		Pnd = (tk * m * Pnd1 + tk * Pn1 - tp * Pnd0) / s;

		//Pn=((2*s-1)*m*Pn1-(s-1)*Pn0)/s;
		//Pnd=((2*s-1)*m*Pnd1+(2*s-1)*Pn1-(s-1)*Pnd0)/s;

		Pn0 = Pn1;  Pn1 = Pn;
		Pnd0 = Pnd1; Pnd1 = Pnd;

		sign = -sign;
		btk *= bsquare;

		Gf = 2.0 * sign / fact * btk;
		Gmd = Gf * Pnd;
		Gf *= Pn;
		Gbd = Gf * tk / b;					// Gbd = Gf * tk * b^(tk-1)

		//Gf=2.0*pow(b,(2*s-1))*pow(-1,(s-1))*Pn/fact;                  		
		//Gbd=2.0*(2.0*s-1)*pow(b,(2*s-2))*pow(-1,(s-1))*Pn/fact;				  
		//Gmd=2.0*pow(b,(2*s-1))*pow(-1,(s-1))*Pnd/fact; 

		Gf += Gf1;
		Gbd += Gbd1;
		Gmd += Gmd1;

		Gf1 = Gf; Gbd1 = Gbd; Gmd1 = Gmd;
	}
}

void
Greenf::Asymp1(const double m, const double b, double& b3, double& b4)
{
	int r, s;
	double fact, fact1, Pn, Pn0, Pn1, Pnd;
	double tp, tr, Pnd0, Pnd1, Gbd1, Gmd1, Gf1;

	fact1 = 1.0; Pn0 = 1.0; Pn1 = m;
	Pnd0  = 0.0;  Pnd1 = 1.0;
	Gf1   = -8.0;  Gf1 /= b3;
	Gbd1  = 24.0; Gbd1 /= b4; Gmd1 = 0.0;

	for (r = 3, s = 2; s <= 100; r += 2, s++) {
		fact = fact1 * r; fact1 = fact;
		//itself
		Gf = -2.0; Gf *= fact; Gf *= (2 * s);
		Gf *= pow(2, s); Gf *= Pn1;
		Gf /= pow(b, (2 * s + 1));
		//derivative of beta
		Gbd = 2.0; Gbd *= fact; Gbd *= (2 * s);
		Gbd *= pow(2, s); Gbd *= (2 * s + 1);
		Gbd *= Pn1; Gbd /= pow(b, (2 * s + 2));
		//derivative of mu
		Gmd = -2.0; Gmd *= fact; Gmd *= (2 * s);
		Gmd *= pow(2, s); Gmd *= Pnd1;
		Gmd /= pow(b, (2 * s + 1));

		Gf += Gf1; Gbd += Gbd1; Gmd += Gmd1;
		//Legendre polynomials
		Pn = (2 * s - 1); Pn *= m; Pn *= Pn1;
		tp = (s - 1); tp *= Pn0; Pn -= tp; Pn /= s;
		//derivative of Legendre polynomials
		Pnd = (2 * s - 1); Pnd *= m; Pnd *= Pnd1;
		tp = (2 * s - 1); tp *= Pn1; tr = (s - 1);
		tr *= Pnd0; Pnd += tp; Pnd -= tr; Pnd /= s;

		Pn0 = Pn1; Pn1 = Pn; Pnd0 = Pnd1; Pnd1 = Pnd;

		if (fabs(Gf - Gf1) < EPS)  break;

		Gf1 = Gf; Gbd1 = Gbd; Gmd1 = Gmd;
	}
}
void
Greenf::Asymp(const double b, const double m)
{
	double b5, sq7_4, sq9_4, c_32, s_52, c_12, c_72, c_52, s_32, s_12, s_72;
	double m2, b2, b4, b6, b3, m2_1, sq1_2, sq3_4, sq5_4, sq1_4, b2m, s32, c32, hh, asn;
	double G2, Gb, Gm, tr, tp;

	b2 = b * b;
	b3 = b2 * b; 
	b4 = b3 * b;
	b5 = b4 * b;
	b6 = b5 * b;
	Asymp1(m, b, b3, b4);

	m2 = m * m; m2_1 = (1.0 - m2);
	sq1_2 = pow(m2_1, 0.5);  sq7_4 = pow(m2_1, 1.75);
	sq9_4 = pow(m2_1, 2.25); sq1_4 = pow(m2_1, 0.25);
	sq3_4 = pow(m2_1, 0.75); sq5_4 = pow(m2_1, 1.25);
	b2m = -m * b2;
	hh = b2 * sq1_2 / 4.0;
	asn = asin(m);

	tp = 1.50; tp *= asn;
	tr = hh; tr += tp;
	s32 = sin(tr);
	c32 = cos(tr);
	tr = hh; tr -= tp;
	s_32 = sin(tr);
	c_32 = cos(tr);
	tp = 2.5; tp *= asn;
	tr = hh; tr -= tp;
	s_52 = sin(tr); c_52 = cos(tr);
	tp = 0.5; tp *= asn; tr = hh; tr -= tp;
	c_12 = cos(tr);
	s_12 = sin(tr); tp = 3.5; tp *= asn; tr = hh; tr -= tp;
	c_72 = cos(tr);
	s_72 = sin(tr);
	//intermediate values for itself
	G2 = b; G2 /= sq1_4; G2 *= s32;
	tr = c_12; tr /= 2.0;
	tr /= b; tr /= sq3_4;
	G2 += tr; tr = s_32;
	tr /= b3; tr /= sq3_4;
	G2 += tr; tr = -9.0;
	tr *= s_52; tr /= 8.0;
	tr /= b3; tr /= sq5_4;
	G2 += tr; tr = -9.0; tr *= c_72;
	tr /= b5; tr /= sq3_4;
	G2 += tr; tr = 24.0; tr *= c_52;
	tr /= b5; tr /= sq3_4;
	G2 += tr;
	//intermediate values for beta
	Gb = s32; Gb /= sq1_4;
	tr = c32; tr *= b2; tr *= sq1_4; tr /= 2.0;
	Gb += tr; tr = -c_12; tr /= 2.0; tr /= b2;
	tr /= sq3_4;
	Gb += tr; tr = -s_12; tr /= 4.0; tr /= sq1_4;
	Gb += tr; tr = -3.0; tr *= s_32; tr /= b4; tr /= sq3_4;
	Gb += tr; tr = c_32; tr /= 2.0; tr /= b2; tr /= sq1_4;
	Gb += tr; tr = 27.0; tr *= s_52; tr /= 8.0; tr /= b4; tr /= sq5_4;
	Gb += tr; tr = -9.0; tr *= c_52; tr /= 16.0; tr /= b2; tr /= sq3_4;
	Gb += tr; tr = 45.0; tr *= c_72; tr /= b6; tr /= sq3_4;
	Gb += tr;
	tr = 9.0; tr *= s_72; tr /= 2.0; tr /= b4; tr /= sq1_4;
	Gb += tr;
	tr = -120.0; tr *= c_52; tr /= b6; tr /= sq3_4;
	Gb += tr;
	tr = -12.0; tr *= s_52; tr /= b4; tr /= sq1_4;
	Gb += tr;
	//intermediate values for mu
	Gm = s32; Gm *= b; Gm *= m; Gm /= 2.0; Gm /= sq5_4;
	tr = c32; tr *= b;
	tp = b2m; tp += 6.0; tr *= tp; tr /= 4.0; tr /= sq3_4;
	Gm += tr;
	tr = 3.0; tr *= c_12; tr *= m; tr /= 4.0; tr /= b; tr /= sq7_4;
	Gm += tr;
	tr = -s_12; tp = b2m; tp -= 2.0; tr *= tp; tr /= 8.0; tr /= b; tr /= sq5_4;
	Gm += tr; tr = 3.0; tr *= s_32; tr *= m; tr /= 2.0; tr /= b3; tr /= sq7_4;
	Gm += tr; tr = c_32; tp = b2m; tp -= 6.0; tr *= tp; tr /= 4.0; tr /= b3; tr /= sq5_4;//forgotten value
	Gm += tr; tr = -45.0; tr *= s_52; tr *= m; tr /= 16.0; tr /= b3; tr /= sq9_4;
	Gm += tr; tr = -9.0; tr *= c_52; tp = b2m; tp -= 10.0; tr *= tp; tr /= 32.0;//forgoten value
	tr /= b5; tr /= sq7_4;
	Gm += tr; tr = -27.0; tr *= c_72; tr *= m; tr /= 2.0;
	tr /= b5; tr /= sq7_4;
	Gm += tr; tr = 9.0; tr *= s_72; tp = b2m; tp -= -14.0;
	tr *= tp; tr /= 4.0; tr /= b5; tr /= sq5_4;
	Gm += tr; tr = 72.0; tr *= c_52;
	tr *= m; tr /= 2.0; tr /= b5; tr /= sq7_4;
	Gm += tr; tr = -6.0; tr *= s_52;
	tp = b2m; tp -= 10.0; tr *= tp; tr /= b5; tr /= sq5_4;
	Gm += tr;

	//itself, beta, mu	
	tr = -b2; tr *= m; tr /= 4.0; tr = exp(tr); tr *= sqrt(2.0);
	G2 *= tr; Gf += G2;
	Gb *= tr; tp = -b; tp *= m; tp *= G2; tp /= 2.0;
	Gb += tp; Gbd += Gb;
	Gm *= tr; tp = -b2; tp *= G2; tp /= 4.0;
	Gm += tp; Gmd += Gm;

	//dG[0]=Gf; dG[1]=Gbd; dG[2]=Gmd;

	//++aint;
}


double Greenf::Bessel(const double k)
{
	const double tr  = 0.5 * k * k;
	const double tr2 = tr * tr;

	double bes  = 1.0 - tr2;   
	double term = tr2;        
	double fac  = 1.0;       
	double sign = -1.0;

	double bes_new;

	for (int s = 2; s <= 98; ++s)
	{
		fac  *= s;
		term *= tr2;
		sign = -sign;

		bes_new = bes + sign * term / (fac * fac);

		if (fabs(bes_new - bes) < EPS)
			break;

		bes = bes_new;
	}
	return bes;
}



void
Greenf::Fquad1(double& k, const double m, const double sq1_2, const double sq3_2, double* gk, double* dgm)
{
	int r;
	double bes, fk, tp, ts, tr, tr1;

	tp = PI; tp /= 4.0;
	ts = 2.0; ts /= PI;
	//for(k=0.0,r=0;k<intv;k+=0.05,r++){						  				  		
	for (k = 0.0, r = 0; k < 27; k += 0.05, r++) {
		bes = Bessel(k);
		//fk=k*bes-sqrt(2.0/PI)*cos(k*k-PI/4.0);
		tr = k; tr *= k;
		tr1 = tr; tr1 -= tp;
		tr1 = cos(tr1);
		tr1 *= sqrt(ts);
		tr1 *= -1.0;
		fk = k; fk *= bes; fk += tr1;
		//gk[r]=k*fk*exp(-k*k*m/sq1_2);
		tr1 = -1.0; tr1 *= tr;
		tr1 *= m; tr1 /= sq1_2;
		tr1 = exp(tr1); tr1 *= fk; tr1 *= k;
		gk[r] = tr1;
		//derivative for mu			  		
		//dgm[r]=-pow(k,3)*fk*exp(-k*k*m/sq1_2)/sq3_2;
		tr *= -1.0; tr *= tr1; tr /= sq3_2;
		dgm[r] = tr;
		//if(k!=0.0 && fabs(gk[r])<1.0e-7 && fabs(dgm[r])<1.0e-7) break;						 				 
		//if(k!=0.0 && fabs(tr1)<1.0e-7 && fabs(tr)<1.0e-7) break;						 				 
		if (k != 0.0 && fabs(tr1) < EPS && fabs(tr) < EPS) break;
	}
}
void
Greenf::Fquad2(double& delt, double& a2, double& a2b, double& db, double& a2m, double& dm, double& a3, double& a3b, double& a3m)
{
	double si2, co2, delt3, ddelt, tr, tp;

	tr = 2.0; tr *= delt;
	si2 = sin(tr);
	//si2=sin(2.0*delt);
	co2 = cos(tr);
	//co2=cos(2.0*delt);
	delt3 = pow(delt, 3.0);

	tr = -2.0; tr *= si2;
	a2 = 3.0;  a2 += co2;
	a2 *= delt; a2 += tr; a2 /= delt3;
	//a2=(delt*(3.0+co2)-2.0*si2)/delt3;	
	tr *= delt;
	tp = 1.0; tp -= co2; tp *= 3.0;
	tp += tr; tp /= delt3;
	ddelt = -3.0; ddelt *= a2;
	ddelt /= delt; ddelt += tp;
	//ddelt=(-3.0/delt*a2+(3.0*(1.0-co2)-2.0*delt*si2)/delt3);			
	a2b = ddelt; a2b *= db;
	//a2b=ddelt*db;
	a2m = ddelt; a2m *= dm;
	//a2m=ddelt*dm;		
	tr = sin(delt); a3 = -1.0;
	a3 *= cos(delt); a3 *= delt;
	a3 += tr; a3 *= 4.0; a3 /= delt3;
	//a3=4.0*(sin(delt)-delt*cos(delt))/delt3;		
	tr *= 4.0;
	tp = -3.0; tp *= a3; tp *= delt;
	tp += tr;	 tp /= pow(delt, 2.0);
	//tp=(-3.0*delt*a3+4.0*sin(delt))/pow(delt,2.0);					
	a3b = tp; a3b *= db;
	//a3b=tp*db;		
	a3m = tp; a3m *= dm;
	//a3m=tp*dm;
}
void
Greenf::Fquad3(double& delt, double& a2, double& a2b, double& db, double& a2m, double& dm, double& a3, double& a3b, double& a3m)
{
	int s, r, r1;
	//double fact[110],al10,al20,al30,al40,al50,al60,alf,dal,tr,ts,tp;
	double al10, al20, al30, al40, al50, al60, alf, dal, tr, ts, tp;

	//factor[1]=1.0;   	
	//for(r=2,s=2;s<=100;r++,s++) 
	  //  factor[s]=factor[s-1]*r;

	al10 = 2.0 / 3.0; al20 = 4.0 / 3.0;
	al30 = 0.0; al40 = 0.0; al50 = 0.0; al60 = 0.0;

	//for(r=4,s=5,r1=1;r1<intv;r1++,r+=2,s+=2){
	for (r = 4, s = 5, r1 = 1; r1 <= 80; r1++, r += 2, s += 2) {
		//itself
		//alf=pow(-1,(r1+1))*pow(2,r)/factor[r]*(1.0-4.0/s);
		alf = -4.0; alf /= s; alf += 1.0;
		alf *= pow(2.0, r); tr = 1.0;
		tr += r1; alf *= pow(-1.0, tr);
		alf /= factor[r];
		//dal=pow(-1,r1)*4.0/factor[r]*(1.0-1.0/s);
		dal = -1.0; dal /= s;
		dal += 1.0; dal *= 4.0;
		dal *= pow(-1.0, r1);
		dal /= factor[r];
		//a2=alf*pow(delt,(r-2));		
		//a2+=al10;			
		tr = -2.0; tr += r;
		tp = tr; tr = pow(delt, tr);
		a2 = alf; a2 *= tr;	a2 += al10;
		//a3=dal*pow(delt,(r-2));		
		//a3+=al20;			
		a3 = dal; a3 *= tr;	a3 += al20;
		//derivative of beta
		//a2b=alf*(r-2)*db*pow(delt,(r-3));		
		//a2b+=al30;		
		ts = tp; ts -= 1.0;
		ts = pow(delt, ts);
		ts *= db; ts *= tp;
		a2b = alf; a2b *= ts; a2b += al30;
		//a3b=dal*(r-2)*db*pow(delt,(r-3));		
		//a3b+=al40;			
		a3b = dal; a3b *= ts; a3b += al40;
		//derivative of mu
		//a2m=alf*(r-2)*dm*pow(delt,(r-3));		
		//a2m+=al50;			
		ts /= db; ts *= dm; a2m = alf;
		a2m *= ts; a2m += al50;
		//a3m=dal*(r-2)*dm*pow(delt,(r-3));		
		//a3m+=al60;
		a3m = dal; a3m *= ts; a3m += al60;

		//if(fabs(a2-al10)<1.0e-12 && fabs(a3-al20)<1.0e-12 && fabs(a2b-al30)<1.0e-12 &&
		  // fabs(a3b-al40)<1.0e-12 && fabs(a2m-al50)<1.0e-12 && fabs(a3m-al60)<1.0e-12) 	break;		

		if (fabs(a2 - al10) < EPS && fabs(a3 - al20) < EPS && fabs(a2b - al30) < EPS &&
			fabs(a3b - al40) < EPS && fabs(a2m - al50) < EPS && fabs(a3m - al60) < EPS) 	break;

		al10 = a2; al20 = a3; al30 = a2b; al40 = a3b; al50 = a2m; al60 = a3m;
	}
}
void
Greenf::Fsummation(double& thet, double* gk, double* dgm, int& tp, double& delt, double& kmin, double& s2n, const double sq1_4, double& s2nb, const double b, const double m, const double sq5_4, double& s2nm)
{
	double tn, tr, ts, tk, sg;

	tn = gk[tp];
	tr = thet; tr *= kmin;
	ts = cos(tr);	tr = sin(tr);
	sg = tn; sg *= tr; s2n += sg;
	//sgk=gk[2*s]*sin(kmin*thet);
	//s2n+=sgk;
	//for derivative of beta	    		
	sg = tn; sg *= kmin; sg *= ts;
	sg /= sq1_4; s2nb += sg;
	//sgkb=gk[2*s]*kmin*cos(kmin*thet)/sq1_4;    
	//s2nb+=sgkb;
	//for derivative of mu	    				
	sg = dgm[tp];
	sg *= tr; tk = tn; tk *= b;
	tk *= kmin; tk *= m; tk *= ts;
	tr = 2.0; tr *= sq5_4;
	tk /= tr; sg += tk; s2nm += sg;
	//sgkm=dgm[2*s]*sin(kmin*thet);
	//sgkm+=gk[2*s]*b*kmin*m/(2.0*sq5_4)*cos(kmin*thet);
	//s2nm+=sgkm;
}
void
Greenf::Fquad(const double b, const double m)
{
	int s, tp;
	double kmin, k, tr, ts, tn, tk, sg;
	double gk[95], h1, a2, a3, s2n1, s2n;
	double dgm[95], a2m, a3m, s2nm, s2n1m, G3;
	double a2b, a3b, s2nb, s2n1b;
	double s32, c32;
	double h2, h3;
	double delt, thet;
	double db, dm;

	const double m2    = m * m;
	const double m2_1  = 1.0 - m2;
	const double sq1_2 = sqrt(m2_1);
	const double sq1_4 = sqrt(sq1_2);
	const double sq3_2 = m2_1 * sq1_2;
	const double sq3_4 = sq1_2 * sq1_4; 
	const double sq5_4 = m2_1 * sq1_4;
	const double b2    = b * b;
	const double b3    = b2 * b;
	const double b2m   = -m * b2;

	Fquad1(k, m, sq1_2, sq3_2, gk, dgm);

	tr = 1.5 * asin(m);
	s32 = b2; s32 *= sq1_2; s32 /= 4.0; s32 += tr;
	c32 = cos(s32);
	s32 = sin(s32);

	h1 = 0.05;
	h2 = h1 * h1;
	h3 = h2 * h1;
	thet = b; thet /= sq1_4;
	delt = thet; delt *= h1;
	db = h1; db /= sq1_4;
	dm = b; dm *= h1; dm *= m; tr = 2.0; tr *= sq5_4; dm /= tr;

	if (delt > 0.35)
		Fquad2(delt, a2, a2b, db, a2m, dm, a3, a3b, a3m);
	else
		Fquad3(delt, a2, a2b, db, a2m, dm, a3, a3b, a3m);

	//for even summation				 	
	s2n = 0.0; s2nb = 0.0; s2nm = 0.0;
	for (s = 0, kmin = 0; kmin <= k; kmin += 2 * h1, s++) {
		tp = 2; tp *= s;
		Fsummation(thet, gk, dgm, tp, delt, kmin, s2n, sq1_4, s2nb, b, m, sq5_4, s2nm);
	}
	//for odd summation
	s2n1 = 0.0; s2n1b = 0.0; s2n1m = 0.0;
	for (s = 1, kmin = h1; kmin <= k; kmin += 2 * h1, s++) {
		tp = 2; tp *= s; tp -= 1;
		Fsummation(thet, gk, dgm, tp, delt, kmin, s2n1, sq1_4, s2n1b, b, m, sq5_4, s2n1m);
	}
	//itself
   //Gf=sqrt(2)*b*exp(b2m/4.0)*s32/sq1_4;                  				  			
   //Gf+=G3;	
   //G3=4.0*h1/sq3_4*(a2*s2n+a3*s2n1);				   
   //beta
   //Gbd=4.0*h1/sq3_4*(a2b*s2n+a2*s2nb+a3b*s2n1+a3*s2n1b);                   	
   //Gbd+=sqrt(2)/sq1_4*exp(b2m/4.0)*((1.0+b2m/2.0)*s32+b2*sq1_2/2.0*c32);				   				   						
   //mu
   //Gmd=1.5*m/m2_1*G3;
   //Gmd+=4.0*h1/sq3_4*(a2m*s2n+a2*s2nm+a3m*s2n1+a3*s2n1m);	
   //Gmd+=sqrt(2)*b/sq1_4*exp(b2m/4.0)*((m/(2.0*m2_1)-b2/4.0)*s32+(b2m+6.0)/(4.0*sq1_2)*c32);				   				   				   			   		

   //itself
	ts = 4.0; ts *= h1; ts /= sq3_4;
	//ts=4.0*h1/sq3_4
	tn = b2m; tn /= 4.0; tn = exp(tn);
	tn *= sqrt(2.0); tn /= sq1_4;
	sg = tn; tn *= b;
	//tn=sqrt(2)*b*exp(b2m/4.0)/sq1_4; 
	tr = a3; tr *= s2n1; G3 = a2; G3 *= s2n;
	G3 += tr;	G3 *= ts;
	//G3=4.0*h1/sq3_4*(a2*s2n+a3*s2n1);				   
	Gf = tn; Gf *= s32;	Gf += G3;
	//Gf=sqrt(2)*b*exp(b2m/4.0)*s32/sq1_4;                  				  			
	//Gf+=G3;	
	//beta
	tr = a3; tr *= s2n1b; Gbd = a3b;
	Gbd *= s2n1; Gbd += tr; tr = a2;
	tr *= s2nb; Gbd += tr; tr = a2b;
	tr *= s2n; Gbd += tr; Gbd *= ts;
	//Gbd=4.0*h1/sq3_4*(a2b*s2n+a2*s2nb+a3b*s2n1+a3*s2n1b);                   	
	tr = b2; tr *= c32; tr *= sq1_2;
	tr /= 2.0; tk = b2m; tk /= 2.0;
	tk += 1.0; tk *= s32; tk += tr;
	tk *= sg;	 Gbd += tk;
	//Gbd+=sqrt(2)/sq1_4*exp(b2m/4.0)*((1.0+b2m/2.0)*s32+b2*sq1_2/2.0*c32);				   				   						
	//mu
	Gmd = 1.5; Gmd *= m;
	Gmd *= G3; Gmd /= m2_1;
	//Gmd=1.5*m/m2_1*G3;
	tr = a3; tr *= s2n1m; tk = a3m;
	tk *= s2n1; tk += tr; tr = a2;
	tr *= s2nm; tk += tr; tr = a2m;
	tr *= s2n; tk += tr; tk *= ts; Gmd += tk;
	//Gmd+=4.0*h1/sq3_4*(a2m*s2n+a2*s2nm+a3m*s2n1+a3*s2n1m);	
	tr = b2m; tr += 6.0; tr *= c32;
	tr /= 4.0; tr /= sq1_2; tk = m;
	tk /= 2.0; tk /= m2_1; k = -0.25;
	k *= b2; k += tk; k *= s32;
	k += tr; k *= tn; Gmd += k;
	//Gmd+=sqrt(2)*b/sq1_4*exp(b2m/4.0)*((m/(2.0*m2_1)-b2/4.0)*s32+(b2m+6.0)/(4.0*sq1_2)*c32);				   				   				   			   		
	//dG[0]=Gf; dG[1]=Gbd; dG[2]=Gmd;	

	//++fint;
}
void
Greenf::Cplxerfbes1(double& x, double& Im, double& Imb, const double m, const double b, double& Imm, const std::vector<double>& cp)
{
	int i, n;
	double tp, ts, tr, tr1;

	n = cp.size();
	tp = x * x;
	for (i = 0; i < n - 1; i += 2) {
		ts = cp[i + 1];
		ts *= -1;
		ts += tp;
		tr = cp[i];
		tr1 = tr;
		tr1 /= ts;
		Im += tr1;//itself
		Imb += tr1;
		ts *= ts;
		tr /= ts;
		tr *= tp;
		tr *= -2.0;
		Imb += tr;//derivative			
	}
	Im *= x;//itself
	Imm = Imb;//der mu
	Imb /= 2.0; Imb /= sqrt(m);//der beta		
	Imm *= -1.0;
	Imm *= b;
	Imm /= 4.0;
	Imm /= pow(m, 1.5);
}
void
//Greenf::Cplxerfbes(double&b,double&m,Vector&Rw)
Greenf::Cplxerfbes(const double b, const double m)
{
	int i, n;
	double Im = 0.0, Imb = 0.0, Imm, x, y, tb, tm, tp, tr, tr1, ts;

	x = b;
	x /= 2.0;
	x /= sqrt(m);
	y = 0.0;
	tr = -1.0;
	tr *= (x * x);
	tr = exp(tr);
	Rw[0] = tr;
	tb = x;
	tb *= tr;
	tm = tb;
	tb *= -1.0;
	tb /= sqrt(m);
	Rw[1] = tb;
	tm *= b;
	tm /= 2.0;
	tm /= pow(m, 1.5);//forgotten value
	Rw[2] = tm;

	if (x > 6.0)
		Cplxerfbes1(x, Im, Imb, m, b, Imm, cp1);
	else if (x > 5.0 && x <= 6.0)
		Cplxerfbes1(x, Im, Imb, m, b, Imm, cp2);
	else {
		n = anh.size();
		tp = x;
		tp /= 5.0;
		tr1 = tp;
		tr1 *= tp;
		tr1 *= -1.0;
		tr1 += 1.0;
		tr1 = sqrt(tr1);
		tr1 *= 5.0;
		tp = acos(tp);
		for (i = 0; i < n; i++) {
			ts = 2.0;
			ts *= i;
			ts += 1.0;
			tr = ts;
			ts *= tp;
			tm = ts;
			tm = cos(tm);
			tm *= anh[i];
			Im += tm;//itself
			ts = sin(ts);
			ts *= anh[i];
			ts *= tr;
			ts /= tr1;
			Imb += ts;//derivative			
		}
		Im *= 2.0; Im /= sqrt(PI);//itself 
		Imm = Imb;//der mu
		tp = PI;
		tp *= m;
		Imb /= sqrt(tp);//der beta
		Imm *= -1.0;
		Imm *= b;
		Imm /= 2.0;
		Imm /= sqrt(PI);
		Imm /= pow(m, 1.5);
	}
	Rw[3] = Im; Rw[4] = Imb; Rw[5] = Imm;
}
void
Greenf::Besexp1(int r, int k, int& i, int& j, double& tr, double* Re,
	double* Im, const double b, double& tt, double* Im1)
{
	double ti;

	tr = Re[j];
	tr *= i;
	ti = Im[i];
	ti *= b;
	if (r == 0) {
		tr -= ti;
		if (k == 0)
			tr -= Im1[i];
	}
	else if (r == 1) {
		tr += ti;
		if (k == 1)
			tr += Im1[i];
	}
	tr /= tt;
}
void
Greenf::Besexp2(int k, int& i, int& j, double& tr, double* Re, double* Im,
	const double m, const double b, double& tm, double* Im1, double* Re1)
{
	double ti;
	tr = m;
	tr *= Re[j];
	tr -= Re1[j];
	tr *= i;
	ti = m;
	ti *= Im[i];
	ti -= Im1[i];
	ti *= b;
	if (k == 0) tr -= ti;
	else tr += ti;
	tr /= tm;
}
void
Greenf::Betabessel(double& Gf, double& tr, double& sq1_2, int& s, double* ImI)
{
	int i;
	double ti;

	Gf = sq1_2;
	Gf /= 3.0;
	ti = 2.0;
	ti *= s;
	Gf = pow(Gf, ti);
	i = (int)tr;
	i *= s;
	i += 2;
	ti = ImI[i];
	Gf *= ti;
	Gf *= an[s];
	Gf *= tr;
}
void
Greenf::Besexp(const double b, const double m)
{
	int i, j, s, v, n;
	double ReI[24], ImI[24];
	double ReIb[24], ImIb[24];
	double ReIm[24], ImIm[24];
	double m2, m2_1, sq1_2;
	double Gf0, Gbd0, Gmd0, tp, tr, ti, tt, tm;
	//Vector Rw(6);


	m2 = m * m;
	m2_1 = (1.0 - m2);
	m2_1 = (m2_1 < 1.0e-6 ? 0.0 : m2_1);
	sq1_2 = pow(m2_1, 0.5);

	Cplxerfbes(b, m);

	//real part
	tp = PI;
	tp /= m;
	tp = sqrt(tp);
	tp *= 0.5;
	tr = Rw[0];
	ti = Rw[3];
	tr *= tp;
	ReI[0] = tr;
	//ReI[0]=0.5*sqrt(PI/m)*Rw[0];
	ti *= tp;
	ImI[0] = ti;
	//ImI[0]=0.5*sqrt(PI/m)*Rw[3];
	tr = Rw[1];
	ti = Rw[4];
	tr *= tp;
	ReIb[0] = tr;
	//ReIb[0]=0.5*sqrt(PI/m)*Rw[1];	
	ti *= tp;
	ImIb[0] = ti;
	//ImIb[0]=0.5*sqrt(PI/m)*Rw[4];	
	tt = 2.0;
	tt *= m;
	tp = sqrt(PI);
	tp /= 4.0;
	tp /= pow(m, 1.5);
	tr = tt;
	ti = tr;
	tr *= Rw[3];
	tr -= Rw[0];
	tr *= tp;
	ReIm[0] = tr;
	//ReIm[0]=sqrt(PI)*(2.0*m*Rw[3]-Rw[0])/(4.0*pow(m,1.5));	
	ti *= Rw[5];
	ti -= Rw[3];
	ti *= tp;
	ImIm[0] = ti;
	//ImIm[0]=sqrt(PI)*(2.0*m*Rw[5]-Rw[3])/(4.0*pow(m,1.5));

	//imaginary part	
	tr = -1.0;
	tr *= b;
	tr *= ImI[0];
	tr += 1.0;
	tr /= tt;
	ReI[1] = tr;
	//ReI[1]=(1.0-b*ImI[0])/(2.0*m);	
	tr = b;
	tr *= ReI[0];
	tr /= tt;
	ImI[1] = tr;
	//ImI[1]=b*ReI[0]/(2.0*m);
	tr = b;
	ti = tr;
	tr *= ImIb[0];
	tr += ImI[0];
	tr /= tt;
	tr *= -1.0;
	ReIb[1] = tr;
	//ReIb[1]=-(ImI[0]+b*ImIb[0])/(2.0*m);		
	tr *= -1.0;
	ti *= ReIb[0];
	ti += ReI[0];
	ti /= tt;
	ImIb[1] = ti;
	//ImIb[1]=(ReI[0]+b*ReIb[0])/(2.0*m);	
	tm = tt;
	tm *= m;
	tr = b;
	ti = tr;
	tr *= ImI[0];
	ti *= m;
	ti *= ImIm[0];
	tr -= ti;
	tr -= 1.0;
	tr /= tm;
	ReIm[1] = tr;
	//ReIm[1]=(b*ImI[0]-b*m*ImIm[0]-1.0)/(2.0*m*m);	
	tr = m;
	tr *= ReIm[0];
	tr -= ReI[0];
	tr *= b;
	tr /= tm;
	ImIm[1] = tr;
	//ImIm[1]=b*(m*ReIm[0]-ReI[0])/(2.0*m*m);	

	n = anh.size();
	for (v = 2; v < n; v++) {
		i = v; i -= 1;
		j = i; j -= 1;
		Besexp1(0, 1, i, j, tr, ReI, ImI, b, tt, ReIb);
		ReI[v] = tr;
		//ReI[v]=((v-1)*ReI[v-2]-b*ImI[v-1])/(2.0*m);		
		Besexp1(1, 0, i, j, tr, ImI, ReI, b, tt, ImIb);
		ImI[v] = tr;
		//ImI[v]=((v-1)*ImI[v-2]+b*ReI[v-1])/(2.0*m);	
		Besexp1(0, 0, i, j, tr, ReIb, ImIb, b, tt, ImI);
		ReIb[v] = tr;
		//ReIb[v]=((v-1)*ReIb[v-2]-b*ImIb[v-1]-ImI[v-1])/(2.0*m);
		Besexp1(1, 1, i, j, tr, ImIb, ReIb, b, tt, ReI);
		ImIb[v] = tr;
		//ImIb[v]=((v-1)*ImIb[v-2]+b*ReIb[v-1]+ReI[v-1])/(2.0*m);

		Besexp2(0, i, j, tr, ReIm, ImIm, m, b, tm, ImI, ReI);
		ReIm[v] = tr;
		//ReIm[v]=-b*(m*ImIm[v-1]-ImI[v-1])/(2.0*m*m)+(v-1)*(m*ReIm[v-2]-ReI[v-2])/(2.0*m*m);

		Besexp2(1, i, j, tr, ImIm, ReIm, m, b, tm, ReI, ImI);
		ImIm[v] = tr;
		//ImIm[v]=b*(m*ReIm[v-1]-ReI[v-1])/(2.0*m*m)+(v-1)*(m*ImIm[v-2]-ImI[v-2])/(2.0*m*m);
	}
	//calculation of source potential
	tr = 4.0;
	Gf0 = tr;
	Gf0 *= ImI[2];
	Gbd0 = tr;
	Gbd0 *= ImIb[2];
	Gmd0 = tr;
	Gmd0 *= ImIm[2];
	if (m2_1 == 0.0) {
		Gf = Gf0;
		Gbd = Gbd0;
		ti = 8.0 / 9.0;
		ti *= an[1];
		ti *= ImI[6];
		Gmd0 -= ti;
		Gmd = Gmd0;
		//Gmd=4.0*ImIm[2]-8.0/9.0*an[1]*ImI[6];	
	}
	else {
		for (s = 1; s < 6; s++) {
			Betabessel(Gf, tr, sq1_2, s, ImI);
			//Gf=4.0*an[s]*pow((sq1_2/3.0),2*s)*ImI[4*s+2];				      			
			Betabessel(Gbd, tr, sq1_2, s, ImIb);
			//Gbd=4.0*an[s]*pow((sq1_2/3.0),2*s)*ImIb[4*s+2];					  
			tm = sq1_2;
			tm /= 3.0;
			tp = 2;
			tp *= s;
			ti = pow(tm, tp);
			i = (int)tr;
			i *= s;
			i += 2;
			ti *= ImIm[i];
			Gmd = ImI[i];
			tp -= 1.0;
			Gmd *= pow(tm, tp);
			Gmd *= -2.0;
			Gmd *= s;
			Gmd *= m;
			Gmd /= 3.0;
			Gmd /= sq1_2;
			Gmd += ti;
			Gmd *= an[s];
			Gmd *= tr;
			//Gmd=4.0*an[s]*(-2.0*s*m/(3.0*sq1_2)*pow((sq1_2/3.0),(2*s-1))*ImI[4*s+2]+pow((sq1_2/3.0),(2*s))*ImIm[4*s+2]);				  

			Gf += Gf0; Gbd += Gbd0; Gmd += Gmd0;
			Gf0 = Gf;  Gbd0 = Gbd;  Gmd0 = Gmd;
		}
	}
	//dG[0]=Gf; dG[1]=Gbd; dG[2]=Gmd;	

	//++bint;
}
void
Greenf::AsBexp1(const double b, const double m, double& s2, double& s2b, double& s2m)
{
	int r, s;
	double s21, s2b1, s2m1, gam1, gam;

	double b3 = b * b * b;
	double b4 = b3 * b;
	s21 = 2.0; s21 /= b3;
	s2b1 = -6.0; s2b1 /= b4;
	s2m1 = 0.0; gam1 = 1.0;
	//for(r=3,s=2;s<intv;r+=2,s++){ 
	for (r = 3, s = 2; s <= 33; r += 2, s++) {
		gam = gam1 * r; gam1 = gam;
		//itself
		s2 = gam; s2 *= (2 * s);
		s2 *= pow((2.0 * m), (s - 1));
		s2 /= (pow(b, (2 * s + 1)));
		//derivative of beta
		s2b = -gam; s2b *= (2 * s); s2b *= (2 * s + 1);
		s2b *= pow((2.0 * m), (s - 1));
		s2b /= (pow(b, (2 * s + 2)));
		//derivative of mu
		s2m = gam; s2m *= (4.0 * s); s2m *= (s - 1);
		s2m *= pow((2.0 * m), (s - 2));
		s2m /= (pow(b, (2 * s + 1)));

		s2 += s21; s2b += s2b1; s2m += s2m1;
		//if(fabs(s2m-s2m1)<1.0e-10) break;
		if (fabs(s2m - s2m1) < EPS) break;
		s21 = s2; s2b1 = s2b; s2m1 = s2m;
	}
}
void
Greenf::AsBexp2(const double b, const double m, double& s6, double& s6b, double& s6m)
{
	int r, s, kmax, v;
	double s21, s2b1, s2m1, gam1, gam;

	s21 = 720.0; s21 /= pow(b, 7.0);
	s2b1 = -5040.0; s2b1 /= pow(b, 8.0);
	s2m1 = 0.0; gam1 = 15.0;

	//for(kmax=3,v=3,r=7,s=4;s<intv;kmax+=v,v++,r+=2,s++){
	for (kmax = 3, v = 3, r = 7, s = 4; s <= 33; kmax += v, v++, r += 2, s++) {
		gam = gam1 * r; gam1 = gam;
		//itself
		s6 = 16.0; s6 *= gam; s6 *= s;
		s6 *= pow((2.0 * m), (s - 3));
		s6 /= (pow(b, (2 * s + 1))); s6 *= kmax;
		//derivative of beta
		s6b = -16.0; s6b *= gam; s6b *= s;
		s6b *= (2 * s + 1); s6b *= pow(2.0 * m, (s - 3));
		s6b /= pow(b, (2 * s + 2)); s6b *= kmax;
		//derivative of mu
		s6m = 32.0; s6m *= gam; s6m *= s;
		s6m *= (s - 3); s6m *= pow(2.0 * m, (s - 4));
		s6m /= pow(b, (2 * s + 1)); s6m *= kmax;

		s6 += s21; s6b += s2b1; s6m += s2m1;
		//if(fabs(s6m-s2m1)<1.0e-10) break;
		if (fabs(s6m - s2m1) < EPS) break;
		s21 = s6; s2b1 = s6b; s2m1 = s6m;
	}
}
void
Greenf::AsBexp3(const double b, const double m, double& s10, double& s10b, double& s10m)
{
	int r, s, kmax, v, q, qq;
	double s21, s2b1, s2m1, gam1, gam;

	s21 = 3628800.0; s21 /= pow(b, 11.0);
	s2b1 = -39916800.0; s2b1 /= pow(b, 12.0);
	s2m1 = 0.0; gam1 = 945.0;

	//for(kmax=3,v=3,q=10,qq=5,r=11,s=6;s<intv;kmax+=v,v++,q+=qq,qq++,r+=2,s++){
	for (kmax = 3, v = 3, q = 10, qq = 5, r = 11, s = 6; s <= 33; kmax += v, v++, q += qq, qq++, r += 2, s++) {
		gam = gam1 * r; gam1 = gam;
		//itself
		s10 = 128.0; s10 *= gam; s10 *= s;
		s10 *= pow((2.0 * m), (s - 5));
		s10 /= (pow(b, (2 * s + 1))); s10 *= q; s10 *= kmax;
		//derivative of beta
		s10b = -128.0; s10b *= gam; s10b *= s;
		s10b *= (2 * s + 1); s10b *= pow((2.0 * m), (s - 5));
		s10b /= (pow(b, (2 * s + 2))); s10b *= q; s10b *= kmax;
		//derivative of mu
		s10m = 256.0; s10m *= gam; s10m *= s;
		s10m *= (s - 5); s10m *= pow((2.0 * m), (s - 6));
		s10m /= (pow(b, (2 * s + 1))); s10m *= q; s10m *= kmax;

		s10 += s21; s10b += s2b1; s10m += s2m1;
		//if(fabs(s10m-s2m1)<1.0e-10) break;
		if (fabs(s10m - s2m1) < EPS) break;
		s21 = s10; s2b1 = s10b; s2m1 = s10m;
	}
}
void
Greenf::AsBexp(const double b, const double m)
{
	//double Gf,Gbd,Gmd,s2,s6,s10,s2b,s6b,s10b,s2m,s6m,s10m;	
	double s2, s6, s10, s2b, s6b, s10b, s2m, s6m, s10m;
	double m2, m2_1, sq1_2, tp, tr, tp1, tr1;

	m2 = m * m;
	m2_1 = (1.0 - m2);
	m2_1 = ((m2_1 < 1.0e-6) ? (0.0) : (m2_1));
	sq1_2 = sqrt(m2_1);

	AsBexp1(b, m, s2, s2b, s2m);
	AsBexp2(b, m, s6, s6b, s6m);
	AsBexp3(b, m, s10, s10b, s10m);
	//calculation of source potential
	if (m2_1 == 0.0) {
		Gf = -4.0;  Gf *= s2;
		Gbd = -4.0; Gbd *= s2b;
		tp = -4.0;  tp *= s2m;
		Gmd = 8.0 / 9.0;
		Gmd *= an[1]; Gmd *= s6; Gmd += tp;
	}
	else {
		Gf = -4.0; Gf *= s2;//itsels
		tp = -4.0; tp *= an[1];
		tp *= pow((sq1_2 / 3.0), 2.0); tp *= s6;
		tr = -4.0; tr *= an[2];
		tr *= pow((sq1_2 / 3.0), 4.0); tr *= s10;
		Gf += tp; Gf += tr;
		//derivative of beta
		Gbd = -4.0; Gbd *= s2b;
		tp = -4.0;  tp *= an[1];
		tp *= pow((sq1_2 / 3.0), 2); tp *= s6b;
		tr = -4.0; tr *= an[2];
		tr *= pow((sq1_2 / 3.0), 4.0); tr *= s10b;
		Gbd += tp; Gbd += tr;
		//derivative of mu
		Gmd = -4.0; Gmd *= s2m;
		tp = 8.0;   tp *= an[1];
		tp *= m / (3.0 * sq1_2); tp *= sq1_2 / 3.0;
		tp *= s6; tr = -4.0; tr *= an[1];
		tr *= pow((sq1_2 / 3.0), 2.0); tr *= s6m;
		tp1 = 8.0; tp1 *= an[2];
		tp1 *= 2.0; tp1 *= m / (3.0 * sq1_2);
		tp1 *= pow((sq1_2 / 3.0), 3.0); tp1 *= s10;
		tr1 = -4.0; tr1 *= an[2];
		tr1 *= pow((sq1_2 / 3.0), 4.0); tr1 *= s10m;
		Gmd += tp; Gmd += tr; Gmd += tp1; Gmd += tr1;
	}
	//dG[0]=Gf; dG[1]=Gbd; dG[2]=Gmd;

	//++abint;
}
void
Greenf::GreenFunctionCal(const double b, const double m)
{
	//if (m > 0.001 && m < 0.7 && b < 8.8)			Series1(b, m);
	//else if (m <= 0.001 && b < 8.8)					Series2(b, m);
	//else if (b >= 8.8 && m < 0.7)					Asymp(b, m);
	//else if (m >= 0.7 && m < 0.97)					Fquad(b, m);
	//else if ((m >= 0.97 || m == 1.0) && b < 9.5)	Besexp(b, m);
	//else   											AsBexp(b, m);

	auto green = TDGF_ba(b,m);
	Gf = green[0];
	Gbd = green[1];
	Gmd = green[2];
}
