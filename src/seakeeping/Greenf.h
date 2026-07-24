#pragma once

#include <array>
#include <vector>

class Greenf
{
public:
	Greenf();
	void GreenFunctionCal(const double, const double);

	//Gf:格林函数本身，Gbd:格林函数对beta的导数，Gmd:格林函数对mu的导数
	double Gf, Gbd, Gmd;

private:
	void Series1(const double b, const double m);
	void Series2(const double b, const double m);
	void Asymp  (const double b, const double m);
	void Fquad  (const double b, const double m);
	void Besexp (const double b, const double m);
	void AsBexp (const double b, const double m);

	double Bessel   (const double k);
	void AsBexp1    (const double b, const double m, double&, double&, double&);
	void AsBexp2    (const double b, const double m, double&, double&, double&);
	void AsBexp3    (const double b, const double m, double&, double&, double&);
	void Asymp1     (const double m, const double b, double&, double&);
	void Fquad1     (double& k, const double m, const double, const double, double*, double*);
	void Fquad2     (double&, double&, double&, double&, double&, double&, double&, double&, double&);
	void Fquad3     (double&, double&, double&, double&, double&, double&, double&, double&, double&);
	void Fsummation (double&, double*, double*, int&, double&, double&, double&, const double, double&, const double b, const double m, const double, double&);

	void Cplxerfbes (const double b, const double m);
	void Cplxerfbes1(double&, double&, double&, const double m, const double b, double&, const std::vector<double>&);

	void Besexp1    (int, int, int&, int&, double&, double*, double*, const double b, double&, double*);
	void Besexp2    (int, int&, int&, double&, double*, double*, const double b, const double m, double&, double*, double*);
	void Betabessel (double&, double&, double&, int&, double*);

	//void GreenfHummerPoints();

	static const std::array<double, 47> factor;
	static const std::array<double, 23> anh;
	static const std::array<double, 6 > an;
	static const std::vector<double> cp1;
	static const std::vector<double> cp2;

	std::array<double, 6 > Rw;

	//std::vector<double> cp1;
	//std::vector<double> cp2;

};

