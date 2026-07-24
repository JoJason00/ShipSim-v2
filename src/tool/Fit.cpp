#include "Fit.h"
#include "../const/Const.h"
#include <iostream>
#include <iomanip>
#include <fstream>

Fit::Fit() :matrix(), dt(0), omiga(0), Amp(0), cosfit()
{
}

void Fit::setdata(const Eigen::MatrixXd& data, double dt, double we, double Amp, std::vector<std::string> name)
{
	matrix      = &data;
	this->dt    = dt;
	this->omiga = we;
	this->Amp   = Amp;
	this->names = name;

	int n = data.cols();
	cosfit.resize(n);
}


int Fit::cut(const Eigen::VectorXd& signal, double dt, double omiga)
{
	const int N = signal.size();
	const int period_points = int(2.0 * PI / (omiga * dt));

	const int keep_periods = 4;  
	const int keep_points = keep_periods * period_points;

	if (keep_points >= N)
		return N / 3;

	return N - keep_points;
}


void Fit::run()
{
    const int TS = (*matrix).rows();
    const int nDOF = (*matrix).cols();

    const double w = omiga;
    const double dt_local = dt;
    //const double invAmp = 1.0 / Amp;

	double t, y, wt, s, c;

    for (int j = 0; j < nDOF; ++j)
    {
        const Eigen::VectorXd& col = (*matrix).col(j);
        const int i0 = cut(col, dt_local, w);

		std::cout << "\ni0==" << i0 << "\t" << "n==" << TS << std::endl;

		double ss = 0, cc = 0, sc = 0, ys = 0, yc = 0;

        for (int i = i0; i < TS; ++i)
        {
            t = i * dt_local;
            y = col(i);

            wt = w * t;
            s = sin(wt);
            c = cos(wt);

			ss += s * s; cc += c * c; sc += s * c;
			ys += y * s; yc += y * c;
        }

        // =====       С   ˽  =====
		const double det = ss * cc - sc * sc;
		const double A = (ys * cc - yc * sc) / det;
		const double B = (yc * ss - ys * sc) / det;

        //const double amplitude = std::sqrt(A * A + B * B) * invAmp;
		const double amplitude = std::sqrt(A * A + B * B);

        double phase = std::atan2(-A, B);
        if (phase < 0.0) phase += 2.0 * PI;

        cosfit[j] = { names[j], amplitude, phase };
    }
}


void Fit::writeFile(double L, std::string file)
{
	for (int j = 0; j < cosfit.size(); j++) {
		std::cout << std::fixed << std::setprecision(10)
			<< "\n" << cosfit[j].name
			<< ",	 Amp = " << cosfit[j].amplitude
			<< ",	 Phase = " << cosfit[j].phase << " rad\n";
	}

	std::ofstream out2(file);
	out2 << cosfit[0].name << "," << cosfit[0].amplitude << "," << cosfit[0].phase << std::endl;
	out2 << cosfit[1].name << "," << cosfit[1].amplitude << "," << cosfit[1].phase << std::endl;
	double t;
	for (int i = 0; i < (*matrix).rows(); ++i)
	{
		t = (i + 1) * dt;
		out2 << (i + 1) * dt * sqrt(G / L) << "," 
			<< cosfit[0].amplitude * cos(omiga * t + cosfit[0].phase) 
			<< "," << cosfit[1].amplitude * cos(omiga * t + cosfit[1].phase) 
			<< std::endl;
	}
}

