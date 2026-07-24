#pragma once

#include <string>
#include <vector>
#include <Eigen/Dense>

struct cosFit {
    std::string name;
    double amplitude;
    double phase;
};

class Fit
{
public:
    Fit();
    void setdata(const Eigen::MatrixXd& data, double dt, double we, double Amp, std::vector<std::string> name);
    void run();
    void writeFile(double L, std::string file);
    const std::vector<cosFit>& results() const { return cosfit; }

private:
    const Eigen::MatrixXd* matrix{ nullptr };
    double dt, omiga, Amp;
    std::vector<std::string> names;
    std::vector<cosFit> cosfit;

    int cut(const Eigen::VectorXd& signal, double dt, double omiga);
};

