#pragma once
#include <Eigen/Dense>
#include <optional>

struct FKphi
{
    Eigen::VectorXd df;
    Eigen::VectorXd fk;

    FKphi(int size)
        : df(Eigen::VectorXd::Zero(size)),
          fk(Eigen::VectorXd::Zero(size)) {
    }
};

struct fkpData
{
    int    NE;
    double U;
    Eigen::VectorXd A31, A32, A33;
    Eigen::VectorXd xcr, ycr, zcr;

    fkpData()
        :NE(0), U(0.0) {
    };
};



class WaveBase 
{
public:
    virtual ~WaveBase() = default;
    virtual double Eta(double t) const = 0;
    virtual void Exciting(double tn, FKphi& fkphi) = 0;

    virtual void loadData(const fkpData& Data) = 0;
    virtual double getAmp() = 0;
    virtual double getFreq() = 0;

    // Absolute propagation direction [rad]. Default 0 keeps existing call
    // sites valid; every concrete wave overrides it so LinearCumminsTDGF
    // loops can stay wave-type agnostic.
    virtual double direction() { return 0.0; }

    // Per-wave-entry debug switch ("output_history" in the wave JSON, same
    // level as "type"). When true the solver dumps the incident wave-elevation
    // time history at the fixed body reference point (the η used in the
    // wave-force integration) to a wave_history/ folder. Off by default ->
    // no behaviour or output change.
    bool outputHistory = false;

protected:
    std::optional<fkpData> data;
};

