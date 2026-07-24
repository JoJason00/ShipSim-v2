#pragma once

#include "Eigen/Dense"
#include <memory>
#include <string>
#include <vector>
#include "../config/CaseConfig.h"
#include "../io/CaseLoader.h"
#include "Element.h"
#include "../wave/WaveBase.h"
#include "../wave/WaveField.h"
#include "Gsinteg.h"
#include "../wave/RegularWave.h"
#include "RAO4DTable.h"
#include "SeakeepingDOF.h"
#include "GreenTable.h"


#ifndef GREEN_STORE_FLOAT
#define GREEN_STORE_FLOAT 1
#endif

#if GREEN_STORE_FLOAT
using GScalar = float;
#else
using GScalar = double;
#endif

using GMat = Eigen::Matrix<GScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
using GVec = Eigen::Matrix<GScalar, Eigen::Dynamic, 1>;
using SgMatG = Eigen::Matrix<GScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>; // NE x bufCols

using ODEFunction = std::function<std::vector<double>
    (double t, const std::vector<double>& y, double alpha, const Eigen::MatrixXd& totalForce)>;

struct CaseContext {
    int i_case = 0;
    int DOF = 0;
    int iFn = 0;
    int iDir = 0;
    int iw = 0;

    std::shared_ptr<WaveBase> wave;
    std::shared_ptr<RegularWave> reg;

    // Unified incident-sea description. Regular = 1 component; irregular = N
    // at one direction; cross = the two systems live in subWaves/subFields.
    WaveField waveField;

    // Step 4: a crossing sea is the linear sum of independent sub-waves, each
    // with its own FK kernel (direction/region) and its own incident-elevation
    // field. Regular/Irregular -> exactly one entry (== wave / waveField), so
    // the single-wave path is bit-identical. Cross -> two entries; the per-step
    // excitation is Σ_s ( kernel_s ⊛ η_s ), accumulated into ExcitingForce.
    std::vector<std::shared_ptr<WaveBase>> subWaves;
    std::vector<WaveField>                 subFields;
    int activeSub = 0;   // which sub-wave incidentEtaAtShip / kernel use now

    // Ship start position [m] (Earth frame) for a crossing sea. Adds the
    // spatial phase −k_s(x0 cosθ_s + y0 sinθ_s) per sub-wave component, so the
    // two crossing systems' relative phase at the ship depends on where it
    // starts. 0 for single regular/irregular waves -> path stays bit-identical.
    double crossStartX = 0.0;
    double crossStartY = 0.0;

    // Per-wave-entry debug switch: dump incident η(t) at the fixed body
    // reference point to wave_history/. Set from WaveBase::outputHistory.
    bool outputWaveHistory = false;

    // Step 5: following / stern-quartering seas. ω_e(ω) is non-monotonic, so a
    // broadband sub-wave is split by encounter region. -1 = no filter (head /
    // beam: ω-independent kernel, single pass). 0..3 = only sum components
    // whose wave_force_region::classify == this WaveForceRegion.
    int activeRegion = -1;

    double Amp = 0.0;
    double W = 0.0;
    double we = 0.0;
    double dirRad = 0.0;
    double Fn = 0.0;
    double U = 0.0;
    double UsquareG = 0.0;
    double dt = 0.0;

    bool new_direction = true;
    bool new_Fn = true;

    //used to genarate Impulse function
    double dt_const = 0.0;
    double TG_const = 0.0;
    double tMot_const = 0.0;

    int TS = 0;
    int tMot = 0;
    int TG = 0;
    int NE = 0;
    int n_WL = 0;

    std::array<double,6> a_scale;          // rho*disp*G*Amp/L
    std::string tag;

    Eigen::MatrixXd ExcitingForce; // TS x DOF, total first-order wave force
    Eigen::MatrixXd FKForceHist;   
    Eigen::MatrixXd DiffForceHist; 

    double F3_amplitude = 0;
    double F5_amplitude = 0;
};


struct GreenfData {
    std::unique_ptr<std::vector<Eigen::MatrixXd>> Gz;
    std::unique_ptr<std::vector<Eigen::MatrixXd>> dGz;
    std::unique_ptr<std::vector<Eigen::MatrixXd>> Gz_wl;
    std::unique_ptr<std::vector<Eigen::MatrixXd>> dGz_wl;
};

//   䲨          
struct eForce
{
    Eigen::MatrixXd veSg_d;      // ֲ Դǿ,double  
    SgMatG veSg_g;               //GScalar        ʷ         ÿ   cast  

    Eigen::MatrixXd fkForce;    //FK  ʱ         ڵĲ   
    Eigen::MatrixXd dForce;     //      ʱ         ڵĲ   
    Eigen::VectorXd sPot;

    Eigen::MatrixXd forceIntHist;   // I_k(t) = -ArInt_k · phi
    Eigen::MatrixXd forceIntSmoothHist;  // smoothed I_k(t)
    Eigen::MatrixXd forceDtHist;    // dI_k/dt after smoothed differentiation

    Eigen::RowVectorXd eFdt;
    Eigen::RowVectorXd force;   //FK  +     

    Eigen::MatrixXd veVn_d;
    SgMatG         veVn_g;
    Eigen::MatrixXd vePhi_d;
    SgMatG         vePhi_g;

    eForce() = default;

    void resize(int NE, int TS, int TG)
    {
        // If geometry/time-grid shape is unchanged, keep FK / diffraction impulse rows.
        // LinearCumminsTDGF may skip prepareFKImpulse when only ω_inc changes inside the
        // same Lu bucket — clearing here would zero the wave-force convolution.
        if (veSg_d.rows() == NE && veSg_d.cols() == TG + 1
            && fkForce.rows() == TS && fkForce.cols() == 6
            && dForce.rows() == TS && dForce.cols() == 6)
        {
            return;
        }

        veSg_d = Eigen::MatrixXd::Zero(NE, TG + 1);
        veSg_g = SgMatG::Zero(NE, TG + 1);
        fkForce = Eigen::MatrixXd::Zero(TS, 6);
        dForce = Eigen::MatrixXd::Zero(TS, 6);
        sPot = Eigen::VectorXd::Zero(NE);
        eFdt = Eigen::RowVectorXd::Zero(6);
        force = Eigen::RowVectorXd::Zero(6);

        forceIntHist.resize(TS, 6);
        forceIntHist.setZero();

        forceIntSmoothHist.resize(TS, 6);
        forceIntSmoothHist.setZero();

        forceDtHist.resize(TS, 6);
        forceDtHist.setZero();

        veVn_d = Eigen::MatrixXd::Zero(NE, TG + 1);
        veVn_g = SgMatG::Zero(NE, TG + 1);
        vePhi_d = Eigen::MatrixXd::Zero(NE, TG + 1);
        vePhi_g = SgMatG::Zero(NE, TG + 1);
    }
};

//      
struct rForce
{
    Eigen::MatrixXd vSg_d;      // ֲ Դǿ,double  
    SgMatG vSg_g;               //GScalar        ʷ         ÿ   cast  

    Eigen::MatrixXd rdForce;    //FK  ʱ         ڵĲ   
    Eigen::VectorXd sPot;

    Eigen::RowVectorXd rFdt;
    Eigen::RowVectorXd force;   //      

    rForce() = default;

    void resize(int NE, int TS, int TG)
    {
        vSg_d = Eigen::MatrixXd::Zero(NE, TG + 1);
        vSg_g = SgMatG::Zero(NE, TG + 1);
        rdForce = Eigen::MatrixXd::Zero(TS, 6);
        sPot = Eigen::VectorXd::Zero(NE);
        rFdt = Eigen::RowVectorXd::Zero(6);
        force = Eigen::RowVectorXd::Zero(6);
    }
};

struct addedData
{
    double a33, a43, a53, a34, a44, a54, a35, a45, a55;
    double b33, b43, b53, b34, b44, b54, b35, b45, b55;

    addedData()
        :a33(0), a43(0), a53(0), a34(0), a44(0),
        a54(0), a35(0), a45(0), a55(0),
        b33(0), b43(0), b53(0), b34(0), b44(0),
        b54(0), b35(0), b45(0), b55(0) {
    };
};


//struct RollViscDamping {
//    double B44_lin = 1.3; //     ճ       (N*m*s/rad)
//    double B44_quad = 0.0; //         ϵ   (N*m*s^2/rad^2)   ѡ
//    double B44_cube = 0.0;
//};


class Seakeeping
{
public:
    explicit Seakeeping(const ShipConfig& ShipCfg, std::string casePath, const SeakeepingConfig& SeakeepingCfg);

    void run();

    //GreenfData moveGreenf();

private:
    void ProcessElement();
    void solve();

    void adaptTimeByCircle(double omega);

    RAO4DTable  initRAO4D();
    void initGreenTable();
    CaseContext buildCaseContext(const int i_case, const double Fn, const int iFn);
    void allocCaseBuffers(const CaseContext& ctx);                         //      Gz/forces/motions   
    void computeGreenTables(const CaseContext& ctx);                         // OpenMP    
    void computeExciting(CaseContext& ctx);                               // ʱ     &    ExcitingForce
    void integrateAndFitStore(const CaseContext& ctx, RAO4DTable& tab);        // rk4 + fit +   RAO

    static int findExactIndex(const std::vector<double>& xs, double x, double tol);

    void GreenFunction(int tN, double tn, Gsinteg& green);
    void initialFK(fkpData& fkpdata);
    void ExcitingCal(int tN, double tn, FKphi& fkphi, const CaseContext& ctx);
    void writeFKImplese(const std::string filePath, const CaseContext& ctx);
    void addedBoundary(double tn, int mode, double Amp, double we);
    void RadiationCal(int tN, double tn, const Eigen::VectorXd rVn);
    void AddedMassAndDamping(int mode, Eigen::MatrixXd& Force, int n, double we, double Amp, double& addedMass, double& dampingCoefficient);

    //void SourceConvolution(int&, const Eigen::VectorXd&, Eigen::MatrixXd&, Eigen::VectorXd& sPot);
    void ForceCal(int& tN, Eigen::VectorXd& Pt, Eigen::RowVectorXd& Fdt, Eigen::RowVectorXd& F);
    //void SourcePotential  (int t0, int tN, Eigen::VectorXd& Pt, Eigen::MatrixXd& Sg);
    void SourceConvolution(int& tN, const Eigen::VectorXd& Vn,
        Eigen::MatrixXd& Sg_d, SgMatG& Sg_g,
        Eigen::VectorXd& sPot);
    void SourcePotential(int t0, int tN, Eigen::VectorXd& Pt,
        const Eigen::MatrixXd& Sg_d, const SgMatG& Sg_g);

    void rk4_solve(
        const ODEFunction& f, const Eigen::MatrixXd& ExcitingForce,
        const std::vector<double>& y0, double t0, double t_end, double h);

    std::vector<double> rk4_step(
        const ODEFunction& f, const Eigen::MatrixXd& ExcitingForce, int tN, double t0,
        const std::vector<double>& y0, double h,
        std::vector<double>& k1, std::vector<double>& k2,
        std::vector<double>& k3, std::vector<double>& k4,
        std::vector<double>& y_temp, Eigen::MatrixXd& totalForce, const RollViscDamping& vd);

    void runFreeRollDecay();

private:
    std::string         filePath;
    ShipConfig 	        ShipCfg;
    SeakeepingConfig    SeakeepingCfg;

    std::shared_ptr<Element> element;

    RollViscDamping rollVisc_;

    std::vector<double> betaAxisDeg;
    std::vector<double> omegaAxis;

    double U = 0.0;
    double UsquareG = 0.0;

    HydrostaticsData    hs;
    addedData           added;
    GreenTable          gGreenTable;

    GVec N0sq;

    Eigen::VectorXi PotL_idx;

    Eigen::VectorXd ak;             //   Է      Ҷ   
    Eigen::VectorXd rVn;            //       ߽     

    eForce eforce;
    rForce rforce;

    Eigen::MatrixXd motions;

    std::unique_ptr<std::vector<GMat>> Gz{ nullptr };
    std::unique_ptr<std::vector<GMat>> dGz{ nullptr };
    std::unique_ptr<std::vector<GMat>> Gz_wlT{ nullptr };
    std::unique_ptr<std::vector<GMat>> dGz_wlT{ nullptr };

    int bufCols = 0;
    inline int col(int tn) { return tn % bufCols; }
};
