#pragma once

#include "Eigen/Dense"
#include <memory>
#include <string>
#include <vector>
#include <complex>

#include "../config/CaseConfig.h"
#include "../io/CaseLoader.h"
#include "Element.h"
#include "../wave/WaveBase.h"
#include "../wave/RegularWave.h"
#include "Gsinteg.h"
#include "RAO4DTable.h"
#include "SeakeepingDOF.h"
#include "GreenTable.h"
#include "RollDamping.h"
#include "Seakeeping.h"
#include "RadiationKernelCache.h"
#include "TDGFExcitationKernelIO.h"
#include "DirectPressureFK.h"
#include "TDGFProvider.h"
#include "../io/FKImpulseKernelIO.h"
#include "../tool/SavitzkyGolay.h"

class LinearCumminsTDGF
{
public:
    explicit LinearCumminsTDGF(
        const ShipConfig& shipCfg,
        std::string casePath,
        const SeakeepingConfig& seakeepingCfg);

    void run();

private:
    void assertRigidModesOnly() const;

    void processElement();
    void solve();

    RAO4DTable initRAO4D();
    void initGreenTable();

    void initTDGFProvider();

    void buildEncounterBandForFn(double Fn);
    void setupOfflineKernelGridForFn(double Fn);
    void setupOnlineGreenGrid(double dt, int TG);
    void allocSharedGreenBuffers();
    void writeGreenTable(int i_element);

    CaseContext buildCaseContext(int i_case, double Fn, int iFn) const;

    std::string makeKernelCachePath(double Fn) const;
    bool tryLoadKernelCache(double Fn, RadiationKernelData& kernel) const;
    void buildKernelCacheForFn(double Fn, RadiationKernelData& kernel);
    void writeKernelDebugCsv(double Fn, const RadiationKernelData& kernel) const;

    // ---- debug output ----
    void writeVectorHistoryDebugCsv(
        const std::string& quantityName,
        double Fn,
        int srcMode,
        const std::vector<Eigen::VectorXd>& history) const;

    void writeChiDebugCsv(
        double Fn,
        int srcMode,
        const std::vector<Eigen::VectorXd>& chi_history) const;

    void writeDchiDebugCsv(
        double Fn,
        int srcMode,
        const std::vector<Eigen::VectorXd>& dchi_history) const;

    // ---- post-processing of chi history ----
    std::vector<Eigen::VectorXd> smoothStartupHistory(
        const std::vector<Eigen::VectorXd>& raw_history,
        int fitSteps,
        int blendSteps,
        int polyDeg,
        double jaggedRatioThreshold = 1.35) const;

    std::vector<Eigen::VectorXd> differentiateHistory(
        const std::vector<Eigen::VectorXd>& history) const;

    std::vector<Eigen::MatrixXd> resampleKernelToDt(
        const RadiationKernelData& kernel,
        double dtOnline) const;

    void computeGreenTables();
    /// Recompute shared Green buffers if they do not match `TGShared_` / `dtShared_`
    /// (call sites: radiation kernel build, FK impulse direct build).
    void ensureGreenTablesMatchOnlineGrid();
    void GreenFunction(int tN, double tn, Gsinteg& green);

    void zeroWL(int tN);

    void initialFK(fkpData& fkpdata);

    void rebuildDiffractionForceFromIntegral(const CaseContext& ctx);

    void buildFKImpulseDirect(CaseContext& ctx, FKphi& fkphi);
    void prepareFKImpulse(CaseContext& ctx, FKphi& fkphi);

    void computeExciting(CaseContext& ctx);
    void ExcitingCal(int tN, double tn, FKphi& fkphi, const CaseContext& ctx);
    // dfOnly=true convolves only the diffraction (dForce) part of the impulse
    // table — the DirectPressure path uses this plus directFkOverComponents
    // for the FK part.
    void fkDfCal(int tN, double tn, FKphi& fkphi, const CaseContext& ctx,
                 bool dfOnly = false, bool fkOnly = false);

    // Direct panel-pressure FK summed over the active pass's incident
    // components (same region filter as incidentEtaAtShip), body-frame Φ0 =
    // we_j·t + ε_j (pure seakeeping = straight transit). Adds onto eforce.force.
    void directFkOverComponents(const CaseContext& ctx, double tNow);

    // Unified incident wave elevation seen by the (moving) ship at time t.
    // Step 1: delegates to ctx.wave->Eta(t) (zero behaviour change). Step 2
    // will switch to the encounter-frame sum over ctx.waveField components.
    double incidentEtaAtShip(const CaseContext& ctx, double t) const;
    void writeFKImplese(const std::string filePath, const CaseContext& ctx);



    std::string makeExcitationKernelCachePath(const CaseContext& ctx) const;
    void saveExcitationKernelCache(const CaseContext& ctx) const;

    void RadiationCal(int tN, double tn, const Eigen::VectorXd& rVnLocal);
    void ForceCal(int& tN, Eigen::VectorXd& Pt, Eigen::RowVectorXd& Fdt, Eigen::RowVectorXd& F, const double dt);

    Eigen::VectorXd buildPotentialRhs(
        const Eigen::VectorXd& boundary_condition) const;

    void faiConvolution(
        const int tN, 
        const CaseContext& ctx,
        const Eigen::VectorXd& Vn,
        Eigen::MatrixXd& VnHist_d,
        SgMatG& VnHist_g,
        Eigen::MatrixXd& PhiHist_d,
        SgMatG& PhiHist_g,
        Eigen::VectorXd& phiBody);

    void SourceConvolution(
        int& tN, const Eigen::VectorXd& Vn,
        Eigen::MatrixXd& Sg_d, SgMatG& Sg_g,
        Eigen::VectorXd& sPot);

    void SourcePotential(
        int t0, int tN, Eigen::VectorXd& Pt,
        const Eigen::MatrixXd& Sg_d, const SgMatG& Sg_g);

    // ---- 论文 3.3 / 3.4：刚体船舶运动 ----
    // Instantaneous 1/r source system from Eqs. (3.16a) / (3.17a,b)
    Eigen::VectorXd buildRankineRhs(const Eigen::VectorXd& boundary_condition) const;
    Eigen::VectorXd solveRankineBIE(const Eigen::VectorXd& boundary_condition) const;
    Eigen::VectorXd solveSourceLinearSystem(const Eigen::VectorXd& rhs) const;

    Eigen::VectorXd buildBoundaryCondition_a(int mode) const;
    Eigen::VectorXd buildBoundaryCondition_b(int mode) const;

    double weightedSurfaceIntegral(
        const Eigen::VectorXd& panelWeight,
        const Eigen::VectorXd& phi) const;

    double computeModeIntegral(int mode, const Eigen::VectorXd& phi) const;

    Eigen::VectorXd applyWaterlineBracket1(
        int lagIdx,
        const Eigen::VectorXd& panelValues) const;

    Eigen::VectorXd applyWaterlineBracket2(
        int lagIdx,
        const Eigen::VectorXd& panelValues) const;

    Eigen::VectorXd applyWaterlineN1SqG(
        int lagIdx,
        const Eigen::VectorXd& panelValues) const;

    Eigen::VectorXd applyWaterlineG(
        int lagIdx,
        const Eigen::VectorXd& panelValues) const;

    Eigen::VectorXd solveChiAtTime(
        int tN,
        const Eigen::VectorXd& bc_a,
        const Eigen::VectorXd& bc_b,
        const Eigen::VectorXd& psi1,
        const Eigen::VectorXd& psi2,
        const std::vector<Eigen::VectorXd>& chi_history) const;

    void buildKernelDirectTimeDomain(double Fn, RadiationKernelData& kernel);

    // 非均匀路径：仅在 useNonUniformChiGrid()=true 时被 buildKernelDirectTimeDomain dispatch。
    void buildKernelDirectTimeDomain_NonUniform(double Fn, RadiationKernelData& kernel);

    // 非均匀 chi 求解器：在 chi 节点 iCurr（≥1）处解 χ，按非均匀 left-rectangle 卷积。
    Eigen::VectorXd solveChiAtTime_NonUniform(
        int iCurr,
        const Eigen::VectorXd& bc_a,
        const Eigen::VectorXd& bc_b,
        const Eigen::VectorXd& psi1,
        const Eigen::VectorXd& psi2,
        const std::vector<Eigen::VectorXd>& chi_history) const;

    // 非均匀 3 点中心差分（端点用一阶单边差分）；与 differentiateHistory 同接口但接受 times。
    std::vector<Eigen::VectorXd> differentiateHistoryNonUniform(
        const std::vector<Eigen::VectorXd>& history,
        const std::vector<double>& times) const;

    // ---- 非均匀 chi 时间网格相关 ----
    // 把 SeakeepingCfg.ChiTimeGrid 展开成两套数组：
    //   chiIdxOnFine_[i] = 节点 i 在 dt_fine 细网格上的整数索引
    //   chiTimes_[i]     = chiIdxOnFine_[i] * dt_fine（物理时间）
    // 失败（disabled / 配置非法）时清空两数组。
    void buildChiTimeGridFromConfig();

    // 当前是否启用非均匀 chi 网格（= ChiTimeGrid.enabled 且展开成功）
    bool useNonUniformChiGrid() const
    {
        return !chiIdxOnFine_.empty();
    }

    void solveKernelCase(
        const CaseContext& ctx,
        const RadiationKernelData& kernel,
        RAO4DTable& tab);

    void solveKernelCase_RK4(
        const CaseContext& ctx,
        const RadiationKernelData& kernel,
        RAO4DTable& tab);

    static int findExactIndex(const std::vector<double>& xs, double x, double tol);
    static std::string keyDouble(double x);

    void runFreeRollDecay(const std::string& casePath);

    void TimeToFrequency(const std::string& casePat, double Fn, double beta, RadiationKernelData& kernel);

    void postprocessRadiationKernel(
        double Fn,
        RadiationKernelData& kernel) const;

    // 备选后处理：不做 tail-split（不动 A_inf/B/C_prime），改为把记忆核时间轴
    // 向后延长一倍，延长段补上 K(t) 末端“趋近的常数”。与 postprocessRadiationKernel
    // 二选一，在 solve() 里注释切换。详见 .cpp 实现处的说明。
    void postprocessRadiationKernel2(
        double Fn,
        RadiationKernelData& kernel) const;

    // Temporary, opt-in: replace kernel.Klag with an external K(t) CSV
    // (columns: t,K33,K35,K53,K55). Interpolated onto the existing kernel
    // time grid. A_inf/B/C' and the rest of the pipeline are untouched.
    void applyRadiationKernelOverride(RadiationKernelData& kernel) const;

    void writeKernelDiagnostics(
        double Fn,
        const RadiationKernelData& kernel,
        const std::string& stage = "raw",
        double dtForT = 0.0) const;

private:
    std::string      filePath;
    ShipConfig       ShipCfg;
    SeakeepingConfig SeakeepingCfg;

    std::string Solver;

    std::shared_ptr<Element> element;

    RollViscDamping rollVisc_;

    HydrostaticsData hs;
    GreenTable       gGreenTable;
    TDGFProvider     tdgfProvider_;

    std::vector<double> betaAxisDeg;
    std::vector<double> omegaAxisIncident;

    // Frequency-domain RAO4D table only makes sense for pure regular-wave
    // sweeps. Disabled when any wave is irregular/cross (skip build, fills
    // and CSV write — there is no single ω/β to tabulate).
    bool raoEnabled_ = true;

    GVec            N0sq;
    Eigen::VectorXi PotL_idx;

    // 水线几何
    Eigen::VectorXd wlTx_;
    Eigen::VectorXd wlTy_;
    Eigen::VectorXd wlN1_;
    Eigen::VectorXd wlN2_;
    Eigen::VectorXd wlSegLen_;
    Eigen::VectorXd wlDn1n2Dl_;
    std::vector<int> wlOrder_;

    Eigen::VectorXd ak;
    Eigen::VectorXd rVn;

    eForce eforce;
    rForce rforce;

    Eigen::MatrixXd motions;

    // 历史格林表：只存原始空间积分值，不预乘 dt
    std::shared_ptr<std::vector<GMat>> Green{ nullptr };
    std::shared_ptr<std::vector<GMat>> Green_dnP{ nullptr };
    std::shared_ptr<std::vector<GMat>> Green_dnS{ nullptr };
    std::shared_ptr<std::vector<GMat>> Gw{ nullptr };
    std::shared_ptr<std::vector<GMat>> Gw_dnP{ nullptr };
    std::shared_ptr<std::vector<GMat>> Gw_dnS{ nullptr };
    std::shared_ptr<std::vector<GMat>> Gw_dx{ nullptr };
    std::shared_ptr<std::vector<GMat>> Gw_dl{ nullptr };
    std::shared_ptr<std::vector<GMat>> Gw_dt{ nullptr };

    double U = 0.0;
    double UsquareG = 0.0;

    int bufCols = 0;
    inline int col(int tn) const { return tn % bufCols; }

    double dtShared_ = 0.0;
    int    TGShared_ = 0;
    int    TSShared_ = 0;

    // Green buffers were built for this (dtShared_, TGShared_); must rebuild when
    // setupOnlineGreenGrid changes them without new_direction/new_Fn (same heading, new ω).
    double greenBuiltForDt_{ -1.0 };

    double cfgTimeDtInput_ = 0.0;
    int    cfgGreenStepInput_ = 0;

    // offline kernel grid
    double kernelDtFixed_ = 0.0;
    int    kernelTGFixed_ = 0;

    double omegaMin;
    double omegaMax;

    std::vector<double> F3_amp = { 10.5, 8.19, 7.30, 6.38 };

    //std::vector<double> F5_amp = { 2.18, 2.18, 2.20, 2.18 };

    std::vector<double> F5_amp = { 2.5, 2.55, 2.55, 2.55 };

    // ---- 非均匀 chi 网格展开结果（buildChiTimeGridFromConfig 写入）----
    std::vector<int>    chiIdxOnFine_;   // 节点 i 在 dt_fine 细网格上的整数索引
    std::vector<double> chiTimes_;       // = chiIdxOnFine_[i] * dt_fine
    double              chiDtFine_       = 0.0;
    int                 chiMemoryCutoffLag_ = 0;   // 0/负 => 不截断

    // FK impulse: one library table per (Fn, Dir, Lu bucket H or 1/2/3), not per ω_inc.
    bool fkImpulsePrepared_ = false;
    double fkImpulsePrepFn_ = 0.0;
    double fkImpulsePrepDirRad_ = 0.0;
    WaveForceRegion fkImpulsePrepRegion_ = WaveForceRegion::Head;
    double fkImpulsePrepDtConst_{ -1.0 };
    int    fkImpulsePrepTMot_{ -1 };
};