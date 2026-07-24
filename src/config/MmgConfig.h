#pragma once
#include "../src/wave/WaveBase.h"
#include <string>

struct MmgTimeConfig
{
	double dt{};
	int SeakeepingStep{};
	double v_interavl{};
	double angle_interval{};
};

struct hullCoeff
{
	double Xvv{};
	double Xvr{};
	double Xrr{};
	double Xvvvv{};
	double Xuu{};

	double Yv{};
	double Yr{};
	double Yvvv{};
	double Yvvr{};
	double Yvrr{};
	double Yrrr{};

	double Nv{};
	double Nr{};
	double Nvvv{};
	double Nvvr{};
	double Nvrr{};
	double Nrrr{};
};

struct propellerCoeff
{
	double tp{};
	double np{};
	double Dp{};

	double k0{};
	double k1{};
	double k2{};

	double C1{};
	double C2{};

	double xp{};

	//propeller wake fraction factor
	double wP0{};
};

struct rudderCoeff
{
	double omigaR{};
	double HR{};
	double kappa{};
	double lR{};
	double AR{};
	double f_alpha{};
	double tR{};
	double aH{};
	double xR{};
	double xH{};
	double epsilonP{};

	double kP{};
	double gammaR_plus{};
	double gammaR_minus{};
};

struct addedCoeff
{
	double mx;
	double my;
	double Jz;
};

struct TuringCfg
{
	double TurningCaseCicle;
	double rudder_angle;
	double rudder_rate;

	// 停止判据：
	//   "circle" — 累计航向转够 TurningCaseCicle 圈后停（默认，原行为）
	//   "t_star" — 无量纲时间 t* = t·U0/Lpp 到达 t_star_max 后停
	std::string stopBy{ "circle" };
	double t_star_max{ 0.0 };
};

struct ZigzagCfg
{
	double ZigzagCaseCicle;
	double rudder_angle;
	double rudder_rate;
};

/// Optional parameters for `CoupledMmg3DOFCore` (slow-time 3-DOF MMG in coupled runs).
/// When `enabled` is true, geometry / mass / propeller rpm / rudder rate are read from
/// `case.json` under `Manoeuvring.CoupledMmg3DOF` instead of the legacy hard-coded defaults
/// in `CoupledMmg3DOFCore.cpp`.
struct CoupledMmg3DOFSettings
{
	bool enabled{ false };
	/// When true, hull derivatives and propeller wake/KT polynomial (nondimensional) are
	/// taken from the existing `hullCoeff` / `propellerCoeff` blocks already parsed into `MmgConfig`.
	bool reuseManoeuvringNondimHullProp{ true };

	double rho{ 1025.0 };
	double g{ 9.81 };
	double Lpp{};
	double Lwl{};
	double B{};
	double d{};
	double Cb{ 0.572 };
	double massKg{};
	double xG_m{};
	double kzz_nd{ 0.269 };

	double Dp{};
	double HR{};
	double BR{};
	double AR{};

	double xP_nd{};
	double xR_nd{ -0.5 };
	double xH_nd{ -0.48 };
	double lR_nd{ -1.0 };

	double epsilon{ 0.921 };
	double kappa{ 0.631 };
	double tR{ 0.29 };
	double aH{ 0.237 };
	double gammaR_plus{ 0.139 };
	double gammaR_minus{ 0.088 };

	double deltaRateDegPerS{ 12.0 };
	double maxDeltaDeg{ 35.0 };

	double Fn{ 0.15 };
	double npRps{};

	double mx_nd{ 0.0044 };
	double my_nd{ 0.1299 };
	double Jz_nd{ 0.0077 };

	// ---- Optional 4-DOF (add roll) -------------------------------------------
	// maneuverDOF = 3 (default, unchanged behaviour) or 4 (integrate the slow
	// heel via the roll-coupled MMG roll equation). Roll does NOT feed back into
	// surge/sway/yaw (hull derivatives are f(beta,r') only), so the 3-DOF
	// trajectory is bit-identical whether 3 or 4 DOF is selected.
	int maneuverDOF{ 3 };
	// Added roll moment of inertia Jxx [kg*m^2]. Total roll inertia = Ixx (ship
	// Mass.Ixx) + rollJxxAdd. Set ~ seakeeping A_inf(roll,roll) so the slow heel
	// and the wave-frequency roll share the same natural frequency.
	double rollJxxAdd{ 0.0 };
	// Manoeuvring roll damping (same form as seakeeping: -(lin*p + quad*|p|*p)).
	// Default to the seakeeping-calibrated values so the heel settles realistically.
	double rollB44_lin{ 1.0 };
	double rollB44_quad{ 8.0 };
	// Roll moment arms (non-dimensional ratios, scale-independent):
	//   zH = zH_over_d * d  (hull lateral force depth -> centrifugal heel)
	//   zR = zR_over_d * d  (rudder force depth)
	//   zG = zG_over_d * d  (CG height -> rudder roll moment + inertial coupling)
	// With zG_over_d = 0 the rudder roll moment and the inertial coupling vanish,
	// so only the centrifugal (zH*YH) heel is produced.
	double zH_over_d{ 0.5 };
	double zR_over_d{ 0.7 };
	double zG_over_d{ 0.0 };
};

/// Optional 3-DOF MMG variant aligned with the SJTU thesis (Son–Nomoto) formulation:
/// hull forces via \(V=\sqrt{u^2+v^2}\), \(v'=v/V\), \(r'=rL/V\); \(J=u(1-w_p)/(nD)\);
/// rudder \(v_R\) with \(C_{Rv},C_{Rrv},C_{Rrrv}\); ITTC friction + \(C_r,k_f\) residuary resistance.
/// JSON: `Manoeuvring.CoupledSjtuMmg3DOF`. Select with `Manoeuvring.lowFreqManeuverModel` = `sjtu_s175_nomoto`.
struct CoupledSjtuMmg3DOFSettings
{
	bool enabled{ false };

	double rho{ 1025.0 };
	double g{ 9.81 };
	double Lpp{};
	double Lwl{};
	double B{};
	double d{};
	double massKg{};
	double kzz_nd{ 0.269 };
	/// Wetted surface \(S\) (m²) for ITTC resistance.
	double S_wetted{};

	double Dp{};
	double HR{};
	double BR{};
	double AR{};

	double xP_nd{};
	double xR_nd{ -0.5 };
	double xH_nd{ -0.48 };
	double lR_nd{ -1.0 };

	/// Hull non-inertial forces \(X'_0,Y'_0,N'_0\) in the thesis sense (Eq. 2-10).
	double Xvv{};
	double Xvr{};
	double Xrr{};
	double Yv{};
	double Yr{};
	double Yvvv{};
	double Yvvr{};
	double Yvrr{};
	double Yrrr{};
	double Nv{};
	double Nr{};
	double Nvvv{};
	double Nvvr{};
	double Nvrr{};
	double Nrrr{};

	double mx_nd{};
	double my_nd{};
	double Jz_nd{};

	double Cr{ 0.001 };
	double kf{ 1.04 };
	/// Kinematic viscosity (m²/s) for ITTC \(C_f\).
	double nuWater{ 1.14e-6 };

	double tp{};
	double wP0{};
	double j0{};
	double j1{};
	double j2{};
	/// Wake \(w_p = w_{P0}\exp(-\mathrm{wakeExp}\,\beta_P^2)\) (thesis uses 8).
	double wakeExpBetaP{ 8.0 };

	double epsilon{ 0.921 };
	double kappa{ 0.631 };
	double tR{ 0.29 };
	double aH{ 0.237 };
	double gammaR_plus{ 0.193 };
	double gammaR_minus{ 0.088 };
	/// \(\tau\) in \(u_P = u[(1-w_p)+\tau\beta_P^2]\) (Eq. 2-21); 0 if omitted in literature table.
	double tauPropInflow{ 0.0 };
	double C_Rv{};
	double C_Rrv{};
	double C_Rrrv{};

	double deltaRateDegPerS{ 13.0 };
	double maxDeltaDeg{ 35.0 };
	double Fn{ 0.15 };
	double npRps{};
};

struct MmgConfig
{
	std::vector<std::shared_ptr<WaveBase>> waves;

	MmgTimeConfig Time;
	hullCoeff Hull;
	propellerCoeff Propeller;
	rudderCoeff Rudder;
	addedCoeff Added;

	std::vector<double> Fn;

	TuringCfg Turning;
	ZigzagCfg Zigzag;

	CoupledMmg3DOFSettings coupledMmg3DOF;
	CoupledSjtuMmg3DOFSettings sjtuMmg;

	/// Slow-time manoeuvring integrator: `mmg_default` → `CoupledMmg3DOFCore`;
	/// `sjtu_s175_nomoto` → `CoupledSjtuMmg3DOFCore` (requires `CoupledSjtuMmg3DOF.enabled`).
	std::string lowFreqManeuverModel{ "mmg_default" };

	bool hullCoeff_defined{ false };
	bool propellerCoeff_defined{ false };
	bool rudderCoeff_defined{ false };
	bool addedCoeff_defined{ false };
	bool turn_on_turningcase{ false };
	bool turn_on_zigzagcase{ false };
};