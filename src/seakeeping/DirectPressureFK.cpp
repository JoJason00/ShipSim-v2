#include "DirectPressureFK.h"
#include "../const/Const.h"

#include <cmath>
#include <stdexcept>

Eigen::RowVectorXd DirectPressureFK::computeForce6(
    const Element& element,
    const DirectPressureFKContext& ctx,
    double /*t*/)   // time now folded into ctx.phiOrigin
{
    // Calm water / no incident wave (e.g. coupled run with empty wave list): zero FK.
    if (ctx.amp <= 0.0)
        return Eigen::RowVectorXd::Zero(6);

    if (ctx.omega <= 0.0)
        throw std::runtime_error("DirectPressureFK: incident omega must be positive.");

    const int NE = static_cast<int>(element.xcr.size());
    if (element.ArInt.rows() != NE || element.ArInt.cols() != 6)
        throw std::runtime_error("DirectPressureFK: invalid ArInt size.");

    const double k = ctx.omega * ctx.omega / G;
    const double cb = std::cos(ctx.beta);
    const double sb = std::sin(ctx.beta);

    constexpr double invSqrt3 = 0.5773502691896258;
    const double gaussPts[2] = { invSqrt3, -invSqrt3 };
    Eigen::RowVectorXd force6 = Eigen::RowVectorXd::Zero(6);

    for (int i = 0; i < NE; ++i)
    {
        const ElementMatrix& corners = element.ElementData->at(i);
        const Eigen::RowVectorXd arint = element.ArInt.row(i) / element.Area[i];
        const Eigen::Vector3d c0 = corners.row(0).transpose();
        const Eigen::Vector3d c1 = corners.row(1).transpose();
        const Eigen::Vector3d c2 = corners.row(2).transpose();
        const Eigen::Vector3d c3 = corners.row(3).transpose();

        for (int gi = 0; gi < 2; ++gi)
        {
            const double xi = gaussPts[gi];
            const double one_minus_xi = 1.0 - xi;
            const double one_plus_xi = 1.0 + xi;

            for (int gj = 0; gj < 2; ++gj)
            {
                const double eta = gaussPts[gj];
                const double one_minus_eta = 1.0 - eta;
                const double one_plus_eta = 1.0 + eta;

                const double N1 = 0.25 * one_minus_xi * one_minus_eta;
                const double N2 = 0.25 * one_plus_xi * one_minus_eta;
                const double N3 = 0.25 * one_plus_xi * one_plus_eta;
                const double N4 = 0.25 * one_minus_xi * one_plus_eta;

                const Eigen::Vector3d dx_dxi =
                    (-0.25 * one_minus_eta) * c0
                    + (0.25 * one_minus_eta) * c1
                    + (0.25 * one_plus_eta) * c2
                    + (-0.25 * one_plus_eta) * c3;

                const Eigen::Vector3d dx_deta =
                    (-0.25 * one_minus_xi) * c0
                    + (-0.25 * one_plus_xi) * c1
                    + (0.25 * one_plus_xi) * c2
                    + (0.25 * one_minus_xi) * c3;

                const double jacobian = dx_dxi.cross(dx_deta).norm();
                const double x = N1 * c0[0] + N2 * c1[0] + N3 * c2[0] + N4 * c3[0];
                const double y = N1 * c0[1] + N2 * c1[1] + N3 * c2[1] + N4 * c3[1];
                const double z = N1 * c0[2] + N2 * c1[2] + N3 * c2[2] + N4 * c3[2];

                // Φ0 already carries ω t and the trajectory term −k(xe cosθ +
                // ye sinθ) + ε; only the body-frame spatial term remains.
                const double phase = ctx.phiOrigin - k * (x * cb + y * sb);
                const double p = rho * G * ctx.amp * std::exp(k * z) * std::cos(phase);
                force6 -= arint * (p * jacobian);
            }
        }
    }

    return force6;
}


//#include "DirectPressureFK.h"
//#include "../const/Const.h"
//
//#include <cmath>
//#include <stdexcept>
//
//Eigen::RowVectorXd DirectPressureFK::computeForce6(
//    const Element& element,
//    const DirectPressureFKContext& ctx,
//    double t)
//{
//    if (ctx.amp <= 0.0)
//        throw std::runtime_error("DirectPressureFK: wave amplitude must be positive.");
//    if (ctx.omega <= 0.0)
//        throw std::runtime_error("DirectPressureFK: incident omega must be positive.");
//
//    const int NE = static_cast<int>(element.xcr.size());
//    if (element.ArInt.rows() != NE || element.ArInt.cols() != 6)
//        throw std::runtime_error("DirectPressureFK: invalid ArInt size.");
//
//    const double k = ctx.omega * ctx.omega / G;
//    const double cb = std::cos(ctx.beta);
//    const double sb = std::sin(ctx.beta);
//
//    Eigen::VectorXd pressure(NE);
//
//    for (int i = 0; i < NE; ++i)
//    {
//        const double phase =
//            ctx.omega_e * t
//            - k * (element.xcr[i] * cb + element.ycr[i] * sb);
//
//        pressure[i] = rho * G * ctx.amp * std::exp(k * element.zcr[i])
//            * std::cos(phase);
//    }
//
//    // ArInt(i,j) = Area_i * n_j_i
//    // ������ѹ����Ϊ -p*n����������ȡ����
//    Eigen::RowVectorXd force6 = -pressure.transpose() * element.ArInt;
//
//    return force6;
//}