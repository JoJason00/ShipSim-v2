#pragma once
#include "source_dipole.h"
#include <Eigen/Dense>

using ElementMatrix = Eigen::Matrix<double, 4, 3, Eigen::RowMajor>;
using Eigen::Vector2d, Eigen::Vector3d, Eigen::Vector4d;

class Element
{

public:
	explicit Element(std::string Method, int NumElem,
		std::unique_ptr<std::vector<ElementMatrix>> data);

	void Geometry(const Eigen::Vector3d cg);
	void RankineSource();
	void RankineSource2();

	int waterline();

private:
	void Cornerp(Vector2d&, Eigen::Matrix<double, 4, 2, Eigen::RowMajor>&);
	void NullPoints(double& xn, double& yn, Vector4d&, Vector4d&);
	void NullCentroid(int&, double&, double&, Vector2d&, Vector3d&, Vector3d&,
		Eigen::Matrix<double, 3, 2, Eigen::ColMajor>&);

	void TransforMatrix(int&, Vector3d&, Vector3d&, Vector3d&);
	void MomentInertia(int&, Vector4d&, Vector4d&);
	void LineValues(int&, int&, ElementMatrix&);

	void Rankine(int j, int i);
	void NullCentroidDis(int j, int i, double& dx, double& dy, double& dz);
	void MonopoleCal(int j, int i, double r0);
	void DepoleCal(int j, int i, double r0, double& Vxe, double& Vye, double& Vze);
	void ElementCoor(int j, int i, double& xne, double& yne, double& zne);
	void ExactCal(int j, int i, double& Vxe, double& Vye, double& Vze);
	void ReferenceCoor(int i, double& Vxe, double& Vye, double& Vze);

	void ExactDist(int& i, double& xne, double& yne, double& zne, double* de,
		double* se, double* ce, double* re, double* R);

	void ExactDis(int& i, double& xne, double& yne, double* ce, double* se,
		std::array<std::array<double, 2>, 4>& se1);

	void compute_Sz(int j, int i, double xne, double yne, double zne, double* c_Sz,
		double* h_Sz, double* m_Sz, double* r_Sz);

	void ExactSolu(double* de, double* R, double& dt1);
	void eNormalVelocity(int i, double& Rn);
	void SignChangez(int i);

	void PointtoMesh(const ElementMatrix& point, const int idex);

	void saveRankine(const std::string filePath);

private:
	const int NE;
	bool method;
	std::unique_ptr<std::vector<ElementMatrix>> ElementData;

	AllMesh meshes;

	Eigen::MatrixXd Nvec, xpl, ypl, xpe, ype;
	Eigen::VectorXd Area, xnr, ynr, znr, xcr, ycr, zcr;
	Eigen::VectorXd A11, A12, A13, A21, A22, A23, A31, A32, A33;
	Eigen::VectorXd ixx, ixy, iyy, td;
	Eigen::VectorXi PotL;

	Eigen::MatrixXd ArInt;

	Eigen::MatrixXd Rz;
	Eigen::MatrixXd A_rankine;
	//Eigen::MatrixXd Inv;
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;	//LU ֽ       Է     

	double Vxr, Vyr, Vzr, Ph;

	double SA = 0;	    //    ʪ     
	int	   n_WL = 0;		//ˮ    Ԫ  

	friend class Gsinteg;
	friend class Seakeeping;
	friend class SeakeepingDOF;
	friend class LinearCumminsTDGF;
	friend class DirectPressureFK;
};