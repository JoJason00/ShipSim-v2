#include "Element.h"
#include "../const/Const.h"
#include "source_dipole.h"
#include <iostream>
#include <fstream>

Element::Element(std::string Method, int ne, std::unique_ptr<std::vector<ElementMatrix>> data)
	:NE(ne), ElementData(move(data)), Nvec(ne, 6), xpl(ne, 2), ypl(ne, 2), xpe(ne, 4), ype(ne, 4),
	Area(ne), ArInt(ne, 6), xnr(ne), ynr(ne), znr(ne), xcr(ne), ycr(ne), zcr(ne),
	A11(ne), A12(ne), A13(ne), A21(ne), A22(ne), A23(ne), A31(ne), A32(ne), A33(ne),
	ixx(ne), ixy(ne), iyy(ne), td(ne), PotL(ne), Rz(ne, ne), A_rankine(ne, ne),
	Vxr(0), Vyr(0), Vzr(0), Ph(0), meshes(ne)
{
	if (Method == "Source")
		method = true;
	else if (Method == "Potential")
		method = false;
	else
		throw std::runtime_error("unkonw method: " + Method + "\n");
}


void Element::Geometry(const Eigen::Vector3d cg)
{
	double xn, yn;

	ElementMatrix point;		//四个角点
	Eigen::Vector3d avpoint;		                        //中心点
	Eigen::Vector3d diag1, diag2;		                    //对角线
	Eigen::Vector3d nv, tv, sv;		                        //单位法向量和切向量
	Eigen::Vector3d nv1;
	Eigen::Vector2d cpe;
	Eigen::Vector2d cp;
	Eigen::Vector3d cnt, null;

	Eigen::Matrix<double, 4, 2, Eigen::RowMajor> locp;      //局部坐标下的角点
	Eigen::Matrix<double, 4, 2, Eigen::RowMajor> ce;
	Eigen::Matrix<double, 2, 3, Eigen::RowMajor> tme;
	Eigen::Matrix<double, 3, 2, Eigen::ColMajor> tmeT;

	Eigen::Vector4d xe, ye;

	int k, j;
	for (k = 0, n_WL = 0; k < NE; ++k, ++n_WL)
	{
		point = ElementData->at(k);

		//diag1 = point.row(2) - point.row(0);
		//diag2 = point.row(3) - point.row(1);

		//nv = (diag2.cross(diag1)).normalized();			//这里的法向量nv指向流体域

		 // 1. 先计算一次，检查法向量方向,点的顺序是对的，nV指向物面内，流体外
		avpoint = 0.25 * (point.row(0) + point.row(1) + point.row(2) + point.row(3));
		diag1 = point.row(2) - point.row(0);
		diag2 = point.row(3) - point.row(1);
		nv = (diag2.cross(diag1)).normalized(); 

		// --------------------------
		// 核心修复：如果方向不对，翻转顶点并写回数据
		// --------------------------
		if (nv.dot(cg - avpoint) > 0.0)
		{
			// 交换顶点 1 和 3 (老版本的逻辑)
			Eigen::RowVector3d tmp = point.row(1);
			point.row(1) = point.row(3);
			point.row(3) = tmp;

			ElementData->at(k) = point;

			// 重新计算对角线和法向量（现在方向应该对了）
			diag1 = point.row(2) - point.row(0);
			diag2 = point.row(3) - point.row(1);
			nv = (diag2.cross(diag1)).normalized();
		}

		//std::cout << "element("<<k<<"):\t"<<"nv0:\t" << nv(0)<<"\tnv1:\t" << nv(1) << "\tnv2:\t" << nv(2) << std::endl;

		PointtoMesh(point, k);

		tv = diag1.normalized();
		sv = nv.cross(tv);

		//if (!method) nv = -nv;		                    //势函数法nv指向物体内

		//avpoint = 0.25 * (point.row(0) + point.row(1) + point.row(2) + point.row(3));

		//local coordinates of corner points
		for (j = 0; j < 4; ++j)
		{
			locp.row(j)(0) = (point.row(j).transpose() - avpoint).dot(tv);
			locp.row(j)(1) = (point.row(j).transpose() - avpoint).dot(sv);
		}

		tme.row(0) = tv;
		tme.row(1) = sv;
		tmeT = tme.transpose();

		Cornerp(cpe, locp);
		cnt = tmeT * cpe + avpoint;

		for (j = 0; j < 4; ++j)
		{
			ce.row(j)(0) = locp.row(j)(0) - cpe(0);
			ce.row(j)(1) = locp.row(j)(1) - cpe(1);

			xe(j) = ce.row(j)(0);
			ye(j) = ce.row(j)(1);
		}

		NullPoints(xn, yn, xe, ye);
		NullCentroid(k, xn, yn, cp, null, cnt, tmeT);
		TransforMatrix(k, tv, sv, nv);

		//广义法向量
		//nv1 = null.cross(nv);
		nv1 = (cnt - cg).cross(nv);

		Nvec(k, 0) = nv[0];
		Nvec(k, 1) = nv[1];
		Nvec(k, 2) = nv[2];
		Nvec(k, 3) = nv1[0];
		Nvec(k, 4) = nv1[1];
		Nvec(k, 5) = nv1[2];

		xpe.row(k) = xe;
		ype.row(k) = ye;

		MomentInertia(k, xe, ye);
		LineValues(k, n_WL, point);
	}

	std::cout << "n_WL:\t" << n_WL << std::endl;

	for (int i = 0; i < 6; ++i)
		ArInt.col(i) = Area.cwiseProduct(Nvec.col(i));
}

void Element::RankineSource()
{
	int i, j;
	//double tmp,pz,dpz,Rn;
	double pz, dpz, Rn;

	//tmp=0.5/PI;
	Eigen::MatrixXd A(NE, NE);

	//源点
	for (j = 0; j < NE; j++) {
		//场点
		for (i = 0; i < NE; i++) {
			//------------------------
			//calculation of r values
			//-------------------------
			Rankine(j, i);


			//Rn:method 0在场点求法向速度，method 1在源点求法向速度
			if (!method)
				eNormalVelocity(i, Rn);
			else
				eNormalVelocity(j, Rn);
			// Rn = -Rn;

			pz = Ph;
			dpz = Rn;
			//-------------------------------
			//calculation of r prime values（镜像点）
			//-------------------------------
			//将场面元关于水平面镜像
			SignChangez(i);

			Rankine(j, i);

			pz -= Ph;
			if (!method)
				eNormalVelocity(i, Rn);
			else
				eNormalVelocity(j, Rn);
			// Rn = -Rn;


			dpz -= Rn;
			SignChangez(i);
			//------------------------------------------
			//store values
			//------------------------------------------		   
			//Rz[j][i]=pz*tmp;	dRz[j][i]=dpz*tmp;

			//Rz就是Rankine部分的结果，dRz是求导的结果
			Rz(j, i) = pz;
			A(j, i) = dpz;


			//if (!method)
			//{
			//	if (i == j)
			//		/*A(j, i) = 4.0 * PI - A(i, i);
			//		A(j, i) = 2.0 * PI;*/
			//		continue;
			//	else
			//		A(j, i) *= -1;
			//}

			//std::cout << "A(" << j << ")(" << i << ") = " << A(j, i) << std::endl;
		}
	}
	//------------------------------------------------
	//impulsive values for vertical plane
	//------------------------------------------------
	//矩阵求逆
	//Inv = Inv.inverse();

	//LU分解
	A_rankine = A;
	lu.compute(A_rankine);
}

void Element::RankineSource2()
{
	//场点
	std::array<double, 3> p_point;
	//源点
	std::array<double, 3> q_point;

	//面元单位法向量
	std::array<double, 3> n_vector;
	Source_Dipole sd{};

	//Rankine分布源及分布偶极
	double source_deltav11;
	double source_deltav12;
	double source_r1;
	double dipole_r1;

	//1/r-1/r1的积分
	double source_rr1;

	//-1/r1的导数的积分
	std::array<double, 3> source_v{};

	std::array<double, 3> source_deltav2{};
	std::array<double, 3> source_deltav2_2{};

	//均值
	std::array<double, 3> source_deltav2_avery;

	//每个时间步求解的Green函数
	//TDGF_coeff source_deltavG{};

	//六自由度广义法向量
	double nk = 0.0;

	//线性方程组的系数矩阵
	for (int i = 0; i < NE; i++)
	{
		p_point = meshes.centers[i];
		for (int j = 0; j < NE; j++)
		{
			q_point = meshes.centers[j];
			//面元法向量（源）
			n_vector = meshes.Nor_vectors[j];
			sd.update(p_point, meshes.Meshs.at(j));

			auto source_v = sd.source_v();

			source_deltav11 = sd.get_dipole();

			source_deltav12 = sd.get_source();


			//1/r1的rankine源积分也类似于1/r只是点的z坐标变为负即可
			p_point[2] = -p_point[2];
			sd.update(p_point, meshes.Meshs.at(j));
			source_r1 = -sd.get_source();
			dipole_r1 = -sd.get_dipole();
			p_point[2] = -p_point[2];

			//1/r-1/r1的积分
			source_rr1 = source_deltav12 + source_r1;

			Rz(i, j) = source_rr1;


			//线性方程组系数矩阵把六个方向的最后一列都放在后面
			/*for (int k = 0; k < 6; k++)
			{*/
			//int k = 2;
			//switch (k)
			//{
			//case 0: nk = n_vector[0];
			//	break;
			//case 1: nk = n_vector[1];
			//	break;
			//case 2: nk = n_vector[2];
			//	break;
			//case 3: nk = q_point[1] * n_vector[2] - q_point[2] * n_vector[1];
			//	break;
			//case 4: nk = q_point[2] * n_vector[0] - q_point[0] * n_vector[2];
			//	break;
			//case 5: nk = q_point[0] * n_vector[1] - q_point[1] * n_vector[0];
			//	break;
			//}
			//A_rankine[i][Num + k] += source_rr1 * nk;


			if (i != j)
				A_rankine(i,j) = source_deltav11 + dipole_r1;
			else
				A_rankine(i, j) = 2 * PI + dipole_r1;
		}
	}
	lu.compute(A_rankine);
}


int Element::waterline()
{
	return n_WL;
}


void Element::Cornerp(Vector2d& cp, Eigen::Matrix<double, 4, 2, Eigen::RowMajor>& cpe)
{
	double c1, c2, c3;

	c1 = cpe(1, 1) - cpe(3, 1);
	c2 = cpe(3, 0) * (cpe(0, 1) - cpe(1, 1));
	c3 = cpe(1, 0) * (cpe(3, 1) - cpe(0, 1));

	cp(0) = (c2 + c3) / (3 * c1);
	cp(1) = -cpe(0, 1) / 3;
}


void Element::NullPoints(double& xn, double& yn, Vector4d& xe, Vector4d& ye)
{
	int i, j, k;
	double D, rx, ry, Vx, Vx1, Vy, Vy1, Vxx, Vxy, Vyx, Vyy, detJ,
		Z1, Z2, xp1, yp1, r[4], d, vt;
	double vr, vr1, t1, t2;

	xp1 = 0; yp1 = 0;			   //初始点为（0，0）				
	for (j = 0; j < 20; j++) {	        //最多迭代20次
		for (i = 0; i < 4; i++) {
			t1 = xp1;
			t1 -= xe[i];
			t1 *= t1;
			t2 = yp1;
			t2 -= ye[i];
			t2 *= t2;
			t2 += t1;
			r[i] = sqrt(t2);		//当前点到四个角点的距离	
		}
		Vx = 0.0;  Vy = 0.0;
		Vxx = 0.0; Vxy = 0.0;
		Vyx = 0.0; Vyy = 0.0;
		for (i = 0; i < 4; i++) {
			k = (i < 3) ? (i + 1) : (i - 3);	     //对每条边处理   
			t1 = xe[k];
			t1 -= xe[i];
			t1 *= t1;
			t2 = ye[k];
			t2 -= ye[i];
			t2 *= t2;
			t2 += t1;
			d = sqrt(t2);                //计算边长
			//value D
			t1 = r[k];
			t1 += r[i];
			t1 *= t1;
			t2 = d;
			t2 *= t2;
			t1 -= t2;
			t1 *= 0.5;
			D = t1;              // 计算几何参数：D = 0.5 * ((r[k] + r[i])? - d?);
			//value rx，归一化距离向量
			t1 = xp1;
			t1 -= xe[i];
			t1 /= r[i];
			t2 = xp1;
			t2 -= xe[k];
			t2 /= r[k];
			t2 += t1;
			rx = t2;
			//rx=(xp1-xe[i])/r[i]+(xp1-xe[k])/r[k];	   
			//value ry
			t1 = yp1;
			t1 -= ye[i];
			t1 /= r[i];
			t2 = yp1;
			t2 -= ye[k];
			t2 /= r[k];
			t2 += t1;
			ry = t2;
			//ry=(yp1-ye[i])/r[i]+(yp1-ye[k])/r[k];	   

			//component of velocity potential			  				  	
			vr = r[k];
			vr += r[i];
			vr -= d;
			vr1 = r[k];
			vr1 += r[i];
			vr1 += d;
			if (d == 0.0) { Vx += 0.0; Vy += 0.0; }
			else {
				vt = ye[k]; vt -= ye[i]; vt *= log(vr / vr1); vt /= d;
				Vx += vt;
				vt = xe[i]; vt -= xe[k]; vt *= log(vr / vr1); vt /= d;
				Vy += vt;
			}
			//derivative of velocity potential components
			vt = ye[k];
			vt -= ye[i];
			vt *= rx;
			vt /= D;
			Vxx += vt;
			vt = ye[k];
			vt -= ye[i];
			vt *= ry;
			vt /= D;
			Vxy += vt;
			vt = xe[i];
			vt -= xe[k];
			vt *= rx;
			vt /= D;
			Vyx += vt;
			vt = xe[i];
			vt -= xe[k];
			vt *= ry;
			vt /= D;
			Vyy += vt;
		}
		detJ = Vxx * Vyy;
		detJ -= Vyx * Vxy;
		Z1 = Vyy * Vx;
		Z1 -= Vxy * Vy;
		Z1 /= detJ;
		Z2 = Vxx * Vy;
		Z2 -= Vyx * Vx;
		Z2 /= detJ;

		//Element coordinate system null points in the z=0 plane
		xn = xp1;
		xn -= Z1;
		yn = yp1;
		yn -= Z2;

		//收敛性检查，梯度变化小于阈值
		//if(j>=2 && fabs(Vx-Vx1)<0.0001 && fabs(Vy-Vy1)<0.0001) break;
		if (j >= 2 && fabs(Vx - Vx1) < 0.0001 && fabs(Vy - Vy1) < 0.0001)
		{
			//transformation from ecs to rcs for null points					
			//Vector cp(2);
			//cp[0]=xn;
			//cp[1]=yn;
			//null=tmt*cp;
			//null+=cnt1;

			break;
		}
		xp1 = xn;
		yp1 = yn;
		Vx1 = Vx;
		Vy1 = Vy;
	}
}
void Element::NullCentroid(int& k, double& xn, double& yn, Vector2d& cp, Vector3d& null,
	Vector3d& cnt, Eigen::Matrix<double, 3, 2, Eigen::ColMajor>& tmeT)
{
	//-----------------
	// null points
	//-----------------
	cp(0) = xn; cp(1) = yn;
	null = tmeT * cp;
	null += cnt;           //null是xn在三维中的坐标

	// xnr[k]=null[0];
	// ynr[k]=null[1];
	// znr[k]=null[2];
	xnr[k] = cnt[0];
	ynr[k] = cnt[1];
	znr[k] = cnt[2];
	//-------------------
	// centroid points 
	//-------------------
	xcr[k] = cnt[0];
	ycr[k] = cnt[1];
	zcr[k] = cnt[2];
}

void Element::TransforMatrix(int& k, Vector3d& tv, Vector3d& sv, Vector3d& nv)
{
	A11[k] = tv[0];  A12[k] = tv[1];  A13[k] = tv[2];
	A21[k] = sv[0];  A22[k] = sv[1];  A23[k] = sv[2];
	A31[k] = nv[0];  A32[k] = nv[1];  A33[k] = nv[2];
}
void Element::MomentInertia(int& k, Vector4d& xe, Vector4d& ye)
{
	double Ixx1, Iyy1, Ixy1, s1, s2, s3, st, t1, t2;

	s2 = 0.5 * (xe[2] - xe[0]) * (ye[1] - ye[3]);
	if (k == 0) SA = 2.0 * s2;
	else	    SA += 2.0 * s2;
	Area[k] = s2;

	//maximum diagonals
	s1 = ye[2] - ye[0];
	t1 = xe[2] - xe[0];
	t1 = sqrt(s1 * s1 + t1 * t1);
	//t1=sqrt((xe[2]-xe[0])*(xe[2]-xe[0])+(ye[2]-ye[0])*(ye[2]-ye[0]));		
	s1 = ye[3] - ye[1];
	t2 = xe[3] - xe[1];
	t2 = sqrt(s1 * s1 + t2 * t2);
	//t2=sqrt((xe[3]-xe[1])*(xe[3]-xe[1])+(ye[3]-ye[1])*(ye[3]-ye[1]));		
	td[k] = (t1 > t2) ? t1 : t2;

	//second moment of element
	Ixy1 = 2.0 * ye[0] + ye[1] + ye[3];
	Ixy1 *= ((xe[0] + xe[2]) * (ye[1] - ye[3]));
	Ixy1 -= 2.0 * xe[1] * (ye[0] * ye[0] - ye[1] * ye[1]);
	Ixy1 += 2.0 * xe[3] * (ye[0] * ye[0] - ye[3] * ye[3]);
	Ixy1 *= (xe[2] - xe[0]);
	Ixy1 /= 24.0;
	ixy[k] = Ixy1;

	st = xe[0];
	st += xe[1];
	st += xe[2];
	s2 = st;
	st += xe[3];
	s1 = st;
	st -= xe[1];
	s3 = st;
	t1 = xe[3];
	t1 -= xe[1];
	s1 *= t1;
	//s1*=(xe[3]-xe[1]);
	s1 *= ye[0];
	s2 *= ye[1];
	s2 *= xe[1];
	s3 *= ye[3];
	s3 *= xe[3];
	Ixx1 = xe[0] * xe[0];
	t1 = xe[0];
	t1 *= xe[2];
	Ixx1 += t1;
	//Ixx1+=xe[0]*xe[2];
	Ixx1 += xe[2] * xe[2];
	t1 = ye[1];
	t1 -= ye[3];
	Ixx1 *= t1;
	//Ixx1*=(ye[1]-ye[3]);
	Ixx1 += s1;
	Ixx1 += s2;
	Ixx1 -= s3;
	t1 = xe[2];
	t1 -= xe[0];
	Ixx1 *= t1;
	//Ixx1*=(xe[2]-xe[0]);
	Ixx1 /= 12.0;
	ixx[k] = Ixx1;


	Iyy1 = ye[0];
	Iyy1 += ye[1];
	Iyy1 += ye[3];
	Iyy1 *= Iyy1;
	t1 = ye[0];
	t1 *= ye[1];
	Iyy1 -= t1;
	//Iyy1-=ye[0]*ye[1];
	t1 = ye[0];
	t1 *= ye[3];
	Iyy1 -= t1;
	//Iyy1-=ye[0]*ye[3];
	t1 = ye[1];
	t1 *= ye[3];
	Iyy1 -= t1;
	//Iyy1-=ye[1]*ye[3];
	t1 = xe[2];
	t1 -= xe[0];
	Iyy1 *= t1;
	//Iyy1*=(xe[2]-xe[0]);
	t1 = ye[1];
	t1 -= ye[3];
	Iyy1 *= t1;
	//Iyy1*=(ye[1]-ye[3]);
	Iyy1 /= 12.0;
	iyy[k] = Iyy1;
}

void Element::LineValues(int& k, int& f, ElementMatrix& point)
{
	int i = 0;
	double vt1, vt2, vt3, vt4, ev1, ev2, ev3, ev4;
	double EPS = 1e-4;

	const Vector3d& v1 = point.row(0);
	const Vector3d& v2 = point.row(1);
	const Vector3d& v3 = point.row(2);
	const Vector3d& v4 = point.row(3);

	ev1 = v1[1]; ev2 = v2[1]; ev3 = v3[1]; ev4 = v4[1];

	//检查各顶点z坐标绝对值
	vt1 = fabs(v1[2]);
	vt2 = fabs(v2[2]);
	vt3 = fabs(v3[2]);
	vt4 = fabs(v4[2]);

	//统计z坐标接近0的顶点数量
	if (vt1 <= EPS) ++i;
	if (vt2 <= EPS) ++i;
	if (vt3 <= EPS) ++i;
	if (vt4 <= EPS) ++i;

	//x坐标大的点放在前
	//if (i == 2 && (ev1 > 0.0 || ev2 > 0.0 || ev3 > 0.0 || ev4 > 0.0)) {
	/*if (i == 2) {
		if (vt1 <= EPS && vt2 <= EPS) {
			if (v1[0] > v2[0]) {
				xpl(f, 0) = v1[0];  xpl(f, 1) = v2[0];
				ypl(f, 0) = v1[1];  ypl(f, 1) = v2[1];
				PotL[f] = k;
			}
			else if (v1[0] < v2[0]) {
				xpl(f, 0) = v2[0];  xpl(f, 1) = v1[0];
				ypl(f, 0) = v2[1];  ypl(f, 1) = v1[1];
				PotL[f] = k;
			}
			else  f--;
		}
		if (vt1 <= EPS && vt3 <= EPS) {
			if (v1[0] > v3[0]) {
				xpl(f, 0) = v1[0];  xpl(f, 1) = v3[0];
				ypl(f, 0) = v1[1];  ypl(f, 1) = v3[1];
				PotL[f] = k;
			}
			else if (v1[0] < v3[0]) {
				xpl(f, 0) = v3[0];  xpl(f, 1) = v1[0];
				ypl(f, 0) = v3[1];  ypl(f, 1) = v1[1];
				PotL[f] = k;
			}
			else  f--;
		}
		if (vt1 <= EPS && vt4 <= EPS) {
			if (v1[0] > v4[0]) {
				xpl(f, 0) = v1[0];  xpl(f, 1) = v4[0];
				ypl(f, 0) = v1[1];  ypl(f, 1) = v4[1];
				PotL[f] = k;
			}
			else if (v1[0] < v4[0]) {
				xpl(f, 0) = v4[0];  xpl(f, 1) = v1[0];
				ypl(f, 0) = v4[1];  ypl(f, 1) = v1[1];
				PotL[f] = k;
			}
			else  f--;
		}
		if (vt2 <= EPS && vt3 <= EPS) {
			if (v2[0] > v3[0]) {
				xpl(f, 0) = v2[0];  xpl(f, 1) = v3[0];
				ypl(f, 0) = v2[1];  ypl(f, 1) = v3[1];
				PotL[f] = k;
			}
			else if (v2[0] < v3[0]) {
				xpl(f, 0) = v3[0];  xpl(f, 1) = v2[0];
				ypl(f, 0) = v3[1];  ypl(f, 1) = v2[1];
				PotL[f] = k;
			}
			else  f--;
		}
		if (vt2 <= EPS && vt4 <= EPS) {
			if (v2[0] > v4[0]) {
				xpl(f, 0) = v2[0];  xpl(f, 1) = v4[0];
				ypl(f, 0) = v2[1];  ypl(f, 1) = v4[1];
				PotL[f] = k;
			}
			else if (v2[0] < v4[0]) {
				xpl(f, 0) = v4[0];  xpl(f, 1) = v2[0];
				ypl(f, 0) = v4[1];  ypl(f, 1) = v2[1];
				PotL[f] = k;
			}
			else  f--;
		}
		if (vt3 <= EPS && vt4 <= EPS) {
			if (v3[0] > v4[0]) {
				xpl(f, 0) = v3[0];  xpl(f, 1) = v4[0];
				ypl(f, 0) = v3[1];  ypl(f, 1) = v4[1];
				PotL[f] = k;
			}
			else if (v3[0] < v4[0]) {
				xpl(f, 0) = v4[0];  xpl(f, 1) = v3[0];
				ypl(f, 0) = v4[1];  ypl(f, 1) = v3[1];
				PotL[f] = k;
			}
			else  f--;
		}
	}
	else
		f--;*/


	if (i == 2)
	{
		auto setLineByCross = [&](const Vector3d& a, const Vector3d& b)
		{
			const double ax = a[0], ay = a[1];
			const double bx = b[0], by = b[1];

			// 2D 叉乘 z 分量：a x b
			const double cross = ax * by - ay * bx;
			const double EPS_CROSS = 1e-12;

			// a x b > 0  => b 放前面
			// 否则       => a 放前面
			if (cross > EPS_CROSS)
			{
				xpl(f, 0) = bx;  xpl(f, 1) = ax;
				ypl(f, 0) = by;  ypl(f, 1) = ay;
				PotL[f] = k;
			}
			else if (cross < -EPS_CROSS)
			{
				xpl(f, 0) = ax;  xpl(f, 1) = bx;
				ypl(f, 0) = ay;  ypl(f, 1) = by;
				PotL[f] = k;
			}
			else
			{
				// 两点与原点近乎共线，退化时再用 x/y 做兜底
				if (ax > bx || (std::fabs(ax - bx) <= EPS && ay > by))
				{
					xpl(f, 0) = ax;  xpl(f, 1) = bx;
					ypl(f, 0) = ay;  ypl(f, 1) = by;
					PotL[f] = k;
				}
				else if (ax < bx || (std::fabs(ax - bx) <= EPS && ay < by))
				{
					xpl(f, 0) = bx;  xpl(f, 1) = ax;
					ypl(f, 0) = by;  ypl(f, 1) = ay;
					PotL[f] = k;
				}
				else
				{
					f--;
				}
			}
		};

		if (vt1 <= EPS && vt2 <= EPS) setLineByCross(v1, v2);
		if (vt1 <= EPS && vt3 <= EPS) setLineByCross(v1, v3);
		if (vt1 <= EPS && vt4 <= EPS) setLineByCross(v1, v4);
		if (vt2 <= EPS && vt3 <= EPS) setLineByCross(v2, v3);
		if (vt2 <= EPS && vt4 <= EPS) setLineByCross(v2, v4);
		if (vt3 <= EPS && vt4 <= EPS) setLineByCross(v3, v4);
	}
	else
	{
		f--;
	}
}

void Element::NullCentroidDis(int j, int i, double& dx, double& dy, double& dz)
{
	dx = xnr[j] - xcr[i];
	dy = ynr[j] - ycr[i];
	dz = znr[j] - zcr[i];
}

void Element::MonopoleCal(int j, int i, double r0)
{

	double r01 = r0 * r0 * r0;

	Vxr = Area[i] * (xnr[j] - xcr[i]) / r01;
	Vyr = Area[i] * (ynr[j] - ycr[i]) / r01;
	Vzr = Area[i] * (znr[j] - zcr[i]) / r01;
	Ph = Area[i] / r0;
}

void Element::DepoleCal(int j, int i, double r0, double& Vxe, double& Vye, double& Vze)
{
	double  R, p, q, Wx, Wy, Wz, Wxxx, Wyyy, Wxxy, Wxyy,
		Wxxz, Wxyz, Wyyz, xne, yne, zne;
	double r_2, r_3, r_5, w, wxx, wyy, wxy, tx, ty, tz;

	//计算场点和源点在局部坐标中的距离向量（xne,yne,zne)
	ElementCoor(j, i, xne, yne, zne);

	tx = xne * xne;
	ty = yne * yne;
	tz = zne * zne;

	double r02 = r0 * r0;
	double r03 = r02 * r0;

	R = r03 * r03 * r0;
	r_2 = 1.0 / r02;
	r_3 = 1.0 / r03;
	r_5 = 1.0 / (r03 * r02);
	w = 1.0 / r0;

	p = (ty + tz - 4 * tx);
	q = (tx + tz - 4 * ty);

	Wx = -xne / r03;
	Wy = -yne / r03;
	Wz = -zne / r03;

	wxx = r_3 * (3.0 * tx * r_2 - 1.0);
	wyy = r_3 * (3.0 * ty * r_2 - 1.0);
	wxy = 3.0 * xne * yne * r_5;

	Wxxx = 3 * xne * (3 * p + 10 * xne * xne) / R;
	Wyyy = 3 * yne * (3 * q + 10 * yne * yne) / R;
	Wxxy = 3 * yne * p / R;
	Wxyy = 3 * xne * q / R;
	Wxxz = 3 * zne * p / R;
	Wxyz = -15 * xne * yne * zne / R;
	Wyyz = 3 * zne * q / R;

	Vxe = -(Area[i] * Wx + 0.5 * ixx[i] * Wxxx + ixy[i] * Wxxy + 0.5 * iyy[i] * Wxyy);
	Vye = -(Area[i] * Wy + 0.5 * ixx[i] * Wxxy + ixy[i] * Wxyy + 0.5 * iyy[i] * Wyyy);
	Vze = -(Area[i] * Wz + 0.5 * ixx[i] * Wxxz + ixy[i] * Wxyz + 0.5 * iyy[i] * Wyyz);
	Ph = (Area[i] * w + 0.5 * ixx[i] * wxx + ixy[i] * wxy + 0.5 * iyy[i] * wyy);
}

void Element::ElementCoor(int j, int i, double& xne, double& yne, double& zne)
{
	//将场点和源点间的距离向量从整体坐标转换到四边形局部坐标（tx,ty,tz) to (xne,yne,zne)

	xne = A11[i] * (xnr[j] - xcr[i]) + A12[i] * (ynr[j] - ycr[i]) +
		A13[i] * (znr[j] - zcr[i]);
	yne = A21[i] * (xnr[j] - xcr[i]) + A22[i] * (ynr[j] - ycr[i]) +
		A23[i] * (znr[j] - zcr[i]);
	zne = A31[i] * (xnr[j] - xcr[i]) + A32[i] * (ynr[j] - ycr[i]) +
		A33[i] * (znr[j] - zcr[i]);
}

void Element::ExactCal(int j, int i, double& Vxe, double& Vye, double& Vze)
{
	int k, m;
	double xne, yne, zne;
	double tmp1, ce[4], se[4], R[4], Q, Ji, dt1, sgn;
	std::array<std::array<double, 2>, 4> se1;
	double re[4], de[4], tmp, Ph0;

	//计算面元j null point与面元i质心间的空间距离向量并转到面元坐标系
	ElementCoor(j, i, xne, yne, zne);

	//de:四边形四条边的长度，ce:沿每条边方向的单位向量的x分量
	//se:沿每条边方向的单位向量的y分量
	//re:null point到顶点的距离，R：null point到边的距离有向距离
	ExactDist(i, xne, yne, zne, de, se, ce, re, R);

	//se1:顶点到null point的向量在边方向上的投影长度
	ExactDis(i, xne, yne, ce, se, se1);

	double c_Sz[5], h_Sz[5], m_Sz[5], r_Sz[5];
	compute_Sz(j, i, xne, yne, zne, c_Sz, h_Sz, m_Sz, r_Sz);

	Vxe = Vye = Vze = Ph0 = 0.0;

	double Vze1 = 0;
	for (k = 0; k < 4; k++) {
		m = (k < 3) ? (k + 1) : (k - 3);

		Q = re[k];
		Q += re[m];
		tmp = Q;
		tmp += de[k];
		tmp1 = Q;
		tmp1 -= de[k];
		Q = tmp;
		Q /= tmp1;
		Q = log(Q);
		tmp = -se[k];
		tmp *= Q;
		Vxe += tmp;
		tmp = ce[k];
		tmp *= Q;
		Vye += tmp;
		//Vxe = Σ(-se[k] * ln((r_k + r_{k+1} + d_k)/(r_k + r_{k+1} - d_k)))
		//Vye = Σ(ce[k] * ln((r_k + r_{k+1} + d_k)/(r_k + r_{k+1} - d_k)))

		tmp = re[k];
		tmp *= se1[k][1];
		tmp1 = re[m];
		tmp1 *= se1[k][0];
		tmp -= tmp1;
		tmp *= R[k];
		tmp *= fabs(zne);
		tmp1 = re[k];
		tmp1 *= re[m];
		tmp1 *= pow(R[k], 2.0);
		Ji = pow(zne, 2.0);
		Ji *= se1[k][0];
		Ji *= se1[k][1];
		Ji += tmp1;
		Ji = atan2(tmp, Ji);
		Vze += -Ji;
		//R[k] = (xne-xpe[k])*se[k] - (yne-ype[k])*ce[k]
		//Vze = -Σ(atan2(|z|*R_k*(r_k*s_{k,1} - r_{k+1}*s_{k,0}), 
			  // z?*s_{k,0}*s_{k,1} + r_k*r_{k+1}*R_k?))

		//Vze += (atan((m_Sz[k]*c_Sz[k]-h_Sz[k])/(zne*r_Sz[k]))-atan((m_Sz[k]*c_Sz[m]-h_Sz[m])/(zne*r_Sz[m])));


		tmp = R[k];
		// tmp=-R[k];        //by JS
		tmp *= Q;
		tmp1 = fabs(zne);
		tmp1 *= Ji;
		// tmp1=zne*Vze;
		tmp += tmp1;
		Ph0 += tmp;
		//Ph0 = Σ(R_k * ln((r_k + r_{k+1} + d_k)/(r_k + r_{k+1} - d_k)) + |z| * Ji)
	}
	// Ph0+=zne*Vze;
	// std::cout<<Vze<<","<<Vze1<<std::endl;

	//当场点在源面元内时dt1=2*pi，其余为0
	ExactSolu(de, R, dt1);

	//场点在面元z轴正向或负向
	// sgn=zne>0?1.0:-1.0;
	sgn = zne > 0 ? 1.0 : -1.0;

	//i=j表示同一个面元
	// Vze=(i==j)?(dt1+Vze):sgn*(dt1+Vze);
	Vze = (i == j) ? (dt1 + Vze) : sgn * (dt1 + Vze);
	Ph = Ph0 - fabs(zne) * dt1;



	//  Ph=Ph0;		
}

void Element::ReferenceCoor(int i, double& Vxe, double& Vye, double& Vze)
{
	Vxr = A11[i] * Vxe + A21[i] * Vye + A31[i] * Vze;
	Vyr = A12[i] * Vxe + A22[i] * Vye + A32[i] * Vze;
	Vzr = A13[i] * Vxe + A23[i] * Vye + A33[i] * Vze;
}

void Element::ExactDist(int& i, double& xne, double& yne, double& zne, double* de, double* se, double* ce, double* re, double* R)
{
	//点的索引
	int k, m;

	//一条边的x,y方向长度和总长ti
	double tx, ty, tr, ti;

	//四个顶点（四条边）
	for (k = 0; k < 4; k++) {
		m = (k < 3) ? (k + 1) : (k - 3);      //相邻顶点索引  

		tx = xpe(i, m); tx -= xpe(i, k); tr = tx * tx;
		ty = ype(i, m); ty -= ype(i, k); ti = ty * ty;
		ti += tr; ti = sqrt(ti);
		de[k] = ti;

		if (de[k] == 0) { ce[k] = 0.0; se[k] = 0.0; }
		//单位化，ce,se第k条边单位方向向量（cos和sin）
		else { tx /= de[k]; ty /= de[k]; ce[k] = tx; se[k] = ty; }

		//场点到顶点的距离re
		tx = xne; tx -= xpe(i, k); tr = tx * tx;
		ty = yne; ty -= ype(i, k); ti = ty * ty;
		ti += tr; ti += (zne * zne); ti = sqrt(ti);
		re[k] = ti;

		//计算场点到边的有向距离：R[k] = (xne-xpe[k])*se[k] - (yne-ype[k])*ce[k]
		//如果场点在边左侧则R大于0，在右侧则小于0（相当于叉积）
		tr = tx; tr *= se[k]; ti = ty; ti *= ce[k]; tr -= ti;
		R[k] = tr;
	}
}

void Element::ExactDis(int& i, double& xne, double& yne, double* ce, double* se, std::array<std::array<double, 2>, 4>& se1)
{
	int k, f, m;
	double ti, tr;
	//四条边
	for (k = 0; k < 4; k++) {
		//两个顶点
		for (f = 0; f < 2; f++) {
			if (k == 0) m = f;
			else if (k == 1) m = f + 1;
			else if (k == 2) m = f + 2;
			else if (k == 3 && f == 0) m = f + 3;
			else m = f - 1;

			//顶点到场点的y距离
			tr = ype(i, m); tr -= yne;
			//在边方向投影
			tr *= se[k];
			//顶点到场点的x距离
			ti = xpe(i, m); ti -= xne;
			//在边方向投影
			ti *= ce[k]; ti += tr;
			//存储
			se1[k][f] = ti;
		}
	}
}

void Element::compute_Sz(int j, int i, double xne, double yne, double zne,
	double* c_Sz, double* h_Sz, double* m_Sz, double* r_Sz)
{
	int m;
	for (int k = 0; k < 4; k++)
	{
		m = (k < 3) ? (k + 1) : (k - 3);
		c_Sz[k] = pow(xpe(i, k) - xne, 2) + zne * zne;
		h_Sz[k] = (xpe(i, k) - xne) * (ype(i, k) - yne);
		m_Sz[k] = (ype(i, m) - ype(i, k)) / (xpe(i, m) - xpe(i, k));
		r_Sz[k] = sqrt(c_Sz[k] + pow(ype(i, k) - yne, 2));
	}
	c_Sz[4] = c_Sz[0];
	h_Sz[4] = h_Sz[0];
	m_Sz[4] = m_Sz[0];
	r_Sz[4] = r_Sz[0];
}

void Element::ExactSolu(double* de, double* R, double& dt1)
{
	int i, d, r, r1, n = 4;
	d = r = r1 = 0;
	for (i = 0; i < n; ++i) {
		if (de[i] != 0.0) ++d;   //d:有效边的数量
		if (R[i] > 0.0) ++r;	  //正有向距离边的数量
		if (R[i] >= 0.0) ++r1;
	}
	//四条边的长度均不为0
	if (d == n) {
		//null point在面元内部
		if (r == n) dt1 = 2.0 * PI;
		else dt1 = 0.0;
	}
	else {
		if (r1 == n) dt1 = 2.0 * PI;
		else dt1 = 0.0;
	}
}

void Element::eNormalVelocity(int i, double& Rn)
{
	double ti;

	ti = A33[i];
	ti *= Vzr;
	Rn = A32[i];
	Rn *= Vyr;
	Rn += ti;
	ti = A31[i];
	ti *= Vxr;
	Rn += ti;
}

void Element::SignChangez(int i)
{
	zcr[i] = -zcr[i];
	A13[i] = -A13[i];
	A23[i] = -A23[i];
	A33[i] = -A33[i];
}

void Element::PointtoMesh(const ElementMatrix& point, const int idex)
{
	Mesh singleMesh{};

	for(int i=0;i<3;++i)
	{
		singleMesh.point1[i] = point(0, i);
		singleMesh.point2[i] = point(1, i);
		singleMesh.point3[i] = point(2, i);
		singleMesh.point4[i] = point(3, i);
	}

	meshes.Meshs[idex] = singleMesh;
	meshes.areas[idex] = singleMesh.area();
	meshes.centers[idex] = singleMesh.center_p();
	meshes.Nor_vectors[idex] = singleMesh.Nor_vector();
}

void Element::saveRankine(const std::string filePath)
{
	std::ofstream outFile(filePath);
	if (!outFile.is_open()) {
		std::cerr << "无法打开文件\t"<<filePath<<"\t进行写入！" << std::endl;
		return;
	}
	for (int i = 0; i < NE; ++i)
	{

		for(int j=0;j<NE;++j)
		{
			outFile << A_rankine(i, j) << ",";
		}
		outFile << std::endl;
	}
	outFile.close();
}


void Element::Rankine(int j, int i)
{
	double Vxe, Vye, Vze, r0;

	//xnr[0] = 172.56;
	//ynr[0] = 0.62794;
	//znr[0] = -13.714;

	//源点与场点间的距离
	NullCentroidDis(j, i, Vxe, Vye, Vze);

	r0 = sqrt(Vxe * Vxe + Vye * Vye + Vze * Vze);

	//根据距离选择不同的算法，td为面元最大对角线长度
	if (r0 >= 4 * td[i])
		//远场
		MonopoleCal(j, i, r0);
	else if (r0 >= 2.45 * td[i])
		//中场  	
		DepoleCal(j, i, r0, Vxe, Vye, Vze);
	else if (r0 < 2.45 * td[i])
		//近场
		ExactCal(j, i, Vxe, Vye, Vze);
	if (r0 < 4.0 * td[i])
		//针对近场和中场进行坐标转换到全局，远场直接在全局坐标计算，无需坐标转换
		ReferenceCoor(i, Vxe, Vye, Vze);
};