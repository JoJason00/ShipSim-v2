#include <array>
#include "source_dipole.h"
#include <cmath>
#include <iostream>

namespace
{
	std::array<double, 3> cross_product(std::array<double, 3>& a, std::array<double, 3>& b)
	{
		std::array<double, 3> result;
		result[0] = a[1] * b[2] - a[2] * b[1];
		result[1] = a[2] * b[0] - a[0] * b[2];
		result[2] = a[0] * b[1] - a[1] * b[0];
		return result;
	}

	inline double dot_product(std::array<double, 3>& a, std::array<double, 3>& b)
	{
		double result = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
		return result;
	}
}


std::array<double, 3> Mesh::center_p()
{
	std::array<double, 3> p0;
	p0[0] = (point1[0] + point2[0] + point3[0] + point4[0]) / 4.0;
	p0[1] = (point1[1] + point2[1] + point3[1] + point4[1]) / 4.0;
	p0[2] = (point1[2] + point2[2] + point3[2] + point4[2]) / 4.0;
	return p0;
}

std::array<double, 3> Mesh::Nor_vector()
{
	std::array<double, 3> p13 = { point3[0] - point1[0], point3[1] - point1[1], point3[2] - point1[2] };
	std::array<double, 3> p24 = { point4[0] - point2[0], point4[1] - point2[1], point4[2] - point2[2] };
	std::array<double, 3> nor = cross_product(p13, p24);
	double mod = sqrt(pow(nor[0], 2) + pow(nor[1], 2) + pow(nor[2], 2));
	nor[0] /= mod;
	nor[1] /= mod;
	nor[2] /= mod;
	return nor;
}

//double Mesh::area()
//{
//	std::array<double, 3> p13 = { point3[0] - point1[0], point3[1] - point1[1], point3[2] - point1[2] };
//	std::array<double, 3> p24 = { point4[0] - point2[0], point4[1] - point2[1], point4[2] - point2[2] };
//	std::array<double, 3> nor = cross_product(p13, p24);
//	double mod = sqrt(pow(nor[0], 2) + pow(nor[1], 2) + pow(nor[2], 2)) / 2;
//	return mod;
//}

double Mesh::area()
{
	auto sub = [](const std::array<double, 3>& a, const std::array<double, 3>& b) {
		return std::array<double, 3>{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
	};
	auto norm3 = [](const std::array<double, 3>& v) {
		return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
	};
	auto p12 = sub(point2, point1), p13 = sub(point3, point1);
	auto p14 = sub(point4, point1);
	auto t1 = cross_product(p12, p13);
	auto t2 = cross_product(p13, p14);
	return 0.5 * (norm3(t1) + norm3(t2));
}


//赋值、更新，整体坐标转换到局部坐标，局部坐标xy平面在四边形面元平面上，z轴为其法线方向，指向物体内
void Source_Dipole::update(std::array<double, 3>& qq, Mesh& mesh)
{
	q = qq;
	p = &mesh;
	p0 = mesh.center_p();
	source = 0;
	dipole = 0;
	sv_local = { 0,0,0 };
	sv = { 0,0,0 };
	sv3 = nullptr;
	dv = { 0,0,0 };

	//局部坐标轴基矢量
	//x轴取为点1指向3
	local[0] = { (*p).point3[0] - (*p).point1[0], (*p).point3[1] - (*p).point1[1], (*p).point3[2] - (*p).point1[2] };
	double mod = sqrt(pow(local[0][0], 2) + pow(local[0][1], 2) + pow(local[0][2], 2));
	local[0][0] /= mod;
	local[0][1] /= mod;
	local[0][2] /= mod;

	std::array<double, 3> p24;
	p24 = { (*p).point4[0] - (*p).point2[0], (*p).point4[1] - (*p).point2[1], (*p).point4[2] - (*p).point2[2] };
	mod = sqrt(pow(p24[0], 2) + pow(p24[1], 2) + pow(p24[2], 2));
	p24[0] /= mod;
	p24[1] /= mod;
	p24[2] /= mod;

	//z轴取为13向量和24向量的叉积
	local[2] = cross_product(local[0], p24);
	mod = sqrt(pow(local[2][0], 2) + pow(local[2][1], 2) + pow(local[2][2], 2));
	local[2][0] /= mod;
	local[2][1] /= mod;
	local[2][2] /= mod;

	//则y轴为z轴和x轴的叉积
	local[1] = cross_product(local[2], local[0]);

	//场点和四边形顶点的局部坐标
	std::array<double, 3> p0_q = { q[0] - p0[0],q[1] - p0[1],q[2] - p0[2] };
	q_local[0] = dot_product(p0_q, local[0]);
	q_local[1] = dot_product(p0_q, local[1]);
	q_local[2] = dot_product(p0_q, local[2]);

	std::array<double, 3> p0_p1 = { (*p).point1[0] - p0[0],(*p).point1[1] - p0[1],(*p).point1[2] - p0[2] };
	std::array<double, 3> p0_p2 = { (*p).point2[0] - p0[0],(*p).point2[1] - p0[1],(*p).point2[2] - p0[2] };
	std::array<double, 3> p0_p3 = { (*p).point3[0] - p0[0],(*p).point3[1] - p0[1],(*p).point3[2] - p0[2] };
	std::array<double, 3> p0_p4 = { (*p).point4[0] - p0[0],(*p).point4[1] - p0[1],(*p).point4[2] - p0[2] };
	p_local[0][0] = dot_product(p0_p1, local[0]);
	p_local[0][1] = dot_product(p0_p1, local[1]);
	p_local[1][0] = dot_product(p0_p2, local[0]);
	p_local[1][1] = dot_product(p0_p2, local[1]);
	p_local[2][0] = dot_product(p0_p3, local[0]);
	p_local[2][1] = dot_product(p0_p3, local[1]);
	p_local[3][0] = dot_product(p0_p4, local[0]);
	p_local[3][1] = dot_product(p0_p4, local[1]);

	int j;
	for (int i = 0; i < 4; i++)
	{
		//r_array为场点到四边形顶点的距离
		r_array[i] = { q_local[0] - p_local[i][0],q_local[1] - p_local[i][1] ,q_local[2] };
		r[i] = sqrt(pow(r_array[i][0], 2) + pow(r_array[i][1], 2) + pow(r_array[i][2], 2));

		//l_array为四边形的边长
		j = (i + 1) % 4;
		l_array[i] = { p_local[j][0] - p_local[i][0],p_local[j][1] - p_local[i][1] };
		l[i] = sqrt(pow(l_array[i][0], 2) + pow(l_array[i][1], 2));
	}
}

//返回面元的面积
double Source_Dipole::Gauss_integrate()
{
	//积分结果
	double result = 0;
	//雅可比矩阵行列式
	double Jacobi = 0;

	//2*2高斯点映射坐标
	std::array<double, 2> ksi = { 1.0 / sqrt(3),-1.0 / sqrt(3) };
	std::array<double, 2> eta = { 1.0 / sqrt(3),-1.0 / sqrt(3) };

	//高斯点局部坐标
	/*std::array<double, 3> gauss_local{};
	gauss_local[2] = 0.0;*/

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			Jacobi = ((-(1 - eta[j]) * p_local[0][0] + (1 - eta[j]) * p_local[1][0] + (1 + eta[j]) * p_local[2][0] - (1 + eta[j]) * p_local[3][0])
				* (-(1 - ksi[i]) * p_local[0][1] - (1 + ksi[i]) * p_local[1][1] + (1 + ksi[i]) * p_local[2][1] + (1 - ksi[i]) * p_local[3][1])
				- (-(1 - ksi[i]) * p_local[0][0] - (1 + ksi[i]) * p_local[1][0] + (1 + ksi[i]) * p_local[2][0] + (1 - ksi[i]) * p_local[3][0])
				* (-(1 - eta[j]) * p_local[0][1] + (1 - eta[j]) * p_local[1][1] + (1 + eta[j]) * p_local[2][1] - (1 + eta[j]) * p_local[3][1])) / 16.0;

			/*gauss_local[0] = ((1 - ksi[i]) * (1 - eta[j]) * p_local[0][0] + (1 + ksi[i]) * (1 - eta[j]) * p_local[1][0]
				+ (1 + ksi[i]) * (1 + eta[j]) * p_local[2][0] + (1 - ksi[i]) * (1 + eta[j]) * p_local[3][0]) / 4.0;
			gauss_local[1] = ((1 - ksi[i]) * (1 - eta[j]) * p_local[0][1] + (1 + ksi[i]) * (1 - eta[j]) * p_local[1][1]
				+ (1 + ksi[i]) * (1 + eta[j]) * p_local[2][1] + (1 - ksi[i]) * (1 + eta[j]) * p_local[3][1]) / 4.0;*/

			result += Jacobi;
		}
	}
	return result;
}

//2*2高斯积分（返回一个数）
double Source_Dipole::Gauss_integrate(double(*operation)(std::array<double, 3>&, std::array<double, 3>&))
{
	//积分结果
	double result = 0;
	//雅可比矩阵行列式
	double Jacobi = 0;

	//2*2高斯点映射坐标
	std::array<double, 2> ksi = { 1.0 / sqrt(3),-1.0 / sqrt(3) };
	std::array<double, 2> eta = { 1.0 / sqrt(3),-1.0 / sqrt(3) };

	//高斯点局部坐标及整体坐标
	std::array<double, 2> gauss_local{};
	std::array<double, 3> gauss_global{};
	/*gauss_local[2] = 0.0;*/

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			Jacobi = ((-(1 - eta[j]) * p_local[0][0] + (1 - eta[j]) * p_local[1][0] + (1 + eta[j]) * p_local[2][0] - (1 + eta[j]) * p_local[3][0])
				* (-(1 - ksi[i]) * p_local[0][1] - (1 + ksi[i]) * p_local[1][1] + (1 + ksi[i]) * p_local[2][1] + (1 - ksi[i]) * p_local[3][1])
				- (-(1 - ksi[i]) * p_local[0][0] - (1 + ksi[i]) * p_local[1][0] + (1 + ksi[i]) * p_local[2][0] + (1 - ksi[i]) * p_local[3][0])
				* (-(1 - eta[j]) * p_local[0][1] + (1 - eta[j]) * p_local[1][1] + (1 + eta[j]) * p_local[2][1] - (1 + eta[j]) * p_local[3][1])) / 16.0;

			gauss_local[0] = ((1 - ksi[i]) * (1 - eta[j]) * p_local[0][0] + (1 + ksi[i]) * (1 - eta[j]) * p_local[1][0]
				+ (1 + ksi[i]) * (1 + eta[j]) * p_local[2][0] + (1 - ksi[i]) * (1 + eta[j]) * p_local[3][0]) / 4.0;
			gauss_local[1] = ((1 - ksi[i]) * (1 - eta[j]) * p_local[0][1] + (1 + ksi[i]) * (1 - eta[j]) * p_local[1][1]
				+ (1 + ksi[i]) * (1 + eta[j]) * p_local[2][1] + (1 - ksi[i]) * (1 + eta[j]) * p_local[3][1]) / 4.0;

			gauss_global[0] = local[0][0] * gauss_local[0] + local[1][0] * gauss_local[1] + p0[0];
			gauss_global[1] = local[0][1] * gauss_local[0] + local[1][1] * gauss_local[1] + p0[1];
			gauss_global[2] = local[0][2] * gauss_local[0] + local[1][2] * gauss_local[1] + p0[2];

			result += operation(q, gauss_global) * Jacobi;
		}
	}
	return result;
}

//4*4高斯积分，返回一个数
double Source_Dipole::Gauss_integrate44(double(*operation)(std::array<double, 3>&, std::array<double, 3>&))
{
	// 积分结果
	double result = 0.0;

	// 16×16高斯积分点的位置和权重
	const std::array<double, 16> gauss_points = {
		-0.9894009349916499,
		-0.9445750230732326,
		-0.8656312023878318,
		-0.7554044083550030,
		-0.6178762444026438,
		-0.4580167776572274,
		-0.2816035507792589,
		-0.0950125098376374,
		 0.0950125098376374,
		 0.2816035507792589,
		 0.4580167776572274,
		 0.6178762444026438,
		 0.7554044083550030,
		 0.8656312023878318,
		 0.9445750230732326,
		 0.9894009349916499
	};

	const std::array<double, 16> weights = {
		0.0271524594117541,
		0.0622535239386479,
		0.0951585116824928,
		0.1246289712555339,
		0.1495959888165767,
		0.1691565193950025,
		0.1826034150449236,
		0.1894506104550685,
		0.1894506104550685,
		0.1826034150449236,
		0.1691565193950025,
		0.1495959888165767,
		0.1246289712555339,
		0.0951585116824928,
		0.0622535239386479,
		0.0271524594117541
	};

	// 遍历所有16×16=256个高斯点组合
	for (int i = 0; i < 16; ++i) {
		for (int j = 0; j < 16; ++j) {
			const double ksi = gauss_points[i];
			const double eta = gauss_points[j];
			const double weight = weights[i] * weights[j];

			// 计算雅可比行列式
			double Jacobi = ((-(1 - eta) * p_local[0][0] + (1 - eta) * p_local[1][0] + (1 + eta) * p_local[2][0] - (1 + eta) * p_local[3][0])
				* (-(1 - ksi) * p_local[0][1] - (1 + ksi) * p_local[1][1] + (1 + ksi) * p_local[2][1] + (1 - ksi) * p_local[3][1])
				- (-(1 - ksi) * p_local[0][0] - (1 + ksi) * p_local[1][0] + (1 + ksi) * p_local[2][0] + (1 - ksi) * p_local[3][0])
				* (-(1 - eta) * p_local[0][1] + (1 - eta) * p_local[1][1] + (1 + eta) * p_local[2][1] - (1 + eta) * p_local[3][1])) / 16.0;

			// 计算高斯点在局部坐标系中的坐标
			std::array<double, 2> gauss_local;
			gauss_local[0] = ((1 - ksi) * (1 - eta) * p_local[0][0] + (1 + ksi) * (1 - eta) * p_local[1][0]
				+ (1 + ksi) * (1 + eta) * p_local[2][0] + (1 - ksi) * (1 + eta) * p_local[3][0]) / 4.0;
			gauss_local[1] = ((1 - ksi) * (1 - eta) * p_local[0][1] + (1 + ksi) * (1 - eta) * p_local[1][1]
				+ (1 + ksi) * (1 + eta) * p_local[2][1] + (1 - ksi) * (1 + eta) * p_local[3][1]) / 4.0;

			// 转换到全局坐标系
			std::array<double, 3> gauss_global;
			gauss_global[0] = local[0][0] * gauss_local[0] + local[1][0] * gauss_local[1] + p0[0];
			gauss_global[1] = local[0][1] * gauss_local[0] + local[1][1] * gauss_local[1] + p0[1];
			gauss_global[2] = local[0][2] * gauss_local[0] + local[1][2] * gauss_local[1] + p0[2];

			// 累加积分结果
			result += operation(q, gauss_global) * Jacobi * weight;
		}
	}

	return result;
}

//4*4高斯积分，返回一个向量
std::array<double, 3> Source_Dipole::Gauss_integrate44(std::array<double, 3>(*operation)(std::array<double, 3>&, std::array<double, 3>&))
{
	// 积分结果
	std::array<double, 3> result{};

	// 16×16高斯积分点的位置和权重
	const std::array<double, 16> gauss_points = {
		-0.9894009349916499,
		-0.9445750230732326,
		-0.8656312023878318,
		-0.7554044083550030,
		-0.6178762444026438,
		-0.4580167776572274,
		-0.2816035507792589,
		-0.0950125098376374,
		 0.0950125098376374,
		 0.2816035507792589,
		 0.4580167776572274,
		 0.6178762444026438,
		 0.7554044083550030,
		 0.8656312023878318,
		 0.9445750230732326,
		 0.9894009349916499
	};

	const std::array<double, 16> weights = {
		0.0271524594117541,
		0.0622535239386479,
		0.0951585116824928,
		0.1246289712555339,
		0.1495959888165767,
		0.1691565193950025,
		0.1826034150449236,
		0.1894506104550685,
		0.1894506104550685,
		0.1826034150449236,
		0.1691565193950025,
		0.1495959888165767,
		0.1246289712555339,
		0.0951585116824928,
		0.0622535239386479,
		0.0271524594117541
	};

	//存储operation的中间函数值
	std::array<double, 3> middle{};

	// 遍历所有16×16=256个高斯点组合
	for (int i = 0; i < 16; ++i) {
		for (int j = 0; j < 16; ++j) {
			const double ksi = gauss_points[i];
			const double eta = gauss_points[j];
			const double weight = weights[i] * weights[j];

			// 计算雅可比行列式
			double Jacobi = ((-(1 - eta) * p_local[0][0] + (1 - eta) * p_local[1][0] + (1 + eta) * p_local[2][0] - (1 + eta) * p_local[3][0])
				* (-(1 - ksi) * p_local[0][1] - (1 + ksi) * p_local[1][1] + (1 + ksi) * p_local[2][1] + (1 - ksi) * p_local[3][1])
				- (-(1 - ksi) * p_local[0][0] - (1 + ksi) * p_local[1][0] + (1 + ksi) * p_local[2][0] + (1 - ksi) * p_local[3][0])
				* (-(1 - eta) * p_local[0][1] + (1 - eta) * p_local[1][1] + (1 + eta) * p_local[2][1] - (1 + eta) * p_local[3][1])) / 16.0;

			// 计算高斯点在局部坐标系中的坐标
			std::array<double, 2> gauss_local;
			gauss_local[0] = ((1 - ksi) * (1 - eta) * p_local[0][0] + (1 + ksi) * (1 - eta) * p_local[1][0]
				+ (1 + ksi) * (1 + eta) * p_local[2][0] + (1 - ksi) * (1 + eta) * p_local[3][0]) / 4.0;
			gauss_local[1] = ((1 - ksi) * (1 - eta) * p_local[0][1] + (1 + ksi) * (1 - eta) * p_local[1][1]
				+ (1 + ksi) * (1 + eta) * p_local[2][1] + (1 - ksi) * (1 + eta) * p_local[3][1]) / 4.0;

			// 转换到全局坐标系
			std::array<double, 3> gauss_global;
			gauss_global[0] = local[0][0] * gauss_local[0] + local[1][0] * gauss_local[1] + p0[0];
			gauss_global[1] = local[0][1] * gauss_local[0] + local[1][1] * gauss_local[1] + p0[1];
			gauss_global[2] = local[0][2] * gauss_local[0] + local[1][2] * gauss_local[1] + p0[2];

			// 累加积分结果
			middle = operation(q, gauss_global);

			result[0] += middle[0] * Jacobi * weight;
			result[1] += middle[1] * Jacobi * weight;
			result[2] += middle[2] * Jacobi * weight;
		}
	}

	return result;
}





//2*2高斯积分重载（返回一个向量）
std::array<double, 3> Source_Dipole::Gauss_integrate(std::array<double, 3>(*operation)(std::array<double, 3>&, std::array<double, 3>&))
{
	//积分结果,整体坐标和局部坐标
	std::array<double, 3> result{};
	std::array<double, 3> result_local{};
	//雅可比矩阵行列式
	double Jacobi = 0;

	//2*2高斯点映射坐标
	std::array<double, 2> ksi = { 1.0 / sqrt(3),-1.0 / sqrt(3) };
	std::array<double, 2> eta = { 1.0 / sqrt(3),-1.0 / sqrt(3) };

	//高斯点局部坐标及整体坐标
	std::array<double, 2> gauss_local{};
	std::array<double, 3> gauss_global{};
	/*gauss_local[2] = 0.0;*/

	//存储operation的中间函数值
	std::array<double, 3> middle{};

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			Jacobi = ((-(1 - eta[j]) * p_local[0][0] + (1 - eta[j]) * p_local[1][0] + (1 + eta[j]) * p_local[2][0] - (1 + eta[j]) * p_local[3][0])
				* (-(1 - ksi[i]) * p_local[0][1] - (1 + ksi[i]) * p_local[1][1] + (1 + ksi[i]) * p_local[2][1] + (1 - ksi[i]) * p_local[3][1])
				- (-(1 - ksi[i]) * p_local[0][0] - (1 + ksi[i]) * p_local[1][0] + (1 + ksi[i]) * p_local[2][0] + (1 - ksi[i]) * p_local[3][0])
				* (-(1 - eta[j]) * p_local[0][1] + (1 - eta[j]) * p_local[1][1] + (1 + eta[j]) * p_local[2][1] - (1 + eta[j]) * p_local[3][1])) / 16.0;

			gauss_local[0] = ((1 - ksi[i]) * (1 - eta[j]) * p_local[0][0] + (1 + ksi[i]) * (1 - eta[j]) * p_local[1][0]
				+ (1 + ksi[i]) * (1 + eta[j]) * p_local[2][0] + (1 - ksi[i]) * (1 + eta[j]) * p_local[3][0]) / 4.0;
			gauss_local[1] = ((1 - ksi[i]) * (1 - eta[j]) * p_local[0][1] + (1 + ksi[i]) * (1 - eta[j]) * p_local[1][1]
				+ (1 + ksi[i]) * (1 + eta[j]) * p_local[2][1] + (1 - ksi[i]) * (1 + eta[j]) * p_local[3][1]) / 4.0;

			gauss_global[0] = local[0][0] * gauss_local[0] + local[1][0] * gauss_local[1] + p0[0];
			gauss_global[1] = local[0][1] * gauss_local[0] + local[1][1] * gauss_local[1] + p0[1];
			gauss_global[2] = local[0][2] * gauss_local[0] + local[1][2] * gauss_local[1] + p0[2];

			/*middle = operation(q_local, gauss_local);*/
			middle = operation(q, gauss_global);

			result_local[0] += middle[0] * Jacobi;
			result_local[1] += middle[1] * Jacobi;
			result_local[2] += middle[2] * Jacobi;
		}
	}

	//局部坐标转换到整体坐标
	/*result[0] = result_local[0] * local[0][0] + result_local[1] * local[1][0] + result_local[2] * local[2][0];
	result[1] = result_local[0] * local[0][1] + result_local[1] * local[1][1] + result_local[2] * local[2][1];
	result[2] = result_local[0] * local[0][2] + result_local[1] * local[1][2] + result_local[2] * local[2][2];*/

	return result_local;
}


std::array<double, 3> Source_Dipole::Gauss_integrate(std::array<double, 3>(*operation)(std::array<double, 3>&, std::array<double, 3>&, double), double delta_time)
{
	//积分结果,整体坐标和局部坐标
	std::array<double, 3> result{};
	std::array<double, 3> result_local{};
	//雅可比矩阵行列式
	double Jacobi = 0;

	//2*2高斯点映射坐标
	std::array<double, 2> ksi = { 1.0 / sqrt(3),-1.0 / sqrt(3) };
	std::array<double, 2> eta = { 1.0 / sqrt(3),-1.0 / sqrt(3) };

	//高斯点局部坐标及整体坐标
	std::array<double, 2> gauss_local{};
	std::array<double, 3> gauss_global{};
	/*gauss_local[2] = 0.0;*/

	//存储operation的中间函数值
	std::array<double, 3> middle{};

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			Jacobi = ((-(1 - eta[j]) * p_local[0][0] + (1 - eta[j]) * p_local[1][0] + (1 + eta[j]) * p_local[2][0] - (1 + eta[j]) * p_local[3][0])
				* (-(1 - ksi[i]) * p_local[0][1] - (1 + ksi[i]) * p_local[1][1] + (1 + ksi[i]) * p_local[2][1] + (1 - ksi[i]) * p_local[3][1])
				- (-(1 - ksi[i]) * p_local[0][0] - (1 + ksi[i]) * p_local[1][0] + (1 + ksi[i]) * p_local[2][0] + (1 - ksi[i]) * p_local[3][0])
				* (-(1 - eta[j]) * p_local[0][1] + (1 - eta[j]) * p_local[1][1] + (1 + eta[j]) * p_local[2][1] - (1 + eta[j]) * p_local[3][1])) / 16.0;

			gauss_local[0] = ((1 - ksi[i]) * (1 - eta[j]) * p_local[0][0] + (1 + ksi[i]) * (1 - eta[j]) * p_local[1][0]
				+ (1 + ksi[i]) * (1 + eta[j]) * p_local[2][0] + (1 - ksi[i]) * (1 + eta[j]) * p_local[3][0]) / 4.0;
			gauss_local[1] = ((1 - ksi[i]) * (1 - eta[j]) * p_local[0][1] + (1 + ksi[i]) * (1 - eta[j]) * p_local[1][1]
				+ (1 + ksi[i]) * (1 + eta[j]) * p_local[2][1] + (1 - ksi[i]) * (1 + eta[j]) * p_local[3][1]) / 4.0;

			gauss_global[0] = local[0][0] * gauss_local[0] + local[1][0] * gauss_local[1] + p0[0];
			gauss_global[1] = local[0][1] * gauss_local[0] + local[1][1] * gauss_local[1] + p0[1];
			gauss_global[2] = local[0][2] * gauss_local[0] + local[1][2] * gauss_local[1] + p0[2];

			middle = operation(q, gauss_global, delta_time);

			result_local[0] += middle[0] * Jacobi;
			result_local[1] += middle[1] * Jacobi;
			result_local[2] += middle[2] * Jacobi;
		}
	}

	//局部坐标转换到整体坐标
	/*result[0] = result_local[0] * local[0][0] + result_local[1] * local[1][0] + result_local[2] * local[2][0];
	result[1] = result_local[0] * local[0][1] + result_local[1] * local[1][1] + result_local[2] * local[2][1];
	result[2] = result_local[0] * local[0][2] + result_local[1] * local[1][2] + result_local[2] * local[2][2];*/

	return result_local;
}

TDGF_coeff Source_Dipole::Gauss_integrate(TDGF_coeff(*operation)(std::array<double, 3>&, std::array<double, 3>&, int, int), int l, int m)
{
	TDGF_coeff result{};
	TDGF_coeff result_local{};

	//雅可比矩阵行列式
	double Jacobi = 0;

	//2*2高斯点映射坐标
	std::array<double, 2> ksi = { 1.0 / sqrt(3),-1.0 / sqrt(3) };
	std::array<double, 2> eta = { 1.0 / sqrt(3),-1.0 / sqrt(3) };

	//高斯点局部坐标及整体坐标
	std::array<double, 2> gauss_local{};
	std::array<double, 3> gauss_global{};
	/*gauss_local[2] = 0.0;*/

	//存储operation的中间函数值
	TDGF_coeff middle{};

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			Jacobi = ((-(1 - eta[j]) * p_local[0][0] + (1 - eta[j]) * p_local[1][0] + (1 + eta[j]) * p_local[2][0] - (1 + eta[j]) * p_local[3][0])
				* (-(1 - ksi[i]) * p_local[0][1] - (1 + ksi[i]) * p_local[1][1] + (1 + ksi[i]) * p_local[2][1] + (1 - ksi[i]) * p_local[3][1])
				- (-(1 - ksi[i]) * p_local[0][0] - (1 + ksi[i]) * p_local[1][0] + (1 + ksi[i]) * p_local[2][0] + (1 - ksi[i]) * p_local[3][0])
				* (-(1 - eta[j]) * p_local[0][1] + (1 - eta[j]) * p_local[1][1] + (1 + eta[j]) * p_local[2][1] - (1 + eta[j]) * p_local[3][1])) / 16.0;

			gauss_local[0] = ((1 - ksi[i]) * (1 - eta[j]) * p_local[0][0] + (1 + ksi[i]) * (1 - eta[j]) * p_local[1][0]
				+ (1 + ksi[i]) * (1 + eta[j]) * p_local[2][0] + (1 - ksi[i]) * (1 + eta[j]) * p_local[3][0]) / 4.0;
			gauss_local[1] = ((1 - ksi[i]) * (1 - eta[j]) * p_local[0][1] + (1 + ksi[i]) * (1 - eta[j]) * p_local[1][1]
				+ (1 + ksi[i]) * (1 + eta[j]) * p_local[2][1] + (1 - ksi[i]) * (1 + eta[j]) * p_local[3][1]) / 4.0;

			gauss_global[0] = local[0][0] * gauss_local[0] + local[1][0] * gauss_local[1] + p0[0];
			gauss_global[1] = local[0][1] * gauss_local[0] + local[1][1] * gauss_local[1] + p0[1];
			gauss_global[2] = local[0][2] * gauss_local[0] + local[1][2] * gauss_local[1] + p0[2];

			/*middle = operation(q_local, gauss_local, l, m);*/
			middle = operation(q, gauss_global, l, m);

			result_local.T_A[0] += middle.T_A[0] * Jacobi;
			result_local.T_A[1] += middle.T_A[1] * Jacobi;
			result_local.T_A[2] += middle.T_A[2] * Jacobi;

			result_local.T_B += middle.T_B * Jacobi;
			/*result.T_B += middle.T_B * Jacobi;*/
		}
	}

	//局部坐标转换到整体坐标
	/*result.T_A[0] = result_local.T_A[0] * local[0][0] + result_local.T_A[1] * local[1][0] + result_local.T_A[2] * local[2][0];
	result.T_A[1] = result_local.T_A[0] * local[0][1] + result_local.T_A[1] * local[1][1] + result_local.T_A[2] * local[2][1];
	result.T_A[2] = result_local.T_A[0] * local[0][2] + result_local.T_A[1] * local[1][2] + result_local.T_A[2] * local[2][2];*/

	return result_local;
}

TDGF_coeff Source_Dipole::Gauss_integrate(TDGF_coeff(*operation)(std::array<double, 3>&, std::array<double, 3>&, int), int tN)
{
	TDGF_coeff result{};
	TDGF_coeff result_local{};

	//雅可比矩阵行列式
	double Jacobi = 0;

	//2*2高斯点映射坐标
	std::array<double, 2> ksi = { 1.0 / sqrt(3),-1.0 / sqrt(3) };
	std::array<double, 2> eta = { 1.0 / sqrt(3),-1.0 / sqrt(3) };

	//高斯点局部坐标及整体坐标
	std::array<double, 2> gauss_local{};
	std::array<double, 3> gauss_global{};
	/*gauss_local[2] = 0.0;*/

	//存储operation的中间函数值
	TDGF_coeff middle{};

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			Jacobi = ((-(1 - eta[j]) * p_local[0][0] + (1 - eta[j]) * p_local[1][0] + (1 + eta[j]) * p_local[2][0] - (1 + eta[j]) * p_local[3][0])
				* (-(1 - ksi[i]) * p_local[0][1] - (1 + ksi[i]) * p_local[1][1] + (1 + ksi[i]) * p_local[2][1] + (1 - ksi[i]) * p_local[3][1])
				- (-(1 - ksi[i]) * p_local[0][0] - (1 + ksi[i]) * p_local[1][0] + (1 + ksi[i]) * p_local[2][0] + (1 - ksi[i]) * p_local[3][0])
				* (-(1 - eta[j]) * p_local[0][1] + (1 - eta[j]) * p_local[1][1] + (1 + eta[j]) * p_local[2][1] - (1 + eta[j]) * p_local[3][1])) / 16.0;

			gauss_local[0] = ((1 - ksi[i]) * (1 - eta[j]) * p_local[0][0] + (1 + ksi[i]) * (1 - eta[j]) * p_local[1][0]
				+ (1 + ksi[i]) * (1 + eta[j]) * p_local[2][0] + (1 - ksi[i]) * (1 + eta[j]) * p_local[3][0]) / 4.0;
			gauss_local[1] = ((1 - ksi[i]) * (1 - eta[j]) * p_local[0][1] + (1 + ksi[i]) * (1 - eta[j]) * p_local[1][1]
				+ (1 + ksi[i]) * (1 + eta[j]) * p_local[2][1] + (1 - ksi[i]) * (1 + eta[j]) * p_local[3][1]) / 4.0;

			gauss_global[0] = local[0][0] * gauss_local[0] + local[1][0] * gauss_local[1] + p0[0];
			gauss_global[1] = local[0][1] * gauss_local[0] + local[1][1] * gauss_local[1] + p0[1];
			gauss_global[2] = local[0][2] * gauss_local[0] + local[1][2] * gauss_local[1] + p0[2];

			middle = operation(q, gauss_global, tN);


			result_local.T_A[0] += middle.T_A[0] * Jacobi;
			result_local.T_A[1] += middle.T_A[1] * Jacobi;
			result_local.T_A[2] += middle.T_A[2] * Jacobi;

			result_local.T_B += middle.T_B * Jacobi;
		}
	}

	//局部坐标转换到整体坐标
	/*result.T_A[0] = result_local.T_A[0] * local[0][0] + result_local.T_A[1] * local[1][0] + result_local.T_A[2] * local[2][0];
	result.T_A[1] = result_local.T_A[0] * local[0][1] + result_local.T_A[1] * local[1][1] + result_local.T_A[2] * local[2][1];
	result.T_A[2] = result_local.T_A[0] * local[0][2] + result_local.T_A[1] * local[1][2] + result_local.T_A[2] * local[2][2];*/

	return result_local;
}



TDGF_coeff Source_Dipole::Gauss_integrate44(TDGF_coeff(*operation)(std::array<double, 3>&, std::array<double, 3>&, int), int tN)
{
	TDGF_coeff result_local{};

	//雅可比矩阵行列式
	double Jacobi = 0;

	// 16×16高斯积分点的位置和权重
	const std::array<double, 16> gauss_points = {
		-0.9894009349916499,
		-0.9445750230732326,
		-0.8656312023878318,
		-0.7554044083550030,
		-0.6178762444026438,
		-0.4580167776572274,
		-0.2816035507792589,
		-0.0950125098376374,
		 0.0950125098376374,
		 0.2816035507792589,
		 0.4580167776572274,
		 0.6178762444026438,
		 0.7554044083550030,
		 0.8656312023878318,
		 0.9445750230732326,
		 0.9894009349916499
	};

	const std::array<double, 16> weights = {
		0.0271524594117541,
		0.0622535239386479,
		0.0951585116824928,
		0.1246289712555339,
		0.1495959888165767,
		0.1691565193950025,
		0.1826034150449236,
		0.1894506104550685,
		0.1894506104550685,
		0.1826034150449236,
		0.1691565193950025,
		0.1495959888165767,
		0.1246289712555339,
		0.0951585116824928,
		0.0622535239386479,
		0.0271524594117541
	};

	//存储operation的中间函数值
	TDGF_coeff middle{};

	// 遍历所有16×16=256个高斯点组合
	for (int i = 0; i < 16; ++i) {
		for (int j = 0; j < 16; ++j) {
			const double ksi = gauss_points[i];
			const double eta = gauss_points[j];
			const double weight = weights[i] * weights[j];

			// 计算雅可比行列式
			double Jacobi = ((-(1 - eta) * p_local[0][0] + (1 - eta) * p_local[1][0] + (1 + eta) * p_local[2][0] - (1 + eta) * p_local[3][0])
				* (-(1 - ksi) * p_local[0][1] - (1 + ksi) * p_local[1][1] + (1 + ksi) * p_local[2][1] + (1 - ksi) * p_local[3][1])
				- (-(1 - ksi) * p_local[0][0] - (1 + ksi) * p_local[1][0] + (1 + ksi) * p_local[2][0] + (1 - ksi) * p_local[3][0])
				* (-(1 - eta) * p_local[0][1] + (1 - eta) * p_local[1][1] + (1 + eta) * p_local[2][1] - (1 + eta) * p_local[3][1])) / 16.0;

			// 计算高斯点在局部坐标系中的坐标
			std::array<double, 2> gauss_local;
			gauss_local[0] = ((1 - ksi) * (1 - eta) * p_local[0][0] + (1 + ksi) * (1 - eta) * p_local[1][0]
				+ (1 + ksi) * (1 + eta) * p_local[2][0] + (1 - ksi) * (1 + eta) * p_local[3][0]) / 4.0;
			gauss_local[1] = ((1 - ksi) * (1 - eta) * p_local[0][1] + (1 + ksi) * (1 - eta) * p_local[1][1]
				+ (1 + ksi) * (1 + eta) * p_local[2][1] + (1 - ksi) * (1 + eta) * p_local[3][1]) / 4.0;

			// 转换到全局坐标系
			std::array<double, 3> gauss_global;
			gauss_global[0] = local[0][0] * gauss_local[0] + local[1][0] * gauss_local[1] + p0[0];
			gauss_global[1] = local[0][1] * gauss_local[0] + local[1][1] * gauss_local[1] + p0[1];
			gauss_global[2] = local[0][2] * gauss_local[0] + local[1][2] * gauss_local[1] + p0[2];

			// 累加积分结果
			middle = operation(q, gauss_global, tN);

			result_local.T_A[0] += middle.T_A[0] * Jacobi * weight;
			result_local.T_A[1] += middle.T_A[1] * Jacobi * weight;
			result_local.T_A[2] += middle.T_A[2] * Jacobi * weight;

			result_local.T_B += middle.T_B * Jacobi * weight;
		}
	}

	return result_local;
}



//分布源诱导速度（Rankine源）
std::array<double, 3> Source_Dipole::source_v()
{
	int j;

	for (int i = 0; i < 4; i++)
	{
		j = (i + 1) % 4;
		sv_local[0] += (-l_array[i][1] / l[i] * log((r[i] + r[j] + l[i]) / (r[i] + r[j] - l[i])));
		sv_local[1] += l_array[i][0] / l[i] * log((r[i] + r[j] + l[i]) / (r[i] + r[j] - l[i]));
	}

	if (sv3 == nullptr)
	{
		if (q_local[2] == 0)
		{
			sv_local[2] = 0;
		}
		else
		{
			double ci[2];
			double mi;
			double hi[2];
			double ri[2];
			for (int i = 0; i < 4; i++)
			{
				j = (i + 1) % 4;
				ci[0] = pow(r_array[i][0], 2) + pow(q_local[2], 2);
				ci[1] = pow(r_array[j][0], 2) + pow(q_local[2], 2);
				mi = l_array[i][1] / l_array[i][0];
				hi[0] = r_array[i][0] * r_array[i][1];
				hi[1] = r_array[j][0] * r_array[j][1];
				ri[0] = sqrt(ci[0] + pow(r_array[i][1], 2));
				ri[1] = sqrt(ci[1] + pow(r_array[j][1], 2));
				sv_local[2] += (atan((mi * ci[0] - hi[0]) / (q_local[2] * ri[0])) - atan((mi * ci[1] - hi[1]) / (q_local[2] * ri[1])));
			}
		}
		sv3 = &sv_local[2];
	}
	//局部坐标转换到整体坐标
	sv[0] = sv_local[0] * local[0][0] + sv_local[1] * local[1][0] + sv_local[2] * local[2][0];
	sv[1] = sv_local[0] * local[0][1] + sv_local[1] * local[1][1] + sv_local[2] * local[2][1];
	sv[2] = sv_local[0] * local[0][2] + sv_local[1] * local[1][2] + sv_local[2] * local[2][2];
	return sv;
}

//分布源（Rankine源）
double Source_Dipole::get_source()
{
	double ci[2];
	double mi;
	double hi[2];
	double ri[2];
	int m;
	if (sv3 == nullptr)
	{
		if (q_local[2] == 0)
		{
			sv_local[2] = 0;
		}
		else
		{
			for (int i = 0; i < 4; i++)
			{
				m = (i + 1) % 4;
				ci[0] = pow(r_array[i][0], 2) + pow(q_local[2], 2);
				ci[1] = pow(r_array[m][0], 2) + pow(q_local[2], 2);
				mi = l_array[i][1] / l_array[i][0];
				hi[0] = r_array[i][0] * r_array[i][1];
				hi[1] = r_array[m][0] * r_array[m][1];
				ri[0] = sqrt(ci[0] + pow(r_array[i][1], 2));
				ri[1] = sqrt(ci[1] + pow(r_array[m][1], 2));
				sv_local[2] += (atan((mi * ci[0] - hi[0]) / (q_local[2] * ri[0])) - atan((mi * ci[1] - hi[1]) / (q_local[2] * ri[1])));

				//std::cout << (mi * ci[0] - hi[0]) / (q_local[2] * ri[0]) << "," << (mi * ci[1] - hi[1]) / (q_local[2] * ri[1]) << std::endl;
			}
		}
		sv3 = &sv_local[2];
	}

	//单层位势
	source = 0;
	for (int j = 0; j < 4; j++)
	{
		m = (j + 1) % 4;
		ci[0] = pow(r_array[j][0], 2) + pow(q_local[2], 2);
		ci[1] = pow(r_array[m][0], 2) + pow(q_local[2], 2);
		ri[0] = sqrt(ci[0] + pow(r_array[j][1], 2));
		ri[1] = sqrt(ci[1] + pow(r_array[m][1], 2));
		source += (l_array[j][0] * r_array[j][1] - l_array[j][1] * r_array[j][0]) / l[j]
			* log((ri[0] + ri[1] + l[j]) / (ri[0] + ri[1] - l[j]));
	}
	source += q_local[2] * sv_local[2];

	return source;
}


//分布偶极诱导速度
std::array<double, 3> Source_Dipole::dipole_v()
{
	double ci[2];
	double ri[2];
	double dv_local[3] = { 0, 0, 0 };
	int j;
	for (int i = 0; i < 4; i++)
	{
		j = (i + 1) % 4;
		ci[0] = pow(r_array[i][0], 2) + pow(q_local[2], 2);
		ci[1] = pow(r_array[j][0], 2) + pow(q_local[2], 2);
		ri[0] = sqrt(ci[0] + pow(r_array[i][1], 2));
		ri[1] = sqrt(ci[1] + pow(r_array[j][1], 2));
		dv_local[0] += -2 * q_local[2] * l_array[i][1] * (ri[0] + ri[1]) / (ri[0] * ri[1] * (pow(ri[0] + ri[1], 2) - pow(l[i], 2)));
		dv_local[1] += 2 * q_local[2] * l_array[i][0] * (ri[0] + ri[1]) / (ri[0] * ri[1] * (pow(ri[0] + ri[1], 2) - pow(l[i], 2)));
		dv_local[2] += -2 * (l_array[i][0] * r_array[i][1] - l_array[i][1] * r_array[i][0]) * (ri[0] + ri[1])
			/ (ri[0] * ri[1] * (pow(ri[0] + ri[1], 2) - pow(l[i], 2)));
	}
	//局部坐标转换到整体坐标
	dv[0] = dv_local[0] * local[0][0] + dv_local[1] * local[1][0] + dv_local[2] * local[2][0];
	dv[1] = dv_local[0] * local[0][1] + dv_local[1] * local[1][1] + dv_local[2] * local[2][1];
	dv[2] = dv_local[0] * local[0][2] + dv_local[1] * local[1][2] + dv_local[2] * local[2][2];
	return dv;
}

//分布偶极
double Source_Dipole::get_dipole()
{
	if (q_local[2] == 0)
	{
		return 0;
	}

	if (sv3 == nullptr)
	{
		double ci[2];
		double mi;
		double hi[2];
		double ri[2];
		int j;
		for (int i = 0; i < 4; i++)
		{
			j = (i + 1) % 4;
			ci[0] = pow(r_array[i][0], 2) + pow(q_local[2], 2);
			ci[1] = pow(r_array[j][0], 2) + pow(q_local[2], 2);
			mi = l_array[i][1] / l_array[i][0];
			hi[0] = r_array[i][0] * r_array[i][1];
			hi[1] = r_array[j][0] * r_array[j][1];
			ri[0] = sqrt(ci[0] + pow(r_array[i][1], 2));
			ri[1] = sqrt(ci[1] + pow(r_array[j][1], 2));
			sv_local[2] += (atan((mi * ci[0] - hi[0]) / (q_local[2] * ri[0])) - atan((mi * ci[1] - hi[1]) / (q_local[2] * ri[1])));
		}
		sv3 = &sv_local[2];
	}
	dipole = -sv_local[2];
	return dipole;
}