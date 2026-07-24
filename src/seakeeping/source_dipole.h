#pragma once
#include <array>
#include <vector>

//四边形面元
class Mesh {
public:
	std::array<double, 3> point1;
	std::array<double, 3> point2;
	std::array<double, 3> point3;
	std::array<double, 3> point4;

	//中心点
	std::array<double, 3> center_p();

	//单位法向量
	std::array<double, 3> Nor_vector();

	//面积
	double area();
};

//面元的集合
struct AllMesh
{
	//所有面元
	std::vector<Mesh> Meshs;

	//所有面元的中心点
	std::vector<std::array<double, 3>> centers;

	//所有面元的单位法向量
	std::vector<std::array<double, 3>> Nor_vectors;

	//所有面元的面积
	std::vector <double> areas;

	AllMesh(int NE):
		Meshs(NE), centers(NE), Nor_vectors(NE), areas(NE)
	{}

};

//三维时域格林函数系数矩阵组合的结构体
struct TDGF_coeff
{
	std::array<double, 3> T_A;
	double T_B;
};


//四边形面元上均匀分布奇点（源、偶极），计算诱导速度
class Source_Dipole
{
public:
	//赋值、更新，整体坐标转换到局部坐标，局部坐标xy平面在四边形面元平面上，z轴为其法线方向，指向物体内
	void update(std::array<double, 3>& qq, Mesh& mesh);

	//2*2高斯积分，返回四边形面元的面积
	double Gauss_integrate();

	//2*2高斯积分（返回一个值）,函数指针第一个参数是场点整体坐标，第二个是源点整体坐标
	double Gauss_integrate(double(*operation)(std::array<double, 3>&, std::array<double, 3>&));

	//4*4高斯积分
	double Gauss_integrate44(double(*operation)(std::array<double, 3>&, std::array<double, 3>&));

	std::array<double, 3> Gauss_integrate44(std::array<double, 3>(*operation)(std::array<double, 3>&, std::array<double, 3>&));

	//2*2高斯积分重载（返回一个向量）
	std::array<double, 3> Gauss_integrate(std::array<double, 3>(*operation)(std::array<double, 3>&, std::array<double, 3>&));

	//2*2高斯积分重载（返回一个向量）（与时间有关）
	std::array<double, 3> Gauss_integrate(std::array<double, 3>(*operation)(std::array<double, 3>&, std::array<double, 3>&, double), double delta_time);

	//2*2高斯积分重载（返回三维时域格林函数瞬时项系数矩阵组合的结构体）
	// （与时间无关，仅需求解一次）(l==m和l!=m的取值不同，l==m代表对角线元素）
	TDGF_coeff Gauss_integrate(TDGF_coeff(*operation)(std::array<double, 3>&, std::array<double, 3>&, int, int), int l, int m);

	//2*2高斯积分重载（返回三维时域格林函数自由面记忆项项系数矩阵组合的结构体）（与时间有关，需在每个时间步计算）
	TDGF_coeff Gauss_integrate(TDGF_coeff(*operation)(std::array<double, 3>&, std::array<double, 3>&, int), int tN);

	TDGF_coeff Gauss_integrate44(TDGF_coeff(*operation)(std::array<double, 3>&, std::array<double, 3>&, int), int tN);

	//分布源诱导速度(Rankine源）
	std::array<double, 3> source_v();

	//分布偶极诱导速度(Rankine源）
	std::array<double, 3> dipole_v();

	//分布源(Rankine源）
	double get_source();

	//分布偶极(Rankine源）
	double get_dipole();

private:
	//整体坐标，场点和源点
	std::array<double, 3> q;
	Mesh* p;

	//局部坐标信息
	std::array<double, 3> p0;
	std::array<std::array<double, 3>, 3> local;

	//局部坐标下的场点和源点
	std::array<double, 3> q_local;
	std::array<std::array<double, 2>, 4> p_local;

	//四边形顶点到场点的向量和距离
	std::array<std::array<double, 3>, 4> r_array;
	std::array<double, 4> r;
	//四边形四条边的向量和距离
	std::array<std::array<double, 2>, 4> l_array;
	std::array<double, 4> l;

	//分布源、偶极及其诱导速度(Rankine源）（整体坐标）
	double source;
	double dipole;
	std::array<double, 3> sv_local;
	std::array<double, 3> sv;
	double* sv3;
	std::array<double, 3> dv;
};