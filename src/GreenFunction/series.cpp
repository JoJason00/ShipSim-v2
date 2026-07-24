#include "series.h"
#include <array>
#include <math.h>
#include <iostream>
#include <vector>

constexpr double PI = 3.14159265359;    //圆周率

constexpr double g = 9.8065;            //重力加速度

//第三区域第一部分中间函数g(k)
//static double G3_g(double miu, double k)
//{
//	double J0 = 1;
//	double c = 1.0;
//	double k2 = pow((pow(k, 2) / 2.0), 2);
//
//	int i = 1;
//	while (abs(c) > 1e-10)
//	{
//		c *= -k2 / pow(i, 2);
//		J0 += c;
//		i++;
//	}
//	double g = k * exp(-pow(k, 2) * miu / sqrt(1 - pow(miu, 2))) * (k * J0 - sqrt(2.0 / PI) * cos(pow(k, 2) - PI / 4.0));
//	return g;
//}


inline double G3_g(double miu, double k)
{
	double J0 = 1;
	double c1 = 1.0;
	double c2 = 0.0;

	double k2 = pow((k * k / 2.0), 2);

	int i = 1;

	while (abs(c1) > 1e-6)
	{
		c2 = c1 * (-k2 / (i * i));
		c1 = c2 * (-k2 / ((i + 1) * (i + 1)));

		J0 += (c1 + c2);

		i += 2;
	}
	double g = k * exp(-(k * k) * miu / sqrt(1 - (miu * miu))) * (k * J0 - sqrt(2.0 / PI) * cos((k * k) - PI / 4.0));
	return g;
}

//复数结构体
struct complex
{
	double real;
	double imag;
};


std::array<double, 3> TDGF_ba(double belta, double miu)
{
	if (miu < 0 || miu > 1 || belta < 0)
	{
		std::cerr << "miu或belta超出取值范围" << std::endl;
		return { 0, 0, 0 };
	}

	//Green函数及其导函数
	std::array<double, 5> Green{};

	//无量纲化格林函数及其导函数
	double Green_bar = 0;
	double Green_belta = 0;
	double Green_miu = 0;

	if (miu >= 0 && miu < 0.7)
	{
		//第一区域
		if (belta >= 0 && belta < 10.0)
		{
			//每一项前可迭代修改的系数（乘上一个数），进而可减少重复计算
			double cn = 2 * belta;
			double cn_belta = 2;
			/*double cn_miu = belta;*/

			//legendre多项式及其导数，每次只存储三项，下标对应循环中真实下标取余的结果
			std::array<double, 3> Pn = { 1.0, miu, 0 };
			std::array<double, 3> Pn_miu = { 0.0, 1.0, 0 };

			//格林函数及其导数循环加上的项,这里每次存储两项
			std::array<double, 2> delta_bar = { 1.0, cn * Pn.at(1) };
			std::array<double, 2> delta_belta = { 1.0, cn_belta * Pn.at(1) };
			std::array<double, 2> delta_miu = { 1.0, cn * Pn_miu.at(1) };

			Green_bar += delta_bar.at(1);
			Green_belta += delta_belta.at(1);
			Green_miu += delta_miu.at(1);

			int i = 1;

			//计算精度设置为10^-10
			while (delta_bar.at(0) * delta_bar.at(0) + delta_bar.at(1) * delta_bar.at(1) > 1e-10 ||
				delta_belta.at(0) * delta_belta.at(0) + delta_belta.at(1) * delta_belta.at(1) > 1e-10 ||
				delta_miu.at(0) * delta_miu.at(0) + delta_miu.at(1) * delta_miu.at(1) > 1e-10)
			{
				i++;

				Pn.at(i % 3) = ((2 * i - 1) * miu * Pn.at((i - 1) % 3) - (i - 1) * Pn.at((i - 2) % 3)) / i;
				Pn_miu.at(i % 3) = ((2 * i - 1) * miu * Pn_miu.at((i - 1) % 3) + (2 * i - 1) * Pn.at((i - 1) % 3)
					- (i - 1) * Pn_miu.at((i - 2) % 3)) / i;

				cn *= -i / ((2.0 * i - 1.0) * (2.0 * i - 2.0)) * belta * belta;
				delta_bar.at(i % 2) = cn * Pn.at(i % 3);
				Green_bar += delta_bar.at(i % 2);

				cn_belta *= -i / ((2.0 * i - 2.0) * (2.0 * i - 3.0)) * belta * belta;
				delta_belta.at(i % 2) = cn_belta * Pn.at(i % 3);
				Green_belta += delta_belta.at(i % 2);

				/*cn_miu *= -i / ((2.0 * i - 1.0) * (2.0 * i - 2.0)) * pow(belta, 2);*/
				delta_miu.at(i % 2) = cn * Pn_miu.at(i % 3);
				Green_miu += delta_miu.at(i % 2);
			}
		}

		//第二区域，格林函数分两部分
		else
		{
			//预先求出多次用到的中间变量
			//belta的次方
			std::array<double, 6> belta_pow;
			belta_pow[0] = belta * belta;
			belta_pow[1] = belta_pow[0] * belta;
			belta_pow[2] = belta_pow[1] * belta;
			belta_pow[3] = belta_pow[2] * belta;
			belta_pow[4] = belta_pow[3] * belta;
			belta_pow[5] = belta_pow[4] * belta;

			//legendre多项式及其导数，每次只存储三项，下标对应循环中真实下标取余的结果
			std::array<double, 3> Pn = { 1.0, miu, 0.0 };
			std::array<double, 3> Pn_miu = { 0.0, 1.0, 0.0 };

			double cn = 8.0 / belta_pow[1];

			std::array<double, 2> delta_bar = { 1.0, -cn * Pn.at(0) };
			std::array<double, 2> delta_belta = { 1.0, 24.0 / belta_pow[2] * Pn.at(0) };
			std::array<double, 2> delta_miu = { 1.0, 0.0 };

			Green_bar += delta_bar.at(1);
			Green_belta += delta_belta.at(1);
			Green_miu += delta_miu.at(1);

			int i = 1;

			//Green_bar = 2.0 * (-4.0 / belta_pow[1] - 48.0 * miu / belta_pow[3] - 360.0 / belta_pow[5] * (3.0 * pow(miu, 2) - 1));

			//格林函数第一部分及其导数
			while (delta_bar.at(0) * delta_bar.at(0) + delta_bar.at(1) * delta_bar.at(1) > 1e-12 ||
				delta_belta.at(0) * delta_belta.at(0) + delta_belta.at(1) * delta_belta.at(1) > 1e-12 ||
				delta_miu.at(0) * delta_miu.at(0) + delta_miu.at(1) * delta_miu.at(1) > 1e-12)
			{
				i++;

				Pn.at(i % 3) = ((2 * i - 1) * miu * Pn.at((i - 1) % 3) - (i - 1) * Pn.at((i - 2) % 3)) / i;

				Pn_miu.at(i % 3) = ((2 * i - 1) * miu * Pn_miu.at((i - 1) % 3) + (2 * i - 1) * Pn.at((i - 1) % 3)
					- (i - 1) * Pn_miu.at((i - 2) % 3)) / i;

				cn *= (((2.0 * i - 1.0) * 2.0 / belta_pow[0]));

				delta_bar.at(i % 2) = -cn * i * Pn.at((i - 1) % 3);

				Green_bar += delta_bar.at(i % 2);

				delta_belta.at(i % 2) = cn * i * (2 * i + 1) / belta * Pn.at((i - 1) % 3);

				Green_belta += delta_belta.at(i % 2);

				delta_miu.at(i % 2) = -cn * i * Pn_miu.at((i - 1) % 3);
				Green_miu += delta_miu.at(i % 2);
			}

			//格林函数第二部分

			//预先求出多次用到的中间变量

			double theta = asin(miu);

			double c0 = sqrt(2.0) * exp(-belta_pow.at(0) * miu / 4.0);

			double c1 = 1.0 - miu * miu;

			//1-miu^2的次方
			std::array<double, 11> miu_pow =
			{ pow(c1, 0.25), sqrt(c1), pow(c1, 0.75), pow(c1, 1.25), pow(c1, 1.5),
			pow(c1, 2.25), pow(c1, 3.0 / 4.0), pow(c1, 7.0 / 3.0), pow(c1, 5.0 / 6.0), pow(c1, 11.0 / 6.0) ,pow(c1, 1.75) };

			double c2 = belta_pow[0] * miu_pow[1] / 4.0;

			//三角函数
			std::array<double, 10> trig =
			{ sin(c2 + 1.5 * theta), cos(c2 + 1.5 * theta), sin(c2 - 0.5 * theta), cos(c2 - 0.5 * theta),
			sin(c2 - theta * 1.5), cos(c2 - theta * 1.5), sin(c2 - 2.5 * theta), cos(c2 - 2.5 * theta),
			sin(c2 - 3.5 * theta), cos(c2 - 3.5 * theta) };

			//格林函数第二部分及其导数
			std::array<double, 3> G_part2;

			G_part2[0] = c0 * (trig[0] * belta / miu_pow[0] + trig[3] / (2.0 * belta * miu_pow[2])
				+ trig[4] / (miu_pow[2] * belta_pow[1]) - trig[6] * 9.0 / (8.0 * belta_pow[1] * miu_pow[3])
				- trig[9] * 9.0 / (belta_pow[3] * miu_pow[2]) + trig[7] * 24.0 / (belta_pow[3] * miu_pow[2]));

			G_part2[1] = c0 * (trig[0] / miu_pow[0] + trig[1] * belta_pow[0] * miu_pow[0] / 2.0 - trig[3] / (belta_pow[0] * miu_pow[2] * 2)
				- trig[2] / (4 * miu_pow[0]) - trig[4] * 3 / (belta_pow[2] * miu_pow[2]) + trig[5] / (2 * belta_pow[0] * miu_pow[0])
				+ trig[6] * 27 / (8 * belta_pow[2] * miu_pow[3]) - trig[7] * 9 / (16 * belta_pow[0] * miu_pow[2]) + trig[9] * 45 / (belta_pow[4] * miu_pow[6])
				+ trig[8] * 9 / (2 * belta_pow[2] * miu_pow[8]) - trig[7] * 120 / (belta_pow[4] * miu_pow[6]) - trig[6] * 12 / (belta_pow[2] * miu_pow[8]))
				- miu * belta / 2.0 * G_part2[0];

			G_part2[2] = c0 * (trig[0] * belta * miu / (2 * miu_pow[3]) + trig[1] * belta * (-miu * belta_pow[0] + 6) / (4 * miu_pow[2])
				+ trig[3] * 3 * miu / (4 * belta * miu_pow[10]) + trig[2] * (miu * belta_pow[0] + 2) / (8 * belta * miu_pow[3])
				+ trig[4] * 3 * miu / (2 * belta_pow[1] * miu_pow[10]) - trig[5] * (miu * belta_pow[0] + 6) / (4 * belta_pow[1] * miu_pow[3])
				- trig[6] * 45 * miu / (16 * belta_pow[1] * miu_pow[5]) + trig[7] * 9 * (miu * belta_pow[0] + 10) / (32 * belta_pow[1] * miu_pow[10])
				- trig[9] * 24 * miu / (belta_pow[3] * miu_pow[7]) - trig[8] * 9 * (miu * belta_pow[1] + 14) / (4 * belta_pow[3] * miu_pow[9])
				+ trig[7] * 64 * miu / (belta_pow[3] * miu_pow[7]) + trig[6] * 6 * (miu * belta_pow[1] + 10) / (belta_pow[3] * miu_pow[9]))
				- belta_pow.at(0) / 4.0 * G_part2[0];

			Green_bar += G_part2[0];
			Green_belta += G_part2[1];
			Green_miu += G_part2[2];
		}
	}

	else if (miu >= 0.7 && miu < 0.98 && belta >= 10)
	{
		//预先求出多次用到的中间变量
			//belta的次方
		std::array<double, 6> belta_pow;
		belta_pow[0] = belta * belta;
		belta_pow[1] = belta_pow[0] * belta;
		belta_pow[2] = belta_pow[1] * belta;
		belta_pow[3] = belta_pow[2] * belta;
		belta_pow[4] = belta_pow[3] * belta;
		belta_pow[5] = belta_pow[4] * belta;

		//legendre多项式及其导数，每次只存储三项，下标对应循环中真实下标取余的结果
		std::array<double, 3> Pn = { 1.0, miu, 0.0 };
		std::array<double, 3> Pn_miu = { 0.0, 1.0, 0.0 };

		double cn = 8.0 / belta_pow[1];

		std::array<double, 2> delta_bar = { 1.0, -cn * Pn.at(0) };
		std::array<double, 2> delta_belta = { 1.0, 24.0 / belta_pow[2] * Pn.at(0) };
		std::array<double, 2> delta_miu = { 1.0, 0.0 };

		Green_bar += delta_bar.at(1);
		Green_belta += delta_belta.at(1);
		Green_miu += delta_miu.at(1);

		int i = 1;

		//Green_bar = 2.0 * (-4.0 / belta_pow[1] - 48.0 * miu / belta_pow[3] - 360.0 / belta_pow[5] * (3.0 * pow(miu, 2) - 1));

		//格林函数第一部分及其导数
		while (delta_bar.at(0) * delta_bar.at(0) + delta_bar.at(1) * delta_bar.at(1) > 1e-12 ||
			delta_belta.at(0) * delta_belta.at(0) + delta_belta.at(1) * delta_belta.at(1) > 1e-12 ||
			delta_miu.at(0) * delta_miu.at(0) + delta_miu.at(1) * delta_miu.at(1) > 1e-12)
		{
			i++;

			Pn.at(i % 3) = ((2 * i - 1) * miu * Pn.at((i - 1) % 3) - (i - 1) * Pn.at((i - 2) % 3)) / i;

			Pn_miu.at(i % 3) = ((2 * i - 1) * miu * Pn_miu.at((i - 1) % 3) + (2 * i - 1) * Pn.at((i - 1) % 3)
				- (i - 1) * Pn_miu.at((i - 2) % 3)) / i;

			cn *= (((2.0 * i - 1.0) * 2.0 / belta_pow[0]));

			delta_bar.at(i % 2) = -cn * i * Pn.at((i - 1) % 3);

			Green_bar += delta_bar.at(i % 2);

			delta_belta.at(i % 2) = cn * i * (2 * i + 1) / belta * Pn.at((i - 1) % 3);

			Green_belta += delta_belta.at(i % 2);

			delta_miu.at(i % 2) = -cn * i * Pn_miu.at((i - 1) % 3);
			Green_miu += delta_miu.at(i % 2);
		}

		//格林函数第二部分

		//预先求出多次用到的中间变量

		double theta = asin(miu);

		double c0 = sqrt(2.0) * exp(-belta_pow.at(0) * miu / 4.0);

		double c1 = 1.0 - miu * miu;

		//1-miu^2的次方
		std::array<double, 11> miu_pow =
		{ pow(c1, 0.25), sqrt(c1), pow(c1, 0.75), pow(c1, 1.25), pow(c1, 1.5),
		pow(c1, 2.25), pow(c1, 3.0 / 4.0), pow(c1, 7.0 / 3.0), pow(c1, 5.0 / 6.0), pow(c1, 11.0 / 6.0) ,pow(c1, 1.75) };

		double c2 = belta_pow[0] * miu_pow[1] / 4.0;

		//三角函数
		std::array<double, 10> trig =
		{ sin(c2 + 1.5 * theta), cos(c2 + 1.5 * theta), sin(c2 - 0.5 * theta), cos(c2 - 0.5 * theta),
		sin(c2 - theta * 1.5), cos(c2 - theta * 1.5), sin(c2 - 2.5 * theta), cos(c2 - 2.5 * theta),
		sin(c2 - 3.5 * theta), cos(c2 - 3.5 * theta) };

		//格林函数第二部分及其导数
		std::array<double, 3> G_part2;

		G_part2[0] = c0 * (trig[0] * belta / miu_pow[0] + trig[3] / (2.0 * belta * miu_pow[2])
			+ trig[4] / (miu_pow[2] * belta_pow[1]) - trig[6] * 9.0 / (8.0 * belta_pow[1] * miu_pow[3])
			- trig[9] * 9.0 / (belta_pow[3] * miu_pow[2]) + trig[7] * 24.0 / (belta_pow[3] * miu_pow[2]));

		G_part2[1] = c0 * (trig[0] / miu_pow[0] + trig[1] * belta_pow[0] * miu_pow[0] / 2.0 - trig[3] / (belta_pow[0] * miu_pow[2] * 2)
			- trig[2] / (4 * miu_pow[0]) - trig[4] * 3 / (belta_pow[2] * miu_pow[2]) + trig[5] / (2 * belta_pow[0] * miu_pow[0])
			+ trig[6] * 27 / (8 * belta_pow[2] * miu_pow[3]) - trig[7] * 9 / (16 * belta_pow[0] * miu_pow[2]) + trig[9] * 45 / (belta_pow[4] * miu_pow[6])
			+ trig[8] * 9 / (2 * belta_pow[2] * miu_pow[8]) - trig[7] * 120 / (belta_pow[4] * miu_pow[6]) - trig[6] * 12 / (belta_pow[2] * miu_pow[8]))
			- miu * belta / 2.0 * G_part2[0];

		G_part2[2] = c0 * (trig[0] * belta * miu / (2 * miu_pow[3]) + trig[1] * belta * (-miu * belta_pow[0] + 6) / (4 * miu_pow[2])
			+ trig[3] * 3 * miu / (4 * belta * miu_pow[10]) + trig[2] * (miu * belta_pow[0] + 2) / (8 * belta * miu_pow[3])
			+ trig[4] * 3 * miu / (2 * belta_pow[1] * miu_pow[10]) - trig[5] * (miu * belta_pow[0] + 6) / (4 * belta_pow[1] * miu_pow[3])
			- trig[6] * 45 * miu / (16 * belta_pow[1] * miu_pow[5]) + trig[7] * 9 * (miu * belta_pow[0] + 10) / (32 * belta_pow[1] * miu_pow[10])
			- trig[9] * 24 * miu / (belta_pow[3] * miu_pow[7]) - trig[8] * 9 * (miu * belta_pow[1] + 14) / (4 * belta_pow[3] * miu_pow[9])
			+ trig[7] * 64 * miu / (belta_pow[3] * miu_pow[7]) + trig[6] * 6 * (miu * belta_pow[1] + 10) / (belta_pow[3] * miu_pow[9]))
			- belta_pow.at(0) / 4.0 * G_part2[0];

		Green_bar += G_part2[0];
		Green_belta += G_part2[1];
		Green_miu += G_part2[2];
	}


	//第三区域，格林函数分两部分
	else if (miu >= 0.7 && miu < 0.98 && belta < 10)
		//else if (miu >= 0.7 && miu < 0.98 )
	{
		//格林函数第一部分

		//提前存储函数G3_g的结果
		std::vector<double> G3g_middle;
		//G3g_middle.reserve(100);

		double kmax = 0.05;
		double g_kmax = G3_g(miu, kmax);
		while (abs(g_kmax) >= 1e-7)
		{
			G3g_middle.push_back(g_kmax);
			kmax += 0.05;
			g_kmax = G3_g(miu, kmax);
		}
		kmax -= 0.05;

		double h = 0.05;         //步长
		int n = (int)(kmax / (2 * h));

		//std::cout << "n:" << n << "    " << "size:" << G3g_middle.size() << std::endl;

		//存储中间变量1-miu^2的次方
		double c_miu = 1 - miu * miu;
		std::array<double, 5> miu_pow = { pow(c_miu,0.25),sqrt(c_miu),pow(c_miu,0.75),pow(c_miu,1.25),pow(c_miu,1.5) };

		double gama = belta / miu_pow[0];
		double delta = gama * h;
		double delta_min = 0.35;
		double delta_belta = h / miu_pow[0];
		double delta_miu = belta * miu * h / (2.0 * miu_pow[3]);

		double S2n = 0;
		double S2n_belta = 0;
		double S2n_miu = 0;
		double S2n1 = 0;
		double S2n1_belta = 0;
		double S2n1_miu = 0;

		double g2i;                //中间变量g(2*i*h)
		double g2i1;
		for (int i = 1; i < n + 1; i++)
		{
			//g2i = G3_g(miu, 2 * i * h);
			//g2i1 = G3_g(miu, (2 * i - 1) * h);
			g2i = G3g_middle[2 * i - 1];
			g2i1 = G3g_middle[2 * i - 1 - 1];

			S2n += g2i * sin(2 * i * delta);
			S2n_belta += g2i * 2 * i * h * cos(2 * i * delta) / miu_pow[0];
			S2n_miu += (-pow(2 * i * h, 2) * g2i / miu_pow[4] * sin(2 * i * delta)
				+ g2i * belta * miu * i * h / miu_pow[3] * cos(2 * i * delta));

			S2n1 += g2i1 * sin((2 * i - 1) * delta);
			S2n1_belta += g2i1 * (2 * i - 1) * h * cos((2 * i - 1) * delta) / miu_pow[0];
			S2n1_miu += (-pow((2 * i - 1) * h, 2) * g2i1 / miu_pow[4] * sin((2 * i - 1) * delta)
				+ g2i1 * belta * miu * (2 * i - 1) * h / (2 * miu_pow[3]) * cos((2 * i - 1) * delta));
		}

		double a2;
		double a2_delta = 0;        //a2对delta的导数
		double a3;
		double a3_delta = 0;         //a3对delta的导数

		//存储delta的次方，方便调用
		double pow_delta[2] = { delta * delta, delta * delta * delta };

		if (delta > delta_min)
		{
			a2 = (delta * (3 + cos(2 * delta)) - 2 * sin(2 * delta)) / pow_delta[1];
			a2_delta = -3.0 * a2 / delta + (3 - 3 * cos(2 * delta) - 2 * delta * sin(2 * delta)) / pow_delta[1];
			a3 = 4 * (sin(delta) - delta * cos(delta)) / pow_delta[1];
			a3_delta = (-3 * delta * a3 + 4 * sin(delta)) / pow_delta[0];
		}
		else
		{
			a2 = 2.0 / 3.0;
			a3 = 4.0 / 3.0;
			double delta2 = pow_delta[0];

			double c2 = 16.0 / 24.0 * delta2;
			double c2_m = c2 * (1 - 4.0 / 5.0);
			double c2_delta = 16.0 * 2.0 / 24.0 * delta;
			double c2_deltam = c2_delta * (1 - 4.0 / 5.0);

			double c3 = -4.0 / 24.0 * delta2;
			double c3_m = c3 * (1 - 1.0 / 5.0);
			double c3_delta = -8.0 / 24.0 * delta;
			double c3_deltam = c3_delta * (1 - 1.0 / 5.0);
			int i = 1;
			while (abs(c2_m) > 1e-10 || abs(c3_m) > 1e-10 || abs(c2_deltam) > 1e-10 || abs(c3_deltam) > 1e-10)
			{
				a2 += c2_m;
				a3 += c3_m;
				a2_delta += c2_deltam;
				a3_delta += c3_deltam;

				i++;
				c2 *= -4.0 / ((2 * i + 2) * (2 * i + 1)) * delta2;
				c2_m = c2 * (1 - 4.0 / (2.0 * i + 3.0));
				c3 *= -1.0 / ((2 * i + 2) * (2 * i + 1)) * delta2;
				c3_m = c3 * (1 - 1.0 / (2.0 * i + 3.0));
				c2_delta *= -4.0 / ((2 * i + 2) * (2 * i + 1)) * delta2;
				c2_deltam = c2_delta * (1 - 4.0 / (2.0 * i + 3.0)) * i;
				c3_delta *= -1.0 / ((2 * i + 2) * (2 * i + 1)) * delta2;
				c3_deltam = c3_delta * (1 - 1.0 / (2.0 * i + 3.0)) * i;
			}
		}

		double green_c = 4.0 * h / miu_pow[2];

		Green_bar += green_c * (a2 * S2n + a3 * S2n1);
		Green_belta += green_c * (a2_delta * delta_belta * S2n + a2 * S2n_belta
			+ a3_delta * delta_belta * S2n1 + a3 * S2n1_belta);
		Green_miu += 1.5 * miu * Green_bar / c_miu + green_c * (a2_delta * delta_miu * S2n
			+ a2 * S2n_miu + a3_delta * delta_miu * S2n1 + a3 * S2n1_miu);

		//格林函数第二部分
		std::array<double, 3> G_part2{};

		double belta_pow = pow(belta, 2);

		double theta = asin(miu);
		double c0 = sqrt(2) * exp(-belta_pow * miu / 4.0);
		double c0_belta = -miu * belta / 2.0;
		double c0_miu = -belta_pow / 4.0;

		double c_trig = belta_pow * miu_pow[1] / 4.0;

		//三角函数
		std::array<double, 2> trig = { sin(c_trig + 1.5 * theta), cos(c_trig + 1.5 * theta) };

		G_part2[0] = c0 * trig[0] * belta / miu_pow[0];
		G_part2[1] = c0 * (trig[0] / miu_pow[0] + trig[1] * belta_pow * miu_pow[0] / 2.0) + c0_belta * G_part2[0];
		G_part2[2] = c0 * (trig[0] * belta * miu / (2 * miu_pow[3]) + trig[1] * belta * (-miu * belta_pow + 6) / (4 * miu_pow[2]))
			+ c0_miu * G_part2[0];

		Green_bar += G_part2[0];
		Green_belta += G_part2[1];
		Green_miu += G_part2[2];

	}

	else if (miu >= 0.98 && miu <= 1)
	{
		//第四区域
		if (belta < 10)
		{
			double belta_pow = belta * belta;
			double miu_pow = miu * miu;

			std::array<complex, 23> In{};
			In[0].real = sqrt(PI / miu) / 2.0 * exp(-belta_pow / (4 * miu));
			//In[0].imag = In[0].real * 2.0 / sqrt(PI) * romberg_integration(e_function, 0.0, belta / (sqrt(miu) * 2));
			//In[0].imag = In[0].real * 2.0 / sqrt(PI) * trapezoidal_rule(e_function, 0.0, belta / (sqrt(miu) * 2),1000);

			//误差函数计算
			double z = belta / (sqrt(miu) * 2);
			double omiga_imag = 0.0;
			double c_omiga = z;
			double delta_omiga = z;
			int omiga_i = 1;
			while (abs(delta_omiga) > 1e-8)
			{
				omiga_imag += delta_omiga;
				c_omiga *= (z * z) / omiga_i;
				delta_omiga = c_omiga / (2 * omiga_i + 1);
				omiga_i++;
			}
			In[0].imag = In[0].real * 2.0 / sqrt(PI) * omiga_imag;

			In[1].real = 1.0 / (2.0 * miu) - belta / (2.0 * miu) * In[0].imag;
			In[1].imag = belta / (2.0 * miu) * In[0].real;

			std::array<complex, 23> In_belta{};
			In_belta[0].real = -belta / (2.0 * miu) * In[0].real;
			In_belta[0].imag = -belta / (2.0 * miu) * In[0].imag + 0.5 / miu;
			In_belta[1].real = -In[0].imag / (2 * miu) - belta * In_belta[0].imag / (2 * miu);
			In_belta[1].imag = In[0].real / (2 * miu) + belta * In_belta[0].real / (2 * miu);

			std::array<complex, 23> In_miu{};
			In_miu[0].real = In[0].real * (-1 + belta_pow / (2 * miu)) / (2 * miu);
			In_miu[0].imag = In[0].imag * (-1 + belta_pow / (2 * miu)) / (2 * miu) - belta / (4.0 * miu_pow);
			In_miu[1].real = In[0].imag * belta / (2 * miu_pow) - In_miu[0].imag * belta / (2 * miu) - 1.0 / (2 * miu_pow);
			In_miu[1].imag = -In[0].real * belta / (2 * miu_pow) + In_miu[0].real * belta / (2 * miu);

			for (int i = 2; i < 23; i++)
			{
				In[i].real = (i - 1) / (2.0 * miu) * In[i - 2].real - belta / (2.0 * miu) * In[i - 1].imag;
				In[i].imag = (i - 1) / (2.0 * miu) * In[i - 2].imag + belta / (2.0 * miu) * In[i - 1].real;

				In_belta[i].real = (i - 1) / (2 * miu) * In_belta[i - 2].real - In[i - 1].imag / (2 * miu) - In_belta[i - 1].imag * belta / (2 * miu);
				In_belta[i].imag = (i - 1) / (2 * miu) * In_belta[i - 2].imag + In[i - 1].real / (2 * miu) + In_belta[i - 1].real * belta / (2 * miu);

				In_miu[i].real = -(i - 1) * In[i - 2].real / (2 * miu_pow) + (i - 1) * In_miu[i - 2].real / (2 * miu)
					+ In[i - 1].imag * belta / (2 * miu_pow) - In_miu[i - 1].imag * belta / (2 * miu);
				In_miu[i].imag = -(i - 1) * In[i - 2].imag / (2 * miu_pow) + (i - 1) * In_miu[i - 2].imag / (2 * miu)
					- In[i - 1].real * belta / (2 * miu_pow) + In_miu[i - 1].real * belta / (2 * miu);
			}

			if (1 - miu < 1e-5)
			{
				Green_bar = 4 * In[2].imag;
				Green_belta = 4 * In_belta[2].imag;
				Green_miu = 8.0 * 2.2499997 * In[6].imag / 9.0 + 4 * In_miu[2].imag;
			}
			else
			{
				std::array<double, 6> a = { 1.0, -2.2499997, 1.2656208, -0.3163866, 0.0444479, 0.0002100 };
				double c_miu = 1;

				for (int i = 0; i < 6; i++)
				{
					Green_bar += a[i] * c_miu * In[4 * i + 2].imag;
					Green_belta += a[i] * c_miu * In_belta[4 * i + 2].imag;
					Green_miu += a[i] * c_miu * (-2 * i * miu / (1 - miu_pow)) * In[4 * i + 2].imag + a[i] * c_miu * In_miu[4 * i + 2].imag;

					c_miu *= (1 - miu_pow) / 9.0;
				}
				Green_bar *= 4;
				Green_belta *= 4;
				Green_miu *= 4;
			}
		}

		//第五区域
		else
		{
			double belta_pow = belta * belta;
			double miu_pow = miu * miu;

			//计算S2及其导数
			double S2 = 0;
			double S2_belta = 0;
			double S2_miu = 0;
			double c_S2 = 1.0 / (belta_pow * belta);

			std::array<double, 3> delta_S2 = { c_S2 * 2, -c_S2 * 6 / belta, 0.0 };
			int n_count = 1;

			while (abs(delta_S2.at(0)) > 1e-7 || abs(delta_S2.at(1)) > 1e-7 || abs(delta_S2.at(2)) > 1e-7)
			{
				S2 += delta_S2.at(0);
				S2_belta += delta_S2.at(1);
				S2_miu += delta_S2.at(2);

				n_count++;
				c_S2 *= (2 * n_count - 1) * 2 * miu / belta_pow;
				delta_S2 = { c_S2 * n_count * 2, -c_S2 * n_count * 2 * (2 * n_count - 1) / belta, c_S2 * 2 * n_count * (n_count - 1) / miu };
			}

			//计算S6及其导数
			double S6 = 0;
			double S6_belta = 0;
			double S6_miu = 0;

			double c_S6 = 1.0 / pow(belta, 7);
			double cn_S6 = 240;             //16 * 5!!
			double S6_sumi = 1;

			std::array<double, 3> delta_S6 = { S6_sumi * c_S6 * cn_S6 * 3, -7 * S6_sumi * c_S6 * cn_S6 * 3 / belta, 0.0 };

			n_count = 1;
			while (abs(delta_S6.at(0)) > 1e-6 || abs(delta_S6.at(1)) > 1e-6 || abs(delta_S6.at(2)) > 1e-6)
			{
				n_count++;
				c_S6 *= 2 * miu / belta_pow;
				cn_S6 *= (2 * n_count + 3);
				S6_sumi += n_count;

				delta_S6 = { S6_sumi * c_S6 * cn_S6 * (2 + n_count),-(2 * n_count + 5) / belta * S6_sumi * c_S6 * cn_S6 * (2 + n_count)
					, S6_sumi * (n_count - 1) / miu * c_S6 * cn_S6 * (2 + n_count) };

				S6 += delta_S6.at(0);
				S6_belta += delta_S6.at(1);
				S6_miu += delta_S6.at(2);
			}

			double a1 = 1.0;
			double a2 = -2.2499997;
			if (1 - miu < 1e-5)           //miu约等于1
			{
				Green_bar = -4 * S2;
				Green_belta = -4 * S2_belta;
				Green_miu = -4 * S2_miu + 8.0 / 9.0 * a1 * S6;
			}
			else
			{
				//计算S10及其导数
				double S10 = 0;
				double S10_belta = 0;
				double S10_miu = 0;

				//double c_S10 = 1.0 / pow(belta, 11);
				double cn_S10 = 120960 / pow(belta, 11);  //128 * 9!! * ( 2 * miu)^(n-1) /belta^(2n-9)
				double S10_sum1 = 1;
				double S10_sum2 = 6;  //1+2+3

				std::array<double, 3> delta_S10 = { S10_sum1 * S10_sum2 * cn_S10 * 5, -11 * S10_sum1 * S10_sum2 * cn_S10 * 5 / belta, 0.0 };

				n_count = 1;
				while (abs(delta_S10.at(0)) > 1e-4 || abs(delta_S10.at(1)) > 1e-4 || abs(delta_S10.at(2)) > 1e-4)
				{
					n_count++;
					//c_S10 *= 2 * miu / pow(belta, 2);
					cn_S10 *= ((2 * n_count + 7) * 2 * miu / belta_pow);
					S10_sum1 += n_count;
					S10_sum2 += (n_count + 2);

					delta_S10 = { S10_sum1 * S10_sum2 * cn_S10 * (4 + n_count),-(2 * n_count + 9) / belta * S10_sum1 * S10_sum2 * cn_S10 * (4 + n_count)
						, S10_sum1 * S10_sum2 * (n_count - 1) / miu * cn_S10 * (4 + n_count) };

					S10 += delta_S10.at(0);
					S10_belta += delta_S10.at(1);
					S10_miu += delta_S10.at(2);
				}

				double c_miu1 = 4 * a1 * (1 - miu_pow) / 9.0;
				double c_miu2 = 4 * a2 * pow(1 - miu_pow, 2) / 81.0;

				Green_bar = -4 * S2 - c_miu1 * S6 - c_miu2 * S10;
				Green_belta = -4 * S2_belta - c_miu1 * S6_belta - c_miu2 * S10_belta;
				Green_miu = -4 * S2_miu + 8.0 / 9.0 * a1 * miu * S6 - c_miu1 * S6_miu + 16.0 / 81.0 * a2 * miu * (1 - miu_pow) * S10 - c_miu2 * S10_miu;
			}
		}
	}

	return { Green_bar,Green_belta,Green_miu };
}

