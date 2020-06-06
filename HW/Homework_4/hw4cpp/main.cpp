#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <ctime>
#include <cmath>
#include "wheels.h"
#include <iomanip>
#include <random>
using namespace std;

const long double pi = 3.14159265358979323846264338327950288;

vec initialize_x0(int n)
{
	// 初始化电子位置，cos(theta)，phi各自随机
	default_random_engine random(time(0));
	uniform_real_distribution<double> dis1(-1.0, 1.0);
	uniform_real_distribution<double> dis2(0., 2*pi);
	
	vec x0(2 * n);
	for (int i = 0; i < n; i++) {
		x0[i] = acos(dis1(random));
		x0[i + n] = dis2(random);
	}

	return x0;
}

double rij(double& thetai, double& thetaj, double& phii, double& phij)
{
	// i,j两点间距（球坐标）
	double r_ij = 2 * ((1 - cos(thetai - thetaj)) + sin(thetai) * sin(thetaj) * (1 - cos(phii - phij)));
	r_ij = sqrt(r_ij);

	return r_ij;
}

double u(vec& x)
{
	// 势能函数
	long double uu = 0;
	int n = x.dim / 2;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < i; j++) {
			uu += 1 / rij(x[i], x[j], x[i + n], x[j + n]);
		}
	}

	return uu;
}

vec g_u(vec& x)
{
	// u的梯度（解析式直接求值）
	int n = x.dim / 2;
	vec g(2 * n);
	for (int i = 0; i < n; i++) {
		double a = 0;
		double b = 0;
		for (int j = 0; j < n; j++) {
			if (j != i) {
				a = a - (sin(x[i] - x[j]) + cos(x[i]) * sin(x[j]) * (1. - cos(x[i + n] - x[j + n]))) / pow(rij(x[i], x[j], x[i + n], x[j + n]), 3.);
				b = b - (sin(x[i]) * sin(x[j]) * sin(x[i + n] - x[j + n])) / pow(rij(x[i], x[j], x[i + n], x[j + n]), 3.);
			}
		}
		g[i] = a;
		g[i + n] = b;
	}

	return g;
}

// 迭代结果类
class result
{
public:
	vec x;
	double u_x;
	int k;

	result() :x(), u_x(0.), k(0) {}
	result(vec& _x, double& _u, int& _k)
		:
		x(_x), u_x(_u), k(_k) {}
	result(const result& a)
	{
		this->x = a.x;
		this->u_x = a.u_x;
		this->k = a.k;
	}

	void operator=(const result& a)
	{
		this->x = a.x;
		this->u_x = a.u_x;
		this->k = a.k;
	}
};

// 黄金分割法线搜索
double kiefer_gold(double (*f)(vec&), vec& x, vec& d, double alpham = 1., double const eps = 1e-3)
{
	double tau = (sqrt(5) - 1.) / 2.;
	double a = 0.;
	double b = alpham;
	double x1 = a + (1 - tau) * (b - a);
	vec xx = x + x1 * d;
	double f1 = f(xx);
	double x2 = a + tau * (b - a);
	xx = x + x2 * d;
	double f2 = f(xx);
	while (true)
	{
		if ((b - a) < eps) {
			if (abs(b - alpham) <= eps) {
				alpham = alpham * 2;
				b = alpham;
				x1 = a + (1 - tau) * (b - a);
				xx = x + x1 * d;
				f1 = f(xx);
				x2 = a + tau * (b - a);
				xx = x + x2 * d;
				f2 = f(xx);
			}
			else {
				break;
			}
		}
		if (f1 > f2) {
			a = x1;
			x1 = x2;
			f1 = f2;
			x2 = a + tau * (b - a);
			xx = x + x2 * d;
			f2 = f(xx);
		}
		else {
			b = x2;
			x2 = x1;
			f2 = f1;
			x1 = a + (1 - tau) * (b - a);
			xx = x + x1 * d;
			f1 = f(xx);
		}
	}

	return (a + b) / 2;
}

// PR 算法共轭梯度
result PRCG(vec& x0, int const maxk = 1e4, double const tol = 1e-6)
{
	vec x(x0);
	int n(x.dim);
	int k = 0;
	vec g_0(g_u(x));
	double g0 = g_0.mode();
	vec g(g_0);
	vec d = -g;
	double mode_g = g.mode();

	while (k < maxk)
	{
		if (mode_g < tol) {
			break;
		}

		double alpha = kiefer_gold(u, x, d, (1. / g0));

		x = x + alpha * d;
		vec g_plus(g_u(x));
		double mode_gp = g_plus.mode();
		double beta_PR = (pow(mode_gp, 2.) - g_plus * g) / pow(mode_g, 2.);
		d = -g_plus + beta_PR * d;
		g = g_plus;
		mode_g = mode_gp;
		k++;
	}

	for (int i = 0; i < n / 2; i++) {
		while (true) {
			if (x[i] >= 0 && x[i] <= pi && x[i + n / 2] >= 0 && x[i + n / 2] < 2 * pi) {
				break;
			}
			if (x[i] < 0) {
				x[i] += pi;
			}
			else if (x[i] > pi) {
				x[i] -= pi;
			}
			if (x[i + n / 2] < 0) {
				x[i + n / 2] += 2 * pi;
			}
			else if (x[i + n / 2] >= 2 * pi) {
				x[i + n / 2] -= 2 * pi;
			}
		}
	}

	double uu = u(x);
	return result(x, uu, k);
}




int main()
{
	clock_t start, end;
	start = clock();
	
	for (int n = 2; n < 65; n++) {
		vec x0 = initialize_x0(n);
		result re = PRCG(x0, 500, 1e-3);
		int k = 0;
		while (true) {
			if (re.k < 500) {
				k += 1;
				if (k == 5) {
					break;
				}
			}
			x0 = initialize_x0(n);
			result rere = PRCG(x0, 500, 1e-3);
			if (rere.u_x < re.u_x) {
				re.x = rere.x;
				re.u_x = rere.u_x;
				re.k = rere.k;
			}
		}
		int kk = re.k;
		re = PRCG(re.x, 10000, 1e-6);
		re.k += kk;
		cout << setprecision(12) << "势能：" <<re.u_x;
		cout << endl << "迭代次数：" << re.k << endl;
	}
	
	end = clock();
	double endtime = (long double)(end - start) / CLOCKS_PER_SEC;
	cout << "Total time:" << endtime * 1000 << "ms" << endl;
	
	return 0;
}