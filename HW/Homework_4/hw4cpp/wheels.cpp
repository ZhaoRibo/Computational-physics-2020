#include <iostream>
#include "wheels.h"
using namespace std;

// 构造函数与析构函数
vec::vec()
{
	data = nullptr;
	dim = 0;
}
vec::vec(int _dim)
{
	data = new double[_dim] {0};
	dim = _dim;
}
vec::vec(const vec& a)
{
	dim = a.dim;
	data = new double[dim] {0};
	for (int i = 0; i < dim; i++) {
		data[i] = a.data[i];
	}
}
vec::vec(double* _data, int _dim)
{
	dim = _dim;
	data = new double[dim] {0};
	for (int i = 0; i < dim; i++) {
		data[i] = _data[i];
	}
}
vec::~vec()
{
	delete[]data;
}

// 索引
double& vec::operator[](int i)
{
	if (i >= dim) {
		cout << "索引超过最大值" << endl;
	}
	return data[i];
}

// 赋值
void vec::operator=(const vec& a)
{
	if (this->dim = a.dim) {
		for (int i = 0; i < this->dim; i++) {
			this->data[i] = a.data[i];
		}
	}
	else
	{
		delete[]this->data;
		this->dim = a.dim;
		this->data = new double[dim] {0};
		for (int i = 0; i < this->dim; i++) {
			this->data[i] = a.data[i];
		}
	}
}

// 一元算符，+-
vec vec::operator+()
{
	return vec(this->data, this->dim);
}
vec vec::operator-()
{
	vec a(this->dim);
	for (int i = 0; i < this->dim; i++) {
		a.data[i] = -this->data[i];
	}
	return a;
}

// 二元算符，+,-,dot,*
const vec vec::operator+(const vec& a)
{
	if (this->dim != a.dim) {
		cout << "维数不等不可相加。" << endl;
	}
	vec b(this->dim);
	for (int i = 0; i < this->dim; i++) {
		b.data[i] = this->data[i] + a.data[i];
	}
	return b;
}
const vec vec::operator-(const vec& a)
{
	if (this->dim != a.dim) {
		cout << "维数不等不可相减。" << endl;
	}
	vec b(this->dim);
	for (int i = 0; i < this->dim; i++) {
		b.data[i] = this->data[i] - a.data[i];
	}
	return b;
}
double operator*(vec& a, vec& b)
{
	if (a.dim != b.dim) {
		cout << "维数不等不可求内积。" << endl;
	}
	double c = 0;
	for (int i = 0; i < a.dim; i++) {
		c += a.data[i] * b.data[i];
	}
	return c;
}
vec operator*(double& a, vec& b)
{
	vec c(b.dim);
	for (int i = 0; i < b.dim; i++) {
		c.data[i] = a * b.data[i];
	}
	return c;
}
vec operator*(vec& a, double& b)
{
	vec c(a.dim);
	for (int i = 0; i < a.dim; i++) {
		c.data[i] = b * a.data[i];
	}
	return c;
}

// 取模
double vec::mode()
{
	double a = 0;
	for (int i = 0; i < this->dim; i++) {
		a += this->data[i] * this->data[i];
	}
	return pow(a, 0.5);
}

// <<
ostream& operator<<(ostream& output, vec& a)
{
	output << "Dim: " << a.dim << '\n' << "Vector: \n";
	for (int i = 0; i < a.dim; i++) {
		output << a.data[i] << '\t';
	}
	output << '\n';
	return output;
}