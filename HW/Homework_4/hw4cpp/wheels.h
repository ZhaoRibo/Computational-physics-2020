#ifndef WHEELS_H_
#define WHEELS_H_
#include <iostream>

using namespace std;

class vec
{
public:
	int dim;		// 维数
	double* data;	// 内容

	// 构造函数与析构函数
	vec();
	vec(int _dim);
	vec(const vec& a);
	vec(double* _data, int _dim);
	~vec();

	// 索引
	double& operator[](int i);

	// 赋值
	void operator=(const vec& a);
	
	// 一元算符，+-
	vec operator+();				// +
	vec operator-();				// -

	// 二元算符，+,-,dot,*
	const vec operator+(const vec& a);	// +
	const vec operator-(const vec& a);	// -
	friend double operator*(vec& a, vec& b);	// dot
	friend vec operator*(double& a, vec& b);
	friend vec operator*(vec& a, double& b);

	// 取模
	double mode();

	// <<
	friend ostream& operator<<(ostream& output, vec& a);

};

#endif
