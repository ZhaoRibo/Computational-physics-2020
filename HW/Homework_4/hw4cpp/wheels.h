#ifndef WHEELS_H_
#define WHEELS_H_
#include <iostream>

using namespace std;

class vec
{
public:
	int dim;		// ά��
	double* data;	// ����

	// ���캯������������
	vec();
	vec(int _dim);
	vec(const vec& a);
	vec(double* _data, int _dim);
	~vec();

	// ����
	double& operator[](int i);

	// ��ֵ
	void operator=(const vec& a);
	
	// һԪ�����+-
	vec operator+();				// +
	vec operator-();				// -

	// ��Ԫ�����+,-,dot,*
	const vec operator+(const vec& a);	// +
	const vec operator-(const vec& a);	// -
	friend double operator*(vec& a, vec& b);	// dot
	friend vec operator*(double& a, vec& b);
	friend vec operator*(vec& a, double& b);

	// ȡģ
	double mode();

	// <<
	friend ostream& operator<<(ostream& output, vec& a);

};

#endif
