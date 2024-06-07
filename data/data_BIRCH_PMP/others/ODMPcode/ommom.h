// ommom.h: interface for the Commom class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_OMMOM_H__DD17074F_EC1C_4CBC_A4DD_E81658BE8659__INCLUDED_)
#define AFX_OMMOM_H__DD17074F_EC1C_4CBC_A4DD_E81658BE8659__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "logerr.h"

class Commom  
{
public:
	int funLH(int p, double lambda, double& grad, double& lb, double& ub, int* sol);
	int LH(int p);
	int optimize();
	int printdata();
	int n, p0, p1;
	int* tmp;
	double** costs;
	int** ps;
	double* values;
	double* lb;
	int** solutions;
	Commom(int n_, double** costs_, int** ps_, int p0_, int p1_);
	virtual ~Commom();
	int DP();
	int DP1();
};

#endif // !defined(AFX_OMMOM_H__DD17074F_EC1C_4CBC_A4DD_E81658BE8659__INCLUDED_)
