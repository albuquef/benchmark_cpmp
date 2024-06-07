// GraphComponent.h: interface for the GraphComponent class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_GRAPHCOMPONENT_H__C48B2675_8C0F_4DD3_B97A_AF21B63DBABE__INCLUDED_)
#define AFX_GRAPHCOMPONENT_H__C48B2675_8C0F_4DD3_B97A_AF21B63DBABE__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

//#include "ilcplex\cplex.h"

#define OPTGAP 0.0001

class ODMProblem;

class Pmed;

struct Arc1;
int compareArcCostArrayAZ( const void *arg1, const void *arg2 );
bool compareArcCostArrayAZ00( const Arc1 &arg1, const Arc1 &arg2 );
int compareArc1CostArrayAZ( const void *arg1, const void *arg2 );

struct Ro
{
	int i;
	char sol;
	char fix;
	double value;
};

struct Buffer
{
	unsigned int size;
	char* data;
};

class GraphComponent  
{
public:
	
	int solsize;
	int packsol(Buffer& buf);
	int unpacksol(Buffer& buf);
	int printg();
	int packprob(Buffer& buf);
	int LH(int p, int high = 0);
	int funLH(int p, double* lambda, double* grad, Ro* ro, 
		double& lb, double& ub, int* sol);
	int greedy(int lh=0, int high = 0);
	int optimize();
	int reopt(int p_);
	int optBB(int p_);
	int arcsorting();
	//int nknapvar;
	double* values;
	double* lb;
	int** solutions;
	//CPXENVptr env;
	//CPXLPptr lp;
	int createLP();
	int optCplex();
	int p1;
	int getp1();
	int p0;
	Arc1** earcs;
	int* nearcs;
	int narcs;
	int* models;
	int cind;
	ODMProblem* prob;
	int nmodels;
	GraphComponent(ODMProblem* prob_, int cind);
	GraphComponent(Buffer& buf);
	GraphComponent(ODMProblem* prob_, int cind_, Arc1 **earcs, int *nearcs);
	virtual ~GraphComponent();
	Pmed* pmed;
};

#endif // !defined(AFX_GRAPHCOMPONENT_H__C48B2675_8C0F_4DD3_B97A_AF21B63DBABE__INCLUDED_)
