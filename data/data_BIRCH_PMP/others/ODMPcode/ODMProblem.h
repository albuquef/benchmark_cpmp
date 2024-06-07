// ODMProblem.h: interface for the ODMProblem class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ODMPROBLEM_H__8FCBD7BF_88E2_4D41_9B79_D772168753F6__INCLUDED_)
#define AFX_ODMPROBLEM_H__8FCBD7BF_88E2_4D41_9B79_D772168753F6__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#define MANYARCS 200000

#include <omp.h>

#include "logerr.h"
#include "CPUTIMER.h"
#include "GraphComponent.h"
#include "Convert.h"


//#include "ilcplex\cplex.h"


class Pmed;


class ODMProblem  
{
public:
	void optBB();
	void reopt();
	int *nearcs;
	Arc1** earcs;
	int optimizeParInit(int lh, int level);
	int optimizeParFinal(int lh, int level);
	int optimize(int lh, int level);
	int optimize();
	int LH(bool high = false);
	int greedy();
	int arcsorting();
	int greedyFull();
	int createKnapsackLP(bool llb=false);
	int nknapvars;
	bool solved;
	int writegraph();
	int writegraph1();
	CPUTIMER* soltimer;
	int writesol(char* name, int classe);
	int writesol1(char* name);
	//CPXENVptr env;
	//CPXLPptr lp;
	int createLP();
	int optCplexFull();
	int optCplex();
	int toobig;
	int narcs;
	int printst();
	GraphComponent** subgraphs;
	double constructioncost;
	int ngraphcomponents;
	int divide();
	ODMProblem(char* probname_, int nmodels_, __int64* modeloptions, double* modelcosts_, 
			   double* modeldemands_, char** modelnames_, int noptions_, 
			   int* hardoptions_, 
			   int dp_, int** solutions_, double* values_, double constructioncost_);

	ODMProblem(char* probname_, int p1_);

	virtual ~ODMProblem();
	char* probname;
	int nmodels;
	int dp;
	char** modelnames;
	__int64* modeloptions;
	__int64* modelhardoptions;
	int* modelcomponents;
	int noptions;
	int nhardoptions;
	int* hardoptions;
	double* modelcosts;
	double* modeldemands;
	int** solutions;
	double* values;
	double* lb;
	CPUTIMER* timer;
	int p0, p1;
};

#endif // !defined(AFX_ODMPROBLEM_H__8FCBD7BF_88E2_4D41_9B79_D772168753F6__INCLUDED_)
