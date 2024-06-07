// Convert.h: interface for the Convert class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CONVERT_H__C8D7D68B_3D1A_4949_BC14_592FA1A5C73E__INCLUDED_)
#define AFX_CONVERT_H__C8D7D68B_3D1A_4949_BC14_592FA1A5C73E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "logerr.h"
#include <string.h>
#include "CPUTIMER.h"


#define MAXOPTIONS  1000
#define MAXNAMELEN  10
#define MAXSTDS     100
#define MAXMODELS   100000

struct Arc1
{
	int u;//, v;
	double cost;
};

struct Arc2
{
	int u, v;
	double cost;
};

#define __int64 unsigned int

class Convert  
{
public:
	double constructioncost;
	String fn1, fn2, fn3, fn4, fn5;
	int writelp();

	int* nearcs;
	
	int writegraph();
	
	int p0;
	int narcs;
	int printprob();

	
	int genergraph();
	char optionnames[MAXOPTIONS][MAXNAMELEN];
	char** realmodelnames;
	int noptions;

	double optioncosts[MAXOPTIONS];
	int hardoptions[MAXOPTIONS];

	int nstds;
	char sdtnames[MAXSTDS][MAXNAMELEN];

	int nmodels;
	String modelnames[MAXMODELS];
	double modelcosts[MAXMODELS];
	int modelstds[MAXMODELS];
	double modeldemands[MAXMODELS];
	__int64 modelmasks[MAXMODELS];
	__int64 modelhardmasks[MAXMODELS];


	int writeprob();
	int getmodels();
	int getcost();
	int getnames();
	int process();
	Convert(char* probname_, char* fn1_, char* fn2, char* fn3, char* fn4_, char* fn5_);
	virtual ~Convert();

	String probname;

};

#endif // !defined(AFX_CONVERT_H__C8D7D68B_3D1A_4949_BC14_592FA1A5C73E__INCLUDED_)
