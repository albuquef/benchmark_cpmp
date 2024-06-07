//
// Copyright (c) 2003 Igor Vasil'ev (igor@diima.unisa.it)
// All rights reserved.
//
// PARAMFILE.h: interface for the PARAMFILE class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_PARAMFILE_H__D964DFE9_6C6F_445B_9815_1F8089FAA08E__INCLUDED_)
#define AFX_PARAMFILE_H__D964DFE9_6C6F_445B_9815_1F8089FAA08E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "logerr.h"

// this structure is used internaly
struct PARAM
{
	String name;
	String value;
	bool used;
};


// class to getting parrametrs from file
class PARAMFILE  
{
public:
	// constructor opens a parameter file with name filename
	// print the list of existing papameters
	// lines started with # are comments
	// other lines
	// parname parvalue
	PARAMFILE(String filename, bool verbose = true);

	// destructor
	virtual ~PARAMFILE();

	// returns int value of papameter with name s
	// prodeces an error when value is out of [lb,ub]
	int intParam(String s, int lb=0x80000000, int ub=0x7fffffff);

	// returns int value of papameter with name s
	// if the parametr is missing, returns def
	// prodeces an error when value is out of [lb,ub]
	int intParam(int def, String s, int lb=0x80000000, int ub=0x7fffffff);

	// returns String value of papameter with name s in
	void strParam(String s, String dest);
	// returns String value of papameter with name s in
	void strParam(const String def, String s, String dest);

	

	// returns double value of papameter with name s
	// prodeces an error when value is out of [lb,ub]
	double doubleParam(String s, double lb=-1e30, double ub=1e30);

	// returns double value of papameter with name s
	// if the parametr is missing, returns def
	// prodeces an error when value is out of [lb,ub]
	double doubleParam(double def, String s, double lb=-1e30, double ub=1e30);
	

private:
	// return the index of parameter s, -1 if such does not exist
	int findParam(String s);
	// deletes redundant spaces in string
	int delSpace(String , const String);
	// tables of params
	PARAM* params;
	// file name
	String parfile;
	// total number of paprams
	int nPar;
};

#endif // !defined(AFX_PARAMFILE_H__D964DFE9_6C6F_445B_9815_1F8089FAA08E__INCLUDED_)
