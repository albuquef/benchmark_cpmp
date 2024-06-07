//
// Copyright (c) 2003 Igor Vasil'ev (igor@diima.unisa.it)
// All rights reserved.
//
// PARAMFILE.cpp: implementation of the PARAMFILE class.
//
//////////////////////////////////////////////////////////////////////

#include "PARAMFILE.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
#ifdef _MSC_VER
using namespace std;
#endif

PARAMFILE::PARAMFILE(String filename, bool verbose)
{
	nPar=0;
#ifdef _MSC_VER
	ifstream* from=new ifstream(filename);
#else
	ifstream* from=new ifstream(filename);
#endif
	strcpy(parfile,filename);
	if( !from || !from->is_open() )
		error("Cannot open parameters file ",filename);
	String s,s1;
	int nl=0;
	while(!from->eof())
	{
		nl++;
		from->getline(s1,200);
		int k=delSpace(s,s1);
		if((k==0)||(s[0]=='#'))
			continue;
		if(k!=2)
		{
			sprintf(s,"Error in param file %s. Line %d",filename,nl);
			error(s);
		}
		nPar++;
	}
	from->close();
	delete from;
	params=new PARAM[nPar];
	if(!params)
	{
		error("Cannot allocate PARAMFILE::params.");
	}
	{
	for(int i=0;i<nPar;i++)
	{
		params[i].used=false;
	}
	}
	from=new ifstream(filename);
	int i=0;
	while(!from->eof())
	{
		from->getline(s1,200);
		int k=delSpace(s,s1);
		if((k==0)||(s[0]=='#'))
			continue;
		int j=0,l=0;
		while(s[j]!=' ')
		{
			params[i].name[j]=s[j];
			j++;
		}
		params[i].name[j]=0;
		j++;
		while(s[j]!=0)
		{
			params[i].value[l]=s[j];
			l++;
			j++;
		}
		params[i].value[l]=0;

		i++;
	
		
	}

	if(verbose)
	{
		print("\n______________________________________________________________________\n");
		print("Parameters table\n");
		{
		for(i=0;i<nPar;i++)
		{
			sprintf(s,"%-25s%20s\n",params[i].name,params[i].value);
			print(s);
		}
		}
		print("______________________________________________________________________\n");
	}
	
	from->close();
	delete from;


}

PARAMFILE::~PARAMFILE()
{
	int nUnused=0;
	{
	for(int i=0;i<nPar;i++)
	{
		if(!params[i].used)
			nUnused++;
	}
	}
	if(nUnused)
	{
		print("\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		print("Warning:\nThere are some unused parameters in parameter file\n",
			parfile,"\n");
		print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
	}
	delete params;

}

int PARAMFILE::delSpace(String dest, const String source)
{
	int b=0;
	while(source[b++]==' ')
	{
	}
	b--;
	int e=strlen(source);
	while(source[(e--)-1]==' ')
	{
	}

	e++;
	if(b>=e)
	{
		strcpy(dest,"");
		return 0;
	}
	int neb=e-b;
	String s1;
	strncpy(s1,source+b,neb);
	s1[neb]=0;
	unsigned int i=0;
	int j=0;
	String s2;
	while(i<strlen(s1))
	{
		s2[j]=s1[i];
		i++;
		j++;
		if(s1[i]==' ')
		{
			while(s1[i+1]==' ')
			{
				i++;
			}
		}

	}
	s2[j]=0;
	strcpy(dest,s2);
	neb=0;
	{
	for(unsigned int i=0;i<strlen(dest);i++)
	{
		if(dest[i]==' ')
			neb++;
	}
	}
	return neb+1;
}

int PARAMFILE::intParam(String s, int lb, int ub)
{

	int ind=findParam(s);
	int val;
	if(ind==-1)
		error("Cannot find parametr ",s," in papameters file.");
	try
	{
		val=atoi(params[ind].value);
		params[ind].used=true;
	}
	catch(...)
	{
		error("Parameter ",s," is not integer.");
	}
	if(val>ub||val<lb)
		error("Parameter ",s," is out of range.");
	return val;
}

int PARAMFILE::intParam(int def, String s, int lb, int ub)
{

	int ind=findParam(s);
	int val;
	if(ind==-1)
	{
		print("Warning!!\n");
		print("Cannot find parametr ",s," in papameters file.\n");
		print("Default value is taken.\n");
		return def;
	}
	try
	{
		val=atoi(params[ind].value);
		params[ind].used=true;
	}
	catch(...)
	{
		error("Parameter ",s," is not integer.");
	}
	if(val>ub||val<lb)
		error("Parameter ",s," is out of range.");
	return val;
}

int PARAMFILE::findParam(String s)
{
	int ind=-1;
	{
	for(int i=0;i<nPar;i++)
	{
		if(!strcmp(s,params[i].name))
		{
			ind=i;
			break;
		}

	}
	}
	return ind;

}
//int doubleParam(String s, double lb=-1e30, double ub=1e30);
//int doubleParam(double def, String s, double lb=-1e30, double ub=1e30);

double PARAMFILE::doubleParam(String s, double lb, double ub)
{

	int ind=findParam(s);
	double val;
	if(ind==-1)
		error("Cannot find parameter ",s," in parameter file.");
	try
	{
		val=atof(params[ind].value);
		params[ind].used=true;
	}
	catch(...)
	{
		error("Parameter ",s," is not double.");
	}
	if(val>ub||val<lb)
		error("Parameter ",s," is out of range.");
	return val;
}

double PARAMFILE::doubleParam(double def, String s, double lb, double ub)
{

	int ind=findParam(s);
	double val;
	if(ind==-1)
	{
		print("Warning!!\n");
		print("Cannot find parametr ",s," in papameters file.\n");
		print("Default value is taken.\n");
		return def;
	}
	try
	{
		val=atof(params[ind].value);
		params[ind].used=true;
	}
	catch(...)
	{
		error("Parameter ",s," is not double.");
	}
	if(val>ub||val<lb)
		error("Parameter ",s," is out of range.");
	return val;
}

void PARAMFILE::strParam(String s, String dest)
{
	int ind=findParam(s);
	if(ind==-1)
		error("Cannot find parameter ",s," in parameter file.");

	strcpy(dest,params[ind].value);
	params[ind].used=true;
}

void PARAMFILE::strParam(const String def, String s, String dest)
{
	int ind=findParam(s);
	if(ind==-1)
	{
		print("Warning!!\n");
		print("Cannot find parametr ",s," in papameters file.\n");
		print("Default value is taken.\n");
		strcpy(dest,def);
	}
	else
	{
		strcpy(dest,params[ind].value);
		params[ind].used=true;
	}
	
}
