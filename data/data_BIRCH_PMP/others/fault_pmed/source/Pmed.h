/********************************************************************
This is complementary software for the paper 
P. Avella, M. Boccia, S. Salerno, I. Vasilyev, 
An aggregation heuristic for large scale p-median problem, 
Computers and Operations Research. 2010, doi:10.1016/j.cor.2011.09.016.


(C) 2008-2011, by Igor Vasilyev (vasilievigor at gmail.com). 
It is free of charge for any way of use, but the reference 
to the author is encouraged.
IT IS PROVIDED WITHOUT ANY WARRANTY,
but comments and bugs can be sent to vasilievigor at gmail.com.
********************************************************************/


#pragma once
#include "main.h"


#define CPX
//#define XPRS

#define DUALGAP 0.00001
#define THETA_MIN 0.01
#define EPSILON   0.00001
#define CORENODES 4.

//#define INTDIST

struct Arc
{
	double d;
	int i;
};

class Pmed
{
public:
	Pmed(void);
	virtual ~Pmed(void);

	void ver();


	void lang_init();
	void lang_final();
	void lang();
	void lang_set1();
	void lang_set1c();
	void lang_set2();
	void lang_set3();

	void readprob(const char *fn, bool firstnumber = true);
	void readprob(int b_, const char *fn);
	void readprob(int b_, double **dd);
	void readprobmatr(char *fn);

	void makedist();
	void writedist();
	void readdist();
	void lagrangean();
	//void lagrangean1();
	void dualvalue();
	void subgradient();
	void step();
	void step1();
	void step2();
	void primal();
	void makelarcs();
	void writelp();
	/*void writecore();
	void writecore1();
	void writecore2();*/
	void writecore3();
	void writecore4();
	void savesol(int pp);

	void savesol1();
	
	void addarcs(int j);

	void greedy();
	void greedy1();

	FILE *fd;

	double Eta;
	int fiter;
	int cfiter;
	double stopphi;
	double EPSCORE;
	int MAXCORE;
	int PRIMALFREQ;
	double MAXNODE;
	double CORETIME;
	int MINARCS;
	long long MEMSIZE;
  int r;  // Number of medians to connect with each node


	inline double getdist(int i, int j)
	{
		double d = 0;
		int l;
		for(l=0; l<n; l++)
		{
			d += (coor[i][l] - coor[j][l]) * (coor[i][l] - coor[j][l]);
		}
		d = (sqrt(d));
#ifdef INTDIST
		d = floor( d );
#endif
		return d;
	}

	int dispInterval;

	double **coor;
	int n;
	int m;
	int p;
	int m1, m2;
	double bub;
	int *bsol;
	Arc **arcs;
	double *c;
	int *narcs;
	int numarcs;

	int *nlarcs;
	Arc **larcs;

	int *lbi;
	int *ubi;
	double *lm;
	double *lmub;
	double *sg;
	double phi, phiFactor, theta, lb, blb;
	Arc *rho;
	Arc *brho;
	double *blm;
	int *flag;
	double normSG;

	bool isfull;

	char probname[200];

};
