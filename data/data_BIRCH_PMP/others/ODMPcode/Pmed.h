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
#include "logerr.h"




#define THETA_MIN 0.01
#define EPSILON   0.00001
#define CORENODES 4.
#define OPTGAP 0.0001
#define SATRY 10
#define INFTY 1.e30

struct Arc
{
	int i;
	double d;
};

struct Arcr
{
	int i;
	double d;
	int fix;
	Arcr():i(0),d(0),fix(0){};
};

struct Arc1;

class Pmed
{
public:
	//int _nc_;// nuber of threads
	Pmed(void);
	virtual ~Pmed(void);
	void globalfix();

	int p0, p1;

	void ver();

	void prim(int *ss);
	void SA();

	double DUALGAP;
	void fix(double lb_);
	void lang_init();
	void lang_final();
	void lang();
	void langBB();
	void lang_set1();
	void lang_set1c();
	void lang_set2();
	void lang_set3();

	void writepp(char * fn);
	void writefix(char * fn);

	void readprob(const char *fn, bool firstnumber = true);
	void readprob(int b_, const char *fn);
	void readprob(int b_, double **dd);
	void readprob(int m_, int *narcs_, Arc1 **arcs_);
	void makedist();
	void writedist();
	void readdist();
	void lagrangean();
	void dualvalue();
	void dualvalueBB();
	void subgradient();
	void step();
	void step1();
	void step2();
	void primal();
	void primalBB();
	void makelarcs();
	void writelp();
	void writelp(char *fn);
	void writelp1(char *fn);
	void writelp2(char *fn);
	int solveCPX(int p_, int *sol_);
	void writecore4();
	void savesol(int pp);
	int optimize(int p_, double ub_, int *sol_);

	void writearcs();

	double ttime;

	void savesol1();
	
	void addarcs(int j);
	void shift(Arc* aa, int nl);
	
	FILE *fd;

	double Eta;
	int fiter;
	int cfiter;
	double stopphi;
	double EPSCORE;
	int MAXCORE;
	int PRIMALFREQ;
	int FIXFREQ;
	double MAXNODE;
	double CORETIME;
	int MINARCS;
	long long MEMSIZE;
	double constf;


	inline double getdist(int i, int j)
	{
		double d = 0;
		int l;
		for(l=0; l<n; l++)
		{
			d += (coor[i][l] - coor[j][l]) * (coor[i][l] - coor[j][l]);
		}
		d = (sqrt(d));
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

	Arc **arcs1;
	int *narcs1;

	int *nlarcs;
	Arc **larcs;

	int *lbi;
	int *ubi;
	double *lm;
	double *lmub;
	double *sg;
	double phi, phiFactor, theta, lb, blb;
	Arcr *rho;
	Arcr *brho;
	double *blm;
	int *flag;
	double normSG;

	int *fixed;
	int nfix1;
	int nfix0;
	int nnfix1;
	int nnfix0;
	int nfixA1;

	bool isfull;

	char probname[200];

	double ** rhop;

	int *nfix1a;
	int *nfix0a;
	Arcr **brhoa;
	double **blma;
	int **lbia;
	int **ubia;
	double *buba;
	double *blba;


	void par1();
	void par2();
	void par3();
	void par4();

	int reopt(int p_, int *sol_);
	int optBB(int p_, int *sol_);

	void bb(double ll, int ii, int pp, int *ss);
	void bbLANG(int ii, int pp, int *ss);

private:
	int Acceptance(double, int, int, int*, Arc**, int*);
	int Acceptance1(double, int, int*);
//		void RhoFixing();
//		void Xfixing();
	boost::random::mt19937 rng;
	double Zopt, ZK;
	int *ycore, *ypot, *ysol;
	int  it_count;

};
