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


#include "Pmed.h"
#include <cplex.h>

extern int _nc_;

inline void swap(Arcr &a1, Arcr &a2)
{
	Arcr tmp;
	tmp = a1;
	a1 = a2;
	a2 = tmp;
}


void checkselection(Arcr* a, int m, int p)
{
	int i;
	for(i=0; i<p; i++)
	{
		if(a[p].d < a[i].d)
			error("E1");
	}
	for(i=p+1; i<m; i++)
	{
		if(a[p].d > a[i].d)
			error("E2");
	}
}

int partition(Arcr* a, int left, int right, int pivoti)
{
	double pivot = a[pivoti].d;
	swap(a[pivoti], a[right]);
	int storei = left;
	int i;
	for(i=left; i<right; i++)
	{
		if(a[i].d < pivot)
		{
			swap(a[storei], a[i]);
			storei++;
		}
	}
	swap(a[right], a[storei]);
	return storei;
}

void select(Arcr* a, int left, int right, int k);

void selectp(Arcr* a, int left, int right, int k)
{
//#pragma omp parallel
	{
//	#pragma omp single nowait
		{
			select(a, left, right, k);
		}
	}
}

void select(Arcr* a, int left, int right, int k)
{
	if(right == left)
		return;
	int pivoti = (right + left) / 2;
	pivoti = partition(a, left, right, pivoti);
	if(k == pivoti)
		return;

	if(right - left < 500)
	{
		if(k < pivoti)
			select(a, left, pivoti - 1, k);
		else
			select(a, pivoti + 1, right, k);
		return;
	}

	if(k < pivoti)
	{
//#pragma omp task 
		{
			select(a, left, pivoti - 1, k);
		}
	}
	else
	{
//#pragma omp task 
		{
			select(a, pivoti + 1, right, k);
		}
	}
	return;
}


struct Pair
{
	int i,j;
};

int comparePairAZ(const void* item1, const void* item2)
{
	Pair *it1 = (Pair*) item1;
	Pair *it2 = (Pair*) item2;
	if(it1->i > it2->i)
		return 1;
	if(it1->i < it2->i)
		return -1;
	if(it1->j > it2->j)
		return 1;
	return -1;
}

int compareDoubleAZ(const void* item1, const void* item2)
{
	double *it1 = (double*) item1;
	double *it2 = (double*) item2;
	if(*it1 > *it2)
		return 1;
	if(*it1 < *it2)
		return -1;
	return 0;
}

int compareArcAZ(const void* item1, const void* item2)
{
	Arc *it1 = (Arc*) item1;
	Arc *it2 = (Arc*) item2;
	if(it1->d > it2->d)
		return 1;
	if(it1->d < it2->d)
		return -1;
	return 0;
}


int compareArcrAZ1(const void* item1, const void* item2)
{
	Arcr *it1 = (Arcr*) item1;
	Arcr *it2 = (Arcr*) item2;
	if(it1->d > it2->d)
		return 1;
	if(it1->d < it2->d)
		return -1;
	return 0;
}

int compareArcrAZ2(const void* item1, const void* item2)
{
	Arcr *it1 = (Arcr*) item1;
	Arcr *it2 = (Arcr*) item2;
	if(it1->fix == it2->fix)
	{
		if(it1->d > it2->d)
			return 1;
		if(it1->d < it2->d)
			return -1;
	}
	if(it1->fix == 1)
		return -1;
	if(it1->fix == 0)
		return 1;
	if(it2->fix == 1)
		return 1;
	if(it2->fix == 0)
		return -1;
	return 0;
}

int compareArcrAZ(const void* item1, const void* item2)
{
	Arcr *it1 = (Arcr*) item1;
	Arcr *it2 = (Arcr*) item2;
	if(it1->fix == 1 || it2->fix == 0)
		return -1;	
	if(it1->fix == 0 || it2->fix == 1)
		return 1;
	if(it1->d > it2->d)
		return 1;
	if(it1->d < it2->d)
		return -1;
	return 0;
}

int compareArcZA(const void* item1, const void* item2)
{
	Arc *it1 = (Arc*) item1;
	Arc *it2 = (Arc*) item2;
	if(it1->d > it2->d)
		return -1;
	if(it1->d < it2->d)
		return 1;
	return 0;
}


void Pmed::ver()
{
	char fn[200];
	//sprintf(fn, "%s.sol", probname);
	sprintf(fn, "%s-%d.sol", probname, p);

	FILE *ff = fopen(fn, "r");

	if(ff == NULL)
		error("Cannot open ", fn);

	int p1 = 0;
	double ub1 = 0.;
	fscanf(ff, "p = %ld\n", &p1);
	fscanf(ff, "bub = %lf\n", &ub1);
	int ch = 0;
	fscanf(ff, "Check = %ld\nMedians:\n", &ch);

	int *meds = new int[p1];
	int i;
	for(i=0; i<p1; i++)
	{
		fscanf(ff, "%ld", &meds[i]);
		//printf("%d\n", meds[i]);
	}
	fclose(ff);
	double ub2 = 0., min, dd;
	int mini;
	for(i=0; i<m; i++)
	{
		min = 1.e30;
		mini = 0;
		for(int j=0; j<p1; j++)
		{
			dd = getdist(i, meds[j]);
			if(min > dd)
			{
				min = dd;
				mini = meds[j];
			}
		}
		ub2 += min;
	}
	bool bb = (p == p1) && (ch) && (fabs(ub1-ub2)<1.e-4);
	printf("p=%d\n%lf\n%lf\ncheck %d\n", p, ub1, ub2, bb);
	ff = fopen("ver.txt", "a");
	fprintf(ff, "%s\t%s\t%d\n", probname, fn, bb);
	fclose(ff);
	delete[] meds;
}

void Pmed::makelarcs()
{
	int i, j;
	if(larcs)
	{
		for(i=0; i<m; i++)
			delete[] larcs[i];
		delete[] larcs;
	}
	if(nlarcs)
		delete[] nlarcs;
	
	nlarcs = new int[m];

	for(i=0; i<m; i++)
	{
		nlarcs[i] = 0;
	}
	for(i=0; i<m; i++)
	{
		for(j=0; j<narcs[i]; j++)
		{
			nlarcs[arcs[i][j].i] ++;
		}
	}
	larcs = new Arc*[m];
	for(i=0; i<m; i++)
	{
		larcs[i] = new Arc[nlarcs[i]];
		nlarcs[i] = 0;
	}
	for(i=0; i<m; i++)
	{
		for(j=0; j<narcs[i]; j++)
		{
			larcs[arcs[i][j].i][nlarcs[arcs[i][j].i]].i = i;
			larcs[arcs[i][j].i][nlarcs[arcs[i][j].i]].d = arcs[i][j].d;
			nlarcs[arcs[i][j].i] ++;
		}
	}
	if(1)
	{
		int na = 0;
		for(i=0; i<m; i++)
		{
			na += nlarcs[i];
		}
		printf("na = %d\n", na);

		na = 0;
		for(i=0; i<m; i++)
		{
			for(j=0; j<nlarcs[i]; j++)
			{
				if( larcs[i][j].d - lm[larcs[i][j].i] < 0. )
					na ++;
			}
		}
		printf("na = %d\n", na);
	}
	
	if(0)
	{
		printf("Entering arcs:\n");
		for(i=0; i<m; i++)
		{
			printf("ff = %d\n", i);
			for(j=0; j<narcs[i]; j++)
			{
				printf("%d %.2f\n", arcs[i][j].i, arcs[i][j].d);
			}
		}
		printf("Leaving arcs:\n");
		for(i=0; i<m; i++)
		{
			printf("ff = %d\n", i);
			for(j=0; j<nlarcs[i]; j++)
			{
				printf("%d %.2f\n", larcs[i][j].i, larcs[i][j].d);
			}
		}
	}

	//error("vv");
}

void Pmed::savesol1()
{
	printf("Saving solution...\n");
	printf("objval = %.4lf\n", bub);
	/*double df = 0.;
	int i, k = 0;
	for(i=0; i<m; i++)
	{
		if(bsol[i] == i)
		{
			k++;
			continue;
		}
		df += getdist(i, bsol[i]);
	}
	printf("objval = %.4lf\n", df);

	bool bb = true;
	if(fabs(df - bub) > 0.001)
		bb = false;
	int *medians = new int[p];
	int nn = 0;
	for(i=0; i<m; i++)
	{
		if(bsol[i] == i)
		{
			medians[nn] = i;
			nn ++;
		}
	}*/

	char fn[200];
	sprintf(fn, "%s-%d.sol", probname, p);
	FILE *ff = fopen(fn, "w");

	fprintf(ff, "p = %d\n", p);
	fprintf(ff, "bub = %lf\n", bub);
	//fprintf(ff, "Check = %d\n", bb);
	fprintf(ff, "Medians:\n");
	int i;
	for(i=0; i<p; i++)
		fprintf(ff, "%d\n", bsol[i]);
	fclose(ff);

	//delete[] medians;
	
}

void Pmed::savesol(int pp)
{
	printf("Cheking solution...\n");
	printf("objval = %.4lf\n", bub);
	double df = 0.;
	int i, k = 0;
	for(i=0; i<m; i++)
	{
		if(bsol[i] == i)
		{
			k++;
			continue;
		}
		df += getdist(i, bsol[i]);
	}
	printf("objval = %.4lf\n", df);
	if(0)
	{
		printf("k = %d\n", k);
		for(i=0; i<m; i++)
		{
			if(bsol[i] != i)
				continue;
			printf("y%d = 1\n", i);
		}
	}

	int *medians = new int[p];
	
	int *mbeg = new int[p+1];
	Pair *sol = new Pair[m];
	for(i=0; i<p+1; i++)
	{
		mbeg[i] = 0;
	}
	int j = 0;
	for(i=0; i<m; i++)
	{
		sol[i].i = bsol[i];
		sol[i].j = i;
	}
	qsort(sol, m, sizeof(Pair), comparePairAZ);

	int ki = 0;
	for(i=0; i<m; i++)
	{
		if(sol[i].i != sol[i].j)
			continue;
		medians[ki] = sol[i].i;
		ki++;
	}

	j=0;
	mbeg[0] = 0;
	for(i=0; i<m; i++)
	{
		if(sol[mbeg[j]].i != sol[i].i)
		{
			j++;
			mbeg[j] = i;
		}
		if(j == p-1)
		{
			j++;
			mbeg[j] = m;
			break;
		}
	}
	//mbeg[p] = m;

	double ** dd = new double*[p];
	for(i=0; i<p; i++)
	{
		dd[i] = new double[p];
	}
	for(i=0; i<p; i++)
	{
		for(j=0; j<p; j++)
		{
			dd[i][j] = 0;
			for(int k=mbeg[j]; k<mbeg[j+1]; k++)
			{
				dd[i][j] += getdist(sol[mbeg[i]].i, sol[k].j);
			}
			//dd[j][i] = dd[i][j];
		}
	}
	


	/*char fn[200];
	sprintf(fn, "%s-%d.dis", probname, p);
	FILE *ff = fopen(fn, "wb");
	
	fwrite(&p, sizeof(int), 1, ff);
	for(i=0; i<p; i++)
	{
		fwrite(dd[i], sizeof(double), p, ff);
	}

	fclose(ff);*/

	delete[] mbeg;
	delete[] sol;

	int oldp = p;

	if(1)
	{
		Pmed pmed;
		pmed.readprob(p, dd);
		pmed.bub = 1.e30;
		pmed.CORETIME = 60;
		pmed.MAXNODE = 3;
		pmed.MAXCORE = 4 * pmed.m;
		pmed.p = pp;
		pmed.MINARCS = 3;
		pmed.lang_set1c();
		pmed.PRIMALFREQ = 10000000;
		
		pmed.lang_init();
		pmed.lang();
		pmed.writecore4();
		
		pmed.lang_set2();
		pmed.lang();
		pmed.writecore4();
		
		pmed.lang_set3();
		pmed.lang();
		pmed.lang_final();

		bub = pmed.bub;

		printf("ub = %lf\n", bub);
		for(i=0; i<m; i++)
		{
			flag[i] = 0;
		}
		for(i=0; i<p; i++)
		{
			if(pmed.bsol[i] != i)
				continue;
			flag[medians[i]] = 1;
		}
		p = pp;
		primal();
		printf("ub = %lf\n", bub);

	}
	delete[] medians;
	for(i=0; i<oldp; i++)
		delete[] dd[i];
	delete[] dd;
}


void Pmed::writecore4()
{

	printf("SA..\n");
	double tt = CPU_time();

	int nnodes = m;
	int nnarcs = 0;
	int i, j;
	for(i=0; i<m; i++)
		nnarcs += narcs[i];

	printf("Core nodes = %d\nCore leaving arcs = %d\n", nnodes, MAXCORE);

	int *ncore = narcs;
	Arc ** core = arcs;

	///////Simulated annealing!!!///////////////////////////////////////////////////////////////

	int *flag2 = new int[m];
	ycore = new int[p];

	int pot = nnodes - p;
	ypot = new int[pot];
	ysol = new int[p];

	for(i=0;i<m;i++)
		flag2[i] = 0;
	for(i=0;i<p;i++)
	{
		flag2[brho[i].i] = 1;
		ycore[i] = brho[i].i;
	}
	int k = 0;
	for(i=p;i<nnodes;i++)
	{
		ypot[k] = brho[i].i;
		k++;
	}

	ZK = 0.;
	int h;
	int nas=0;
	for(k=0; k<m; k++)
	{
		if(flag2[k] == 1)
			continue;
		for(h=0; h<ncore[k]; h++)
		{
			i = core[k][h].i;
			if(flag2[i]==1)
			{
				ZK += core[k][h].d;
				nas++;
				break;
			}
		}
	}
	
	Zopt = ZK;

	printf("SA begins. Upper bound Zopt = %lf\n", Zopt);
	//rng.seed(clock()+time(NULL));//+_getpid());
	rng.seed(0);

	double t = 0.05;//0.005;//effective 0.005;
	double tmin = 0.001;//0.0005;//effective 0.0005;
	double q; 
	int M_in = 500;//nnodes/2; // it was 2*nnodes
	int M_out = 100;//nnodes/4; // nnnodes
	int M_fr = 40*nnodes; //40*nnodes
	int swap_flag;
	double t0 = t;
	it_count = 0;
	q = pow(tmin/t0, (double)1/(M_out-1));

	printf("Iterations number %d\n", M_in*M_out);
	k = 0;
	for(i = 0; i <= M_in*M_out; i++)
	{
		it_count++;
		swap_flag = Acceptance(t, pot, nnodes, flag2, core, ncore);
		if(swap_flag == 1)
		{
			printf("");
		}
		if(i&&(i % M_in == 0))
		{
			t = q * t;
		}
		if(it_count == M_fr)
		{
			break;
		}
		k++;
	}
	

	printf("SA ends %d iterations. Upper bound Zopt = %lf\n", it_count, Zopt);

	if(Zopt < bub)
	{
		bub = Zopt;
		for(i=0; i<m; i++)
		{
			flag2[i] = 0;
		}
		for(i=0; i<p; i++)
		{
			flag2[ycore[i]] = 1;
		}
		//memcpy(bsol, ysol, p*sizeof(int));
		double dff = 0.;
		for(i=0; i<m; i++)
		{
			if(flag2[i])
			{
				bsol[i] = i;
				continue;
			}
			for(j=0; j<narcs[i]; j++)
			{
				int k = arcs[i][j].i;
				if(flag2[k])
				{
					bsol[i] = k;
					dff += arcs[i][j].d;
				}
			}
		}
		if(fabs(bub - dff) > EPSILON)
			error("bub != dff");
	}

	//double GAP = ((bub - blb) / bub * 100.0);

	delete[] ycore;
	delete[] ypot;
	delete[] ysol;

	
	printf("Done in %.2f seconds\n", CPU_time() - tt);

	printf("Solving core (wait %.2f seconds)\n", CORETIME);

	delete[] flag2;
}

void Pmed::SA()
{

	print("SA:");
	double tt = CPU_time();
	int i, j;

	///////Simulated annealing!!!///////////////////////////////////////////////////////////////

	int *flag2 = new int[m];
	ycore = new int[p];

	int pot = m - nfix0 - p;
	ypot = new int[pot];
	ysol = new int[p - nfix1];

	for(i=0;i<m;i++)
		flag2[i] = 0;
	for(i=0;i<p;i++)
	{
		flag2[brho[i].i] = 1;
	}
	for(i=nfix1;i<p;i++)
	{
		ycore[i-nfix1] = brho[i].i;
	}

	int k = 0;
	for(i=p; i<m-nfix0; i++)
	{
		ypot[k] = brho[i].i;
		k++;
	}

	ZK = 0.;
	int h;
	for(k=0; k<m; k++)
	{
		if(flag2[k] == 1)
			continue;
		for(h=0; h<narcs[k]; h++)
		{
			i = arcs[k][h].i;
			if(flag2[i]==1)
			{
				ZK += arcs[k][h].d;
				break;
			}
		}
	}
	
	Zopt = ZK;

	char ss[500];
	sprintf(ss, "%14.2f ->", Zopt);
	print(ss);
	rng.seed(clock()+time(NULL));//+_getpid());
	//rng.seed(0);

	double t = 500.;//0.005;//effective 0.005;
	double tmin = 0.001;//0.0005;//effective 0.0005;
	double q; 
	int M_in = (m - nfix0 - nfix1);//nnodes/2; // it was 2*nnodes
	int M_out = (m - nfix0 - nfix1)/2;//nnodes/4; // nnnodes
	int M_fr = m - nfix0 - nfix1; //40*nnodes
	int swap_flag;
	double t0 = t;
	it_count = 0;
	q = pow(tmin/t0, (double)1/(M_out-1));

	//printf("Iterations number %d\n", M_in*M_out);
	k = 0;
	for(i = 0; i <= M_in*M_out; i++)
	{
		it_count++;
		swap_flag = Acceptance1(t, pot, flag2);
		if(i&&(i % M_in == 0))
		{
			t = q * t;
		}
		if(it_count == M_fr)
		{
			//break;
		}
		k++;
	}
	

	//printf("SA ends %d iterations. Upper bound Zopt = %lf\n", it_count, Zopt);
	sprintf(ss, "%14.2f", Zopt);
	print(ss);
	if(Zopt < bub)
	{
		bub = Zopt;
		for(i=0; i<m; i++)
		{
			flag2[i] = 0;
		}
		for(i=0; i<nfix1; i++)
		{
			flag2[brho[i].i] = 1;
		}
		for(i=nfix1; i<p; i++)
		{
			flag2[ysol[i-nfix1]] = 1;
		}
		//memcpy(bsol, ysol, p*sizeof(int));
		double dff = 0.;
		for(i=0; i<m; i++)
		{
			if(flag2[i])
			{
				bsol[i] = i;
				continue;
			}
			for(j=0; j<narcs[i]; j++)
			{
				int k = arcs[i][j].i;
				if(flag2[k])
				{
					bsol[i] = k;
					dff += arcs[i][j].d;
					break;
				}
			}
		}
		if(fabs(bub - dff) > EPSILON)
		{
			writelp();
			printf("bub = %12.6f dff = %12.6f\n", bub, dff);
			error("bub != dff");
		}
	}

	//double GAP = ((bub - blb) / bub * 100.0);

	delete[] ycore;
	delete[] ypot;
	delete[] ysol;

	
	sprintf(ss,"%10.2f\n", CPU_time() - tt);
	print(ss);

	delete[] flag2;
}


int Pmed::Acceptance(double t, int pot, int nnodes, int *flag2, Arc **core, int *ncore)
{
	int accept;
	boost::random::uniform_int_distribution<> open_(0, pot-1);
	boost::random::uniform_int_distribution<> close_(0,p-1);

	boost::uniform_01<boost::mt19937> zeroone(rng);

	int po, pc;
	po = open_(rng);
	pc = close_(rng);

	flag2[ypot[po]] = 1;
	flag2[ycore[pc]] = 0;


	int i;
	double Zk=0;
	int nas = 0;
#pragma omp parallel for reduction(+:Zk, nas)
	for(i=0;i<m;i++)
	{
		if(flag2[i])
			continue;
		for(int j=0;j<ncore[i];j++)
		{
			if(flag2[core[i][j].i])
			{
				Zk+=core[i][j].d;
				nas++;
				break;
			}
		}
	}
	if(nas == m-p)
	{
		if((Zk - ZK)<0)
		{
			ZK = Zk;
			accept = 1;

			if(Zk<Zopt)
			{
			Zopt = Zk;
			memcpy(ysol, ycore, p*sizeof(int));
			it_count = 0;
			}
	}
	else
	{
		accept = -1;

		double pr = zeroone();
		double qq = exp(-(Zk - ZK)/ t);
		if( qq > pr)
		{
			ZK = Zk;
			accept = 1;
		}
	}
	}
	else
		accept = -1;
	if(accept > 0)
	{
		int inter = ycore[pc];
		ycore[pc] = ypot[po];
		ypot[po] = inter;
	}
	else
	{
		flag2[ypot[po]] = 0;
		flag2[ycore[pc]] = 1;
	}
	return accept;
}

int Pmed::Acceptance1(double t, int pot, int *flag2)
{
	int accept;
	boost::random::uniform_int_distribution<> open_(0, pot-1);
	boost::random::uniform_int_distribution<> close_(0,p-nfix1-1);

	boost::uniform_01<boost::mt19937> zeroone(rng);

	int po, pc;
	po = open_(rng);
	pc = close_(rng);

	flag2[ypot[po]] = 1;
	flag2[ycore[pc]] = 0;


	int i;
	double Zk=0;
	for(i=0;i<m;i++)
	{
		if(flag2[i])
			continue;
		for(int j=0;j<narcs[i];j++)
		{
			if(flag2[arcs[i][j].i])
			{
				Zk+=arcs[i][j].d;
				break;
			}
		}
	}
	
	if((Zk - ZK)<0)
	{
		ZK = Zk;
		accept = 1;
	}
	else
	{
		accept = -1;

		double pr = zeroone();
		double qq = exp(-(Zk - ZK)/ t);
		if( qq > pr)
		{
			ZK = Zk;
			accept = 1;
		}
	}
	
	if(accept > 0)
	{
		//printf("*");
		int inter = ycore[pc];
		ycore[pc] = ypot[po];
		ypot[po] = inter;
		if(Zk<Zopt)
		{
			Zopt = Zk;
			memcpy(ysol, ycore, (p-nfix1)*sizeof(int));
			it_count = 0;
		}
	}
	else
	{
		flag2[ypot[po]] = 0;
		flag2[ycore[pc]] = 1;
	}
	return accept;
}

int iij = 0;


void Pmed::writelp(char *fn)
{
	//char fn[200];
	double tt = CPU_time();
	printf("Writing LP..\n");
	//sprintf(fn, "%s-%2.2d.lp", probname, iij);
	FILE * ff = fopen(fn, "w");

	fprintf(ff, "min\n");
	int nn = 0;
	int i, j;

	for(i=0; i<m; i++)
	{
		fixed[brho[i].i] = i;
	}
	for(i=nfix1; i<m-nfix0; i++)
	{
		if(nn)
			fprintf(ff, " + ");
		fprintf(ff, " 0 y%d", brho[i].i);
		nn ++;
		if(nn%5 == 0) 
			fprintf(ff, "\n");
	}
	for(i=0; i<m; i++)
	{
		if(brho[fixed[i]].fix == 1)
			continue;
		for(j=0; j<narcs[i]; j++)
		{
			if(brho[fixed[arcs[i][j].i]].fix == 0)
				continue;
			fprintf(ff, "+ %g x%d,%d", arcs[i][j].d, arcs[i][j].i, i);
			nn ++;
			if(nn%5 == 0)
				fprintf(ff, "\n");
		}
	}
	fprintf(ff, "\ns.t.\n");
	//printf("\nff\n");
	for(i=0; i<m; i++)
	{
		if(brho[fixed[i]].fix == 1)
			continue;
		nn = 0;
		for(j=0; j<narcs[i]; j++)
		{
			if(brho[fixed[arcs[i][j].i]].fix == 0)
				continue;
			if(nn)
				fprintf(ff, " + ");
			fprintf(ff, "x%d,%d", arcs[i][j].i, i);
			nn ++;
			if(nn%5 == 0)
				fprintf(ff, "\n");
		}
		if(!narcs[i])
		{
			//printf("%d\n", i);
			continue;
		}
		if(brho[fixed[i]].fix == -1)
			fprintf(ff, "+ y%d = 1\n", i);
		else
			fprintf(ff, " = 1\n", i);
	}
	for(i=0; i<m; i++)
	{
		if(brho[fixed[i]].fix == 1)
			continue;
		for(j=0; j<narcs[i]; j++)
		{
			if(brho[fixed[arcs[i][j].i]].fix == -1)
				fprintf(ff, "y%d - x%d,%d >= 0\n", arcs[i][j].i, arcs[i][j].i, i);
		}
	}
	nn = 0;
	for(i=nfix1; i<m-nfix0; i++)
	{
		if(nn)
			fprintf(ff, " + ");
		fprintf(ff, "y%d", brho[i].i);
		nn ++;
		if(nn%5 == 0)
			fprintf(ff, "\n");
	}
	fprintf(ff, " = %d\n", p - nfix1);
	//fprintf(ff, "\nbounds\n");
	//nn = 0;
	/*for(i=0; i<nfix1; i++)
	{
		fprintf(ff, " y%d = 1\n", brho[i].i);
	}
	for(i=m-nfix0; i<m; i++)
	{
		fprintf(ff, " y%d = 0\n", brho[i].i);
	}*/
	fprintf(ff, "\nbinaries\n");
	nn = 0;
	for(i=nfix1; i<m-nfix0; i++)
	{
		fprintf(ff, " y%d", brho[i].i);
		nn ++;
		if(nn%5 == 0) 
			fprintf(ff, "\n");
	}
	fprintf(ff, "\nend\n");
	fclose(ff);

	if(0)
	{
		char fn1[500];
		sprintf(fn1, "%s-%2.2d.txt", probname, iij);
		ff = fopen(fn1, "w");
		nn = 0;
		for(i=0; i<m; i++)
			nn += narcs[i];
		fprintf(ff, "%ld %ld\n", m, nn);
		printf("Done in %.2f seconds\n", CPU_time() - tt);
		for(i=0; i<m; i++)
		{
			for(j=0; j<narcs[i]; j++)
			{
				fprintf(ff, "%d %d %g\n", arcs[i][j].i+1, i+1, arcs[i][j].d);
			}
		}
		fclose(ff);
	}
	iij++;
}

void Pmed::writelp1(char *fn)
{
	//char fn[200];
	double tt = CPU_time();
	printf("Writing LP..\n");
	//sprintf(fn, "%s-%2.2d.lp", probname, iij);
	FILE * ff = fopen(fn, "w");

	fprintf(ff, "min\n");
	int nn = 0;
	int i, j;

	for(i=0; i<m; i++)
	{
		fixed[brho[i].i] = i;
	}
	for(i=nfix1; i<m-nfix0; i++)
	{
		if(nn)
			fprintf(ff, " + ");
		fprintf(ff, " 0 y%d", brho[i].i);
		nn ++;
		if(nn%5 == 0) 
			fprintf(ff, "\n");
	}
	for(i=0; i<m - nfixA1; i++)
	{
		for(j=0; j<narcs[brho[i].i]; j++)
		{
			fprintf(ff, " + %.6f x%d,%d", arcs[brho[i].i][j].d, arcs[brho[i].i][j].i, brho[i].i);
			nn ++;
			if(nn%5 == 0)
				fprintf(ff, "\n");
		}
	}
	fprintf(ff, " + %.6f z\n", constf);
	fprintf(ff, "\ns.t.\n");
	//printf("\nff\n");
	for(i=nfix1; i<m - nfixA1; i++)
	{
		nn = 0;
		for(j=0; j<narcs[brho[i].i]; j++)
		{
			if(nn)
				fprintf(ff, " + ");
			fprintf(ff, "x%d,%d", arcs[brho[i].i][j].i, brho[i].i);
			nn ++;
			if(nn%5 == 0)
				fprintf(ff, "\n");
		}
		if(!narcs[brho[i].i])
		{
			//printf("%d\n", i);
			continue;
		}
		if(brho[i].fix == -1)
			fprintf(ff, "+ y%d = 1\n", brho[i].i);
		else
			fprintf(ff, " = 1\n");
	}
	for(i=nfix1; i<m - nfixA1; i++)
	{
		for(j=0; j<narcs[brho[i].i]; j++)
		{
			if(brho[fixed[arcs[brho[i].i][j].i]].fix == -1)
				fprintf(ff, "y%d - x%d,%d >= 0\n", arcs[brho[i].i][j].i, arcs[brho[i].i][j].i, brho[i].i);
		}
	}
	nn = 0;
	for(i=nfix1; i<m-nfix0; i++)
	{
		if(nn)
			fprintf(ff, " + ");
		fprintf(ff, "y%d", brho[i].i);
		nn ++;
		if(nn%5 == 0)
			fprintf(ff, "\n");
	}
	fprintf(ff, " = %d\n", p - nfix1);
	fprintf(ff, "\nbounds\nz = 1\n");
	//nn = 0;
	/*for(i=0; i<nfix1; i++)
	{
		fprintf(ff, " y%d = 1\n", brho[i].i);
	}
	for(i=m-nfix0; i<m; i++)
	{
		fprintf(ff, " y%d = 0\n", brho[i].i);
	}*/
	fprintf(ff, "\nbinaries\n");
	nn = 0;
	for(i=nfix1; i<m-nfix0; i++)
	{
		fprintf(ff, " y%d", brho[i].i);
		nn ++;
		if(nn%5 == 0) 
			fprintf(ff, "\n");
	}
	fprintf(ff, "\nend\n");
	fclose(ff);

	if(0)
	{
		char fn1[500];
		sprintf(fn1, "%s-%2.2d.txt", probname, iij);
		ff = fopen(fn1, "w");
		nn = 0;
		for(i=0; i<m; i++)
			nn += narcs[i];
		fprintf(ff, "%ld %ld\n", m, nn);
		printf("Done in %.2f seconds\n", CPU_time() - tt);
		for(i=0; i<m; i++)
		{
			for(j=0; j<narcs[i]; j++)
			{
				fprintf(ff, "%d %d %g\n", arcs[i][j].i+1, i+1, arcs[i][j].d);
			}
		}
		fclose(ff);
	}
	iij++;
}

void Pmed::writelp2(char *fn)
{
	//char fn[200];
	double tt = CPU_time();
	printf("Writing LP..\n");
	//sprintf(fn, "%s-%2.2d.lp", probname, iij);
	FILE * ff = fopen(fn, "w");

	fprintf(ff, "min\n");
	int nn = 0;
	int i, j;

	for(i=0; i<m; i++)
	{
		fixed[brho[i].i] = i;
	}
	for(i=nfix1; i<m-nfix0; i++)
	{
		if(nn)
			fprintf(ff, " + ");
		fprintf(ff, " 0 y%d", brho[i].i);
		nn ++;
		if(nn%5 == 0) 
			fprintf(ff, "\n");
	}
	for(i=0; i<m; i++)
	{
		for(j=0; j<narcs[brho[i].i]; j++)
		{
			fprintf(ff, " + %.6f x%d,%d", arcs[brho[i].i][j].d, arcs[brho[i].i][j].i, brho[i].i);
			nn ++;
			if(nn%5 == 0)
				fprintf(ff, "\n");
		}
	}
	//fprintf(ff, " + %.6f z\n", constf);
	fprintf(ff, "\ns.t.\n");
	//printf("\nff\n");
	for(i=nfix1; i<m; i++)
	{
		nn = 0;
		for(j=0; j<narcs[brho[i].i]; j++)
		{
			if(nn)
				fprintf(ff, " + ");
			fprintf(ff, "x%d,%d", arcs[brho[i].i][j].i, brho[i].i);
			nn ++;
			if(nn%5 == 0)
				fprintf(ff, "\n");
		}
		if(narcs[brho[i].i] == 0)
		{
			//printf("%d\n", i);
			continue;
		}
		if(brho[i].fix == -1)
			fprintf(ff, "+ y%d = 1\n", brho[i].i);
		else
			fprintf(ff, " = 1\n");
	}
	for(i=nfix1; i<m; i++)
	{
		for(j=0; j<narcs[brho[i].i]; j++)
		{
			if(brho[fixed[arcs[brho[i].i][j].i]].fix == -1)
				fprintf(ff, "y%d - x%d,%d >= 0\n", arcs[brho[i].i][j].i, arcs[brho[i].i][j].i, brho[i].i);
		}
	}
	nn = 0;
	for(i=nfix1; i<m-nfix0; i++)
	{
		if(nn)
			fprintf(ff, " + ");
		fprintf(ff, "y%d", brho[i].i);
		nn ++;
		if(nn%5 == 0)
			fprintf(ff, "\n");
	}
	fprintf(ff, " = %d\n", p - nfix1);
	//fprintf(ff, "\nbounds\nz = 1\n");
	//nn = 0;
	/*for(i=0; i<nfix1; i++)
	{
		fprintf(ff, " y%d = 1\n", brho[i].i);
	}
	for(i=m-nfix0; i<m; i++)
	{
		fprintf(ff, " y%d = 0\n", brho[i].i);
	}*/
	fprintf(ff, "\nbinaries\n");
	nn = 0;
	for(i=nfix1; i<m-nfix0; i++)
	{
		fprintf(ff, " y%d", brho[i].i);
		nn ++;
		if(nn%5 == 0) 
			fprintf(ff, "\n");
	}
	fprintf(ff, "\nend\n");
	fclose(ff);

	if(0)
	{
		char fn1[500];
		sprintf(fn1, "%s-%2.2d.txt", probname, iij);
		ff = fopen(fn1, "w");
		nn = 0;
		for(i=0; i<m; i++)
			nn += narcs[i];
		fprintf(ff, "%ld %ld\n", m, nn);
		printf("Done in %.2f seconds\n", CPU_time() - tt);
		for(i=0; i<m; i++)
		{
			for(j=0; j<narcs[i]; j++)
			{
				fprintf(ff, "%d %d %g\n", arcs[i][j].i+1, i+1, arcs[i][j].d);
			}
		}
		fclose(ff);
	}
	iij++;
}

void Pmed::writelp()
{
	char fn[200];
	double tt = CPU_time();
	//printf("Writing LP..\n");
	sprintf(fn, "%s-%2.2d.lp", probname, iij);
	FILE * ff = fopen(fn, "w");

	fprintf(ff, "min\n");
	int nn = 0;
	int i, j;
	for(i=0; i<m; i++)
	{
		for(j=0; j<narcs[i]; j++)
		{
			if(nn)
				fprintf(ff, " + ");
			fprintf(ff, "%g x%d,%d", arcs[i][j].d, arcs[i][j].i, i);
			nn ++;
			if(nn%5 == 0)
				fprintf(ff, "\n");
		}
	}
	fprintf(ff, "\ns.t.\n");
	for(i=0; i<m; i++)
	{
		nn = 0;
		for(j=0; j<narcs[i]; j++)
		{
			if(nn)
				fprintf(ff, " + ");
			fprintf(ff, "x%d,%d", arcs[i][j].i, i);
			nn ++;
			if(nn%5 == 0)
				fprintf(ff, "\n");
		}
		if(!narcs[i])
			continue;
		if(brho[fixed[i]].fix == -1)
			fprintf(ff, "+ y%d = 1\n", i);
		else
			fprintf(ff, " = 1\n", i);
	}
	for(i=0; i<m; i++)
	{
		for(j=0; j<narcs[i]; j++)
		{
			if(brho[fixed[arcs[i][j].i]].fix == -1)
				fprintf(ff, "y%d - x%d,%d >= 0\n", arcs[i][j].i, arcs[i][j].i, i);
		}
	}
	nn = 0;
	for(i=nfix1; i<m-nfix0; i++)
	{
		if(nn)
			fprintf(ff, " + ");
		fprintf(ff, "y%d", brho[i].i);
		nn ++;
		if(nn%5 == 0)
			fprintf(ff, "\n");
	}
	fprintf(ff, " = %d\n", p - nfix1);
	//fprintf(ff, "\nbounds\n");
	//nn = 0;
	/*for(i=0; i<nfix1; i++)
	{
		fprintf(ff, " y%d = 1\n", brho[i].i);
	}
	for(i=m-nfix0; i<m; i++)
	{
		fprintf(ff, " y%d = 0\n", brho[i].i);
	}*/
	fprintf(ff, "\nbinaries\n");
	nn = 0;
	for(i=nfix1; i<m-nfix0; i++)
	{
		fprintf(ff, " y%d", brho[i].i);
		nn ++;
		if(nn%5 == 0) 
			fprintf(ff, "\n");
	}
	fprintf(ff, "\nend\n");
	fclose(ff);
	//printf("Done in %.2f seconds\n", CPU_time() - tt);

	if(0)
	{
		sprintf(fn, "%s-%2.2d.txt", probname, iij);
		ff = fopen(fn, "w");
		nn = 0;
		for(i=0; i<m; i++)
		
			nn += narcs[i];
		fprintf(ff, "%ld %ld\n", m, nn);
		printf("Done in %.2f seconds\n", CPU_time() - tt);
		for(i=0; i<m; i++)
		{
			for(j=0; j<narcs[i]; j++)
			{
				fprintf(ff, "%d %d %g\n", arcs[i][j].i+1, i+1, arcs[i][j].d);
			}
		}
		fclose(ff);
	}
	iij++;
}


void Pmed::writefix(char * fn)
{
	FILE *ff = fopen(fn, "w");
	int i = 0;
	fprintf(ff, "%ld %ld\n", m, p);
	fprintf(ff, "nfix1 = %d\n", nfix1);
	for(i=0; i<nfix1; i++)
	{
		fprintf(ff, "%d\n", brho[i].i);
	}
	fprintf(ff, "nfix0 = %d\n", nfix0);
	for(i=m - nfix0; i<m; i++)
	{
		fprintf(ff, "%d\n", brho[i].i);
	}
	fclose(ff);
}

void Pmed::writepp(char * fn)
{
	FILE *ff = fopen(fn, "w");
	int nn = 0, i = 0, j = 0;
	for(i=0; i<m; i++)
		nn += narcs[i];
	fprintf(ff, "%ld %ld\n", m, nn);
	for(i=0; i<m; i++)
	{
		for(j=0; j<narcs[i]; j++)
		{
			fprintf(ff, "%d %d %.6f\n", arcs[i][j].i+1, i+1, arcs[i][j].d);
		}
	}
	fclose(ff);
}

void Pmed::primal()
{
	int *sol = new int[m];
	int k, h, i;
	double ub = 0;
	int nas = 0;

	for(k=0 ; k<m ; k++)
	{
		if (flag[k] == 1) 
		{
			sol[k] = k;
			ub += c[k];
			continue;
		}

		for(h=0; h<narcs[k]; h++)
		{
			i = arcs[k][h].i;
			if(flag[i])
			{
				ub += arcs[k][h].d;
				sol[k] = i;
				nas ++;
				break;
			}
		}
	}

	if(nas == m-p && ub < bub)
	{
		bub = ub;
		memcpy(bsol, sol, m*sizeof(int));
	}
	delete[] sol;
}
void Pmed::primalBB()
{
	int *sol = new int[m];
	int k, h, i;
	double ub = 0;
	int nas = 0;

	for(k=0 ; k<m ; k++)
	{
		if (flag[k] == 1) 
		{
			sol[k] = k;
			ub += c[k];
			continue;
		}

		for(h=0; h<narcs[k]; h++)
		{
			i = arcs[k][h].i;
			if(flag[i])
			{
				ub += arcs[k][h].d;
				sol[k] = i;
				nas ++;
				break;
			}
		}
	}

	if(nas == m-p && ub < bub)
	{
		bub = ub;
		memcpy(bsol, sol, m*sizeof(int));
		printf("BUB = %.2f\n", bub);
	}
	delete[] sol;
}


void Pmed::prim(int *ss)
{
	int *sol = new int[m];
	int k, h, i;
	double ub = 0;
	int nas = 0;
	/*for(i=0; i<m; i++)
	{
		printf("%d",ss[brho[i].i]);
	}
	printf("\n");*/
	for(k=0 ; k<m ; k++)
	{
		if (ss[k] == 1) 
		{
			sol[k] = k;
			ub += c[k];
			continue;
		}

		for(h=0; h<narcs[k]; h++)
		{
			i = arcs[k][h].i;
			if(ss[i])
			{
				ub += arcs[k][h].d;
				sol[k] = i;
				nas ++;
				break;
			}
		}
	}

	if(nas == m-p && ub < bub)
	{
		bub = ub;
		memcpy(bsol, sol, m*sizeof(int));
		printf("\nbub=%12.2f\n", bub);
	}
	delete[] sol;
}




void Pmed::writearcs()
{
	int i, j, k;

	FILE * ff = fopen("ar.txt", "w");

	for(j=0; j<m; j++)
	{
		fprintf(ff, "j=%4d\n", j);
		for(k=0; k<narcs1[j]; k++)
		{
			i = arcs1[j][k].i;
			fprintf(ff, "%4d%10.2f%6d%6d", i, arcs1[j][k].d, rho[fixed[i]].i, rho[fixed[i]].fix);
			if(k>=narcs[j])
			{
				fprintf(ff, "\n");
				continue;
			}
			i = arcs[j][k].i;
			fprintf(ff, " -- %4d%10.2f%6d%6d\n", i, arcs[j][k].d, rho[fixed[i]].i, rho[fixed[i]].fix);
		}
	}
	fclose(ff);

}

void Pmed::fix(double lb_)
{
	//qsort(rho + nfix1, m - nfix0 - nfix1, sizeof(Arcr), compareArcrAZ1);
	int i, j;
	for(i=nfix1; i<p; i++)
	{
		if(lb_ + rho[p].d - rho[i].d < bub + EPSILON)
			break;
		rho[i].fix = 1;
		narcs[rho[i].i] = 0;
		nfix1 ++;
	}
	if(nfix1 == p)
	{
		for(i=p; i<m-nfix0; i++)
		{
			rho[i].fix = 0;
			nfix0 ++;
		}
		return;
	}
	int nf = nfix0;
	for(i=m-nfix0-1; i>=p; i--)
	{
		if(lb_ + rho[i].d - rho[p-1].d < bub + EPSILON)
			break;
		rho[i].fix = 0;
		nfix0 ++;
	}
	if(1 && nf < nfix0)
	{
		int k;
		for(j=0; j<m; j++)
		{
			k = 0;
			while(1)
			{
				if(k >= narcs[j])
					break;
				i = arcs[j][k].i;
				if(rho[fixed[i]].fix == 0)
				{
					if(k < narcs[j] - 1)
						shift(arcs[j] + k, narcs[j] - k);
					narcs[j] --;
				}
				else
				{
					k++;
				}
			}
		}
	}
	//int nfixed = nfix0 + nfix1;
	//if(nfixed)
	//{
	//	printf("\nnfixed = %d\n", nfixed);
	//}
}

void Pmed::dualvalue()
{
	int i, j, h, k;

	for(i=0; i<m; i++)
	{
		flag[i] = false;
		if(0)
		{
			rho[i].i = i;
			rho[i].fix = -1;
			if(narcs[i] == 0)
				rho[i].fix = 1;
			fixed[i] = i;
		}
		fixed[rho[i].i] = i;
		rho[i].d = c[rho[i].i] - lm[rho[i].i];
	}

	double lbl = 0.;
	for(j=0; j<m ; j++)
	{
		lbl += lm[j];
		//printf("LM[%d] = %10.2f lbl = %10.2f\n", j, lm[j], lbl);
		if(rho[fixed[j]].fix == 1)
			continue;
		for(h = 0; h<narcs[j]; h++)
		{
			i = arcs[j][h].i;
			if(arcs[j][h].d - lm[j] >= -EPSILON)
				break;
			if(rho[fixed[i]].fix == 0)
				continue;
			rho[fixed[i]].d += (arcs[j][h].d - lm[j]);
		}
	}

	

	//qsort(rho, m, sizeof(Arcr), compareArcrAZ);
	qsort(rho + nfix1, m - nfix0 - nfix1, sizeof(Arcr), compareArcrAZ1);
	//selectp(rho, nfix1, m-nfix0-1, p-nfix1);
	for(i=0; i<m; i++)
	{
		fixed[rho[i].i] = i;
	}
	if(0)
	{
		for(i=0; i<m; i++)
		{
			printf("%5d %5d %10.1f %10.2f\n", rho[i].i, rho[i].fix, rho[i].d, lm[rho[i].i]);
		}
		writelp();
	}
	//selectp(rho, 0, m-1, p);
	
	lb = lbl;
	
	for(i=0; i<p; i++)
	{
		flag[rho[i].i] = true;
		lb += rho[i].d;
	}
	
	fiter ++;
	if(lb > blb)
	{
		fiter = 0;
		blb = lb;
		memcpy(blm, lm, m*sizeof(double));
		memcpy(brho, rho, m*sizeof(Arcr));
	}
}

void Pmed::dualvalueBB()
{
	int i, j, h, k;

	for(i=0; i<m; i++)
	{
		flag[i] = false;
		fixed[rho[i].i] = i;
		rho[i].d = c[rho[i].i] - lm[rho[i].i];
	}

	double lbl = 0.;
	for(j=0; j<m ; j++)
	{
		lbl += lm[j];
		//printf("LM[%d] = %10.2f lbl = %10.2f\n", j, lm[j], lbl);
		if(rho[fixed[j]].fix == 1)
			continue;
		for(h = 0; h<narcs[j]; h++)
		{
			i = arcs[j][h].i;
			if(arcs[j][h].d - lm[j] >= -EPSILON)
				break;
			if(rho[fixed[i]].fix == 0)
				continue;
			rho[fixed[i]].d += (arcs[j][h].d - lm[j]);
		}
	}

	

	//qsort(rho, m, sizeof(Arcr), compareArcrAZ);
	qsort(rho + nfix1, m - nfix0 - nfix1, sizeof(Arcr), compareArcrAZ1);
	//selectp(rho, nfix1, m-nfix0-1, p-nfix1);
	for(i=0; i<m; i++)
	{
		fixed[rho[i].i] = i;
	}
	if(0)
	{
		for(i=0; i<m; i++)
		{
			printf("%5d %5d %10.1f %10.2f\n", rho[i].i, rho[i].fix, rho[i].d, lm[rho[i].i]);
		}
		writelp();
	}
	//selectp(rho, 0, m-1, p);
	
	lb = lbl;
	
	for(i=0; i<p; i++)
	{
		flag[rho[i].i] = true;
		lb += rho[i].d;
	}

	if(lb > blb)
		blb = lb;
	
	fiter ++;
}


void Pmed::subgradient()
{
	int i, j, k, h;
	for(k=0; k<m; k++)
	{
		j = rho[k].i;
		sg[j] = 1.0 - flag[j];
		if(rho[k].fix == 1)
	
			continue;
		for(h=0; h<narcs[j]; h++)
		{
			i = arcs[j][h].i;
			if(arcs[j][h].d - lm[j] >= -EPSILON) 
				break;
			if(rho[fixed[i]].fix == 0)
				continue;
			if(flag[i] == 1) 
				sg[j]--;
		}
	}
	double nSG = 0.0;
	for(i=0; i<m; i++)
		nSG += sg[i] * sg[i];
	normSG = nSG;
	//printf("\SG\n");
	//for(i=0; i<m; i++)
	//{
	//	printf("%5d %10.2f\n", i, sg[i]);
	//}
}



void Pmed::step()
{
	// update phi
	phi/=phiFactor;
	theta = phi * ( 1.05 * bub - lb) / normSG;
	for(int i=0; i<m; i++)
	{
		lm[i] += theta * sg[i];
		if(lm[i] > lmub[i])
			lm[i] = lmub[i];
	}
}

void Pmed::step1()
{
	// update phi
	if(fiter >= cfiter)
	{
		phi/=phiFactor;
		fiter = 0;
	}
	theta = phi * ( 1.05 * bub - lb) / normSG;

	for(int i=0; i<m; i++)
	{
		lm[i] += theta * sg[i];
		if( ubi[i] >= narcs[i] - 1)
			continue;
		if(lm[i] >= arcs[i][ubi[i]].d)
		{
			lm[i] = arcs[i][ubi[i]].d;
			ubi[i] ++;
			if(ubi[i] >= narcs[i])
				ubi[i] = narcs[i] - 1;
				//addarcs(i);
			//lbi[i] ++;
			continue;
		}
		if(lbi[i] == -1)
		{
			//lm[i] = 0.;
			continue;
		}
		if(lm[i] <= arcs[i][lbi[i]].d)
		{
			lm[i] = arcs[i][lbi[i]].d;
			lbi[i] --;
			//ubi[i] --;
		}
	}

	//printf("\nlm\n");
	//for(int i=0; i<m; i++)
	//{
	//	printf("%5d %10.2f\n", i, lm[i]);
	//}
}

void Pmed::step2()
{
	if(fiter == 30)
	{
		Eta *= 0.95;
	}
	for(int i=0; i<m; i++)
	{
		lm[i] += Eta * sg[i];
		if(lm[i] >= arcs[i][ubi[i]].d)
		{
			lm[i] = arcs[i][ubi[i]].d;
			ubi[i] ++;
			if(ubi[i] >= narcs[i])
				addarcs(i);
			//lbi[i] ++;
			continue;
		}
		if(lbi[i] == -1)
			continue;
		if(lm[i] <= arcs[i][lbi[i]].d)
		{
			lm[i] = arcs[i][lbi[i]].d;
			lbi[i] --;
			//ubi[i] --;
		}
	}
}

void Pmed::addarcs(int j)
{
	error("Exit in addarcs. Sorry.");
	//printf("\n@");
	int i;
	int ml = m2;
	if(ml + narcs[j] > m)
		ml = m - narcs[j];
	
	Arc *col = new Arc[m];

	for(i=0; i<m; i++)
	{
		col[i].i = i;
		if(i==j)
			col[i].d = INFTY;
		else
			col[i].d = getdist(i, j);
	}
	qsort(col, m, sizeof(Arc), compareArcAZ);

	delete[] arcs[j];
	narcs[j] += ml;
	arcs[j] = new Arc[narcs[j]];
	memcpy(arcs[j], col, narcs[j] * sizeof(Arc));
	
	delete[] col;
}

void Pmed::lagrangean()
{
	if(dispInterval)
	{
		printf("*************************************************************\n");
		printf("**                                                         **\n");
		printf("**       Lagrangian heuristic for p-Median promlem         **\n");
		printf("**       authors: Igor Vasiliev (vil@icc.ru)               **\n");
		printf("**                                                         **\n");
		printf("*************************************************************\n");
		printf("%-20s%5d\n", "Number of nodes:", m);
		printf("%-20s%5d\n", "Number of medians:", p);
	}
	printf("Reading distances\n");
	char fn[200];
	sprintf(fn, "%s.dis", probname);
	fd = fopen(fn, "rb");
	if(!fd)
		error("Cannot open file ", fn);

	int i;
	// allocate and initialize lagrangean multipliers
	lm = new double[m];
	blm = new double[m];
	lmub = new double[m];
	lbi = new int[m];
	ubi = new int[m];
	for(i=0; i<m; i++)
	{
		lbi[i] = -1;
		ubi[i] = 1;
		lm[i] = arcs[i][0].d;
		lmub[i] = arcs[i][narcs[i]-1].d;
	}
	// allocate subgradiend
	sg = new double[m];

	// start time
	double btime = CPU_time();

	// set algorithm parameters
	phi = 2.;
	phiFactor = 1.5;
	theta = 1.;
	cfiter = 10;

	// initialize bounds
	blb = -1.e30;
	// initialize best medians
	flag = new int[m];
	// initialize reduced costs
	rho = new Arcr[m];
	brho = new Arcr[m];
	
	// iteration counter
	int it = 0;
	fiter = 0;

	// main loop
	while(true)
	{	
		// evaluate lagrangean function
		dualvalue();
		if(it==0)
		{
			Eta = (bub - blb);
		}
		//primal();
		if(dispInterval && it%dispInterval == 0)
		{
			printf("%6d %15.2f %12.2f %12.2f %12.2f\n", it, blb, bub, 
				(bub - blb) / bub * 100., CPU_time() - btime);
		}
		if((bub - blb) / bub < 3.01)
		{
			phiFactor = 1.1;
			cfiter = 20;
		}
		// check optimality
		if( blb/bub >= 1. - DUALGAP)
		{
			if(dispInterval)
				printf("UB=LB\n");
			break;
		}
		subgradient();
		// check optimality
		if(normSG < EPSILON)
		{
			if(dispInterval)
					printf("nornSG = 0\n");
			//bestUB = ub = lb;
			break;
		}
		// go ahead subgradient
		step1();
		// check stop condition
		if(phi < THETA_MIN)
		{
			if(dispInterval)
				printf("Stop by theta.\n");
			break;
		}
		if(Eta < THETA_MIN)
		{
			if(dispInterval)
				printf("Stop by theta.\n");
			break;
		}
		it++;
	}
	if(dispInterval)
	{
		printf("%6d %15.2f %12.2f %12.2f %12.2f\n", it, blb, bub, 
				(bub - blb) / bub * 100., CPU_time() - btime);
	}

	//writecore();

	fclose(fd);
	//free memory
	delete[] lm;
	delete[] blm;
	delete[] lmub;
	delete[] sg;
	delete[] rho;
	delete[] brho;
	delete[] flag;
	delete[] lbi;
	delete[] ubi;
}

void Pmed::shift(Arc* aa, int nl)
{
	memmove(aa, aa+1, (nl-1) * sizeof(Arc));
}

void Pmed::lang_init()
{
	int i;
	// allocate and initialize lagrangean multipliers
	lm = new double[m];
	blm = new double[m];
	lmub = new double[m];
	lbi = new int[m];
	ubi = new int[m];
	fixed = new int[m];
	rho = new Arcr[m];

	for(i=0; i<m; i++)
	{
		if(narcs[i] > 0)
		{
			lbi[i] = -1;
			ubi[i] = 1;
			lm[i] = arcs[i][0].d;
			blm[i] = arcs[i][0].d;
			lmub[i] = arcs[i][narcs[i]-1].d;
		}
		else
		{
			lbi[i] = -1;
			ubi[i] = -1;
			lm[i] = 0.;
			blm[i] = 0.;
			lmub[i] = 0.;
		}
	}
	// allocate subgradiend
	sg = new double[m];

	// start time
	double btime = CPU_time();

	// initialize bounds
	blb = -1.e30;
	// initialize best medians
	flag = new int[m];
	// initialize reduced costs
	brho = new Arcr[m];
}

void Pmed::lang_final()
{
	
	//free memory
	delete[] lm;
	delete[] blm;
	delete[] lmub;
	delete[] sg;
	delete[] rho;
	delete[] brho;
	delete[] flag;
	delete[] lbi;
	delete[] ubi;
	delete[] fixed;
}

void Pmed::lang_set1()
{
	// set algorithm parameters
	int i, j;
	nnfix1 = 0;
	nnfix0 = 0;
	for(i=0; i<m; i++)
	{
		rho[i].i = i;
		rho[i].fix = -1;
		if(narcs[i] == 0)
		{
			rho[i].fix = 1;
			nnfix1 ++;
		}
	}

	nfix0 = nnfix0;
	nfix1 = nnfix1;



	qsort(rho, m, sizeof(Arcr), compareArcrAZ);
	for(i=0; i<m; i++)
	{
		fixed[rho[i].i] = i;
	}
	delete[] nlarcs;
}

void Pmed::lang_set2()
{
	int i;
	if(narcs)
		delete[] narcs;
	if(arcs)
	{
		for(int i=0; i<m; i++)
			delete[] arcs[i];
		delete[] arcs;
	}
	arcs = new Arc*[m];
	narcs = new int[m];
	memcpy(narcs, narcs1, m * sizeof(int));
	for(i=0; i<m; i++)
	{
		arcs[i] = new Arc[narcs1[i]];
		memcpy(arcs[i], arcs1[i], narcs1[i] * sizeof(Arc));
	}
	
	nfix0 = nnfix0;
	nfix1 = nnfix1;
	for(i=nfix1; i< m-nfix0; i++)
		rho[i].fix = -1;
}
void Pmed::lang_set3()
{
	int i;
	if(narcs)
		delete[] narcs;
	if(arcs)
	{
		for(int i=0; i<m; i++)
			delete[] arcs[i];
		delete[] arcs;
	}
	arcs = new Arc*[m];
	narcs = new int[m];
	memcpy(narcs, narcs1, m * sizeof(int));
	for(i=0; i<m; i++)
	{
		arcs[i] = new Arc[narcs1[i]];
		memcpy(arcs[i], arcs1[i], narcs1[i] * sizeof(Arc));
	}
	
	int k = p - p0;

	nfix0 = nfix0a[k];
	nfix1 = nfix1a[k];
	bub = buba[k];
	blb = blba[k];
	memcpy(rho, brhoa[k], m * sizeof(Arcr));
	memcpy(lm, blma[k], m * sizeof(double));
	memcpy(brho, brhoa[k], m * sizeof(Arcr));
	memcpy(blm, blma[k], m * sizeof(double));
	for(i=0; i<m; i++)
	{
		fixed[rho[i].i] = i;
	}
	
	memcpy(lbi, lbia[k], m * sizeof(int));
	memcpy(ubi, ubia[k], m * sizeof(int));

	if(1)
	{
		for(i=0; i<nfix1; i++)
		{
			narcs[brho[i].i] = 0;
		}
		int j;
		for(j=0; j<m; j++)
		{
			k = 0;
			while(1)
			{
				if(k >= narcs[j])
					break;
				i = arcs[j][k].i;
				if(rho[fixed[i]].fix == 0)
				{
					if(k < narcs[j] - 1)
						shift(arcs[j] + k, narcs[j] - k);
					narcs[j] --;
				}
				else
				{
					k++;
				}
			}
		}
	}
}

void Pmed::lang_set1c()
{
	// set algorithm parameters

	if(!isfull)
	{
		m1 = (int) (MEMSIZE / (long long)sizeof(Arc) / (long long)m);
		if(m1 >= m)
			m1 = m - 1;
		m2 = 100;
	}

	/*phi = 2;
	stopphi = .005;
	phiFactor = 2;
	cfiter = 30;*/

	phi = 1.5;
	stopphi = .001;
	phiFactor = 1.01;
	cfiter = 1;

	EPSCORE = 1.e30;
	
	/*MAXNODE = 8;
	MAXCORE = 12*m;
	MINARCS = 3;*/
	//CORETIME = 60.;

	PRIMALFREQ = 2000000000;
}




void Pmed::lang()
{
	if(dispInterval)
	{
		printf("\n*************************************************************\n");
		printf("**                                                         **\n");
		printf("**       Lagrangian heuristic for p-Median promlem         **\n");
		printf("**       authors: Igor Vasiliev (vil@icc.ru)               **\n");
		printf("**                                                         **\n");
		printf("*************************************************************\n");
		printf("%-20s%5d\n", "Number of nodes:", m);
		printf("%-20s%5d\n", "Number of medians:", p);
	}
	
	double btime = CPU_time();

	// iteration counter
	int it = 0;
	fiter = 0;

	//double stopphi = 0.001;
	//memcpy(lm, blm, m*sizeof(double));

	// main loop
	
	while(true)
	{	
		// evaluate lagrangean function
		dualvalue();
		if(PRIMALFREQ && it%PRIMALFREQ == 0)
		{
			primal();
		}
		if(FIXFREQ && it%FIXFREQ == 0)
		{
			primal();
			fix(lb);
			if(nfix1 == p)
			{
				blb = bub;
				if(dispInterval)
					printf("All variables are fixed.\n");

				break;
			}
		}
		if(it==0)
		{
			Eta = (bub - blb);
		}
		if(dispInterval && it%dispInterval == 0)
		{
			printf("%6d %15.2f %12.2f %12.2f %12.2f %g %d\n", it, blb, bub, 
				(bub - blb) / bub * 100., CPU_time() - btime, phi, nfix1+nfix0);
		}
		// check optimality
		if( blb/bub >= 1. - DUALGAP)
		{
			if(dispInterval)
				printf("UB=LB\n");
			break;
		}
		subgradient();
		// check optimality
		if(normSG < EPSILON)
		{
			if(dispInterval)
					printf("nornSG = 0\n");
			primal();
			dualvalue();
			//bestUB = ub = lb;
			break;
		}
		// go ahead subgradient
		step1();
		// check stop condition
		if(phi < stopphi)
		//if(theta < stopphi)
		{
			if(dispInterval)
				printf("Stop by theta.\n");
			break;
		}
		if(Eta < THETA_MIN)
		{
			if(dispInterval)
				printf("Stop by theta.\n");
			break;
		}
		it++;
	}

	if(dispInterval)
	{
		printf("%6d %15.2f %12.2f %12.2f %12.2f\n", it, blb, bub, 
				(bub - blb) / bub * 100., CPU_time() - btime);
	}

	//ttime = CPU_time() - btime;

	int i;
	for(i=0; i<m; i++)
	{
		brho[i].fix = rho[fixed[brho[i].i]].fix;
	}
	qsort(brho, m, sizeof(Arcr), compareArcrAZ2);
	

	if(1)
	{
		int k = p-p0;
		nfix1a[k] = nfix1;
		nfix0a[k] = nfix0;
		memcpy(brhoa[k], brho, m * sizeof(Arcr));
		memcpy(blma[k], blm, m * sizeof(double));
		memcpy(lbia[k], lbi, m * sizeof(int));
		memcpy(ubia[k], ubi, m * sizeof(int));
		buba[k] = bub;
		blba[k] = blb;
	}

}

void Pmed::langBB()
{
	if(dispInterval)
	{
		printf("\n*************************************************************\n");
		printf("**                                                         **\n");
		printf("**       Lagrangian heuristic for p-Median promlem         **\n");
		printf("**       authors: Igor Vasiliev (vil@icc.ru)               **\n");
		printf("**                                                         **\n");
		printf("*************************************************************\n");
		printf("%-20s%5d\n", "Number of nodes:", m);
		printf("%-20s%5d\n", "Number of medians:", p);
	}
	
	double btime = CPU_time();

	// iteration counter
	int it = 0;
	fiter = 0;

	//double stopphi = 0.001;
	//memcpy(lm, blm, m*sizeof(double));

	// main loop
	
	while(true)
	{	
		// evaluate lagrangean function
		dualvalueBB();
		if(PRIMALFREQ && it%PRIMALFREQ == 0)
		{
			primalBB();
		}
		if(FIXFREQ && it%FIXFREQ == 0)
		{
			primalBB();
			//fix(lb);
			if(nfix1 == p)
			{
				blb = bub;
				if(dispInterval)
					printf("All variables are fixed.\n");
				break;
			}
		}
		if(it==0)
		{
			Eta = (bub - blb);
		}
		if(dispInterval && it%dispInterval == 0)
		{
			printf("%6d %15.2f %12.2f %12.2f %12.2f %g %d\n", it, blb, bub, 
				(bub - blb) / bub * 100., CPU_time() - btime, phi, nfix1+nfix0);
		}
		// check optimality
		if( blb/bub >= 1. - DUALGAP)
		{
			if(dispInterval)
				printf("UB=LB\n");
			break;
		}
		subgradient();
		// check optimality
		if(normSG < EPSILON)
		{
			if(dispInterval)
					printf("nornSG = 0\n");
			primalBB();
			dualvalue();
			//bestUB = ub = lb;
			break;
		}
		// go ahead subgradient
		step1();
		// check stop condition
		if(phi < stopphi)
		//if(theta < stopphi)
		{
			if(dispInterval)
				printf("Stop by theta.\n");
			break;
		}
		if(Eta < THETA_MIN)
		{
			if(dispInterval)
				printf("Stop by theta.\n");
			break;
		}
		it++;
	}

	if(dispInterval)
	{
		printf("%6d %15.2f %12.2f %12.2f %12.2f\n", it, blb, bub, 
				(bub - blb) / bub * 100., CPU_time() - btime);
	}

	//ttime = CPU_time() - btime;


	if(0)
	{
		int k = p-p0;
		nfix1a[k] = nfix1;
		nfix0a[k] = nfix0;
		memcpy(brhoa[k], rho, m * sizeof(Arcr));
		memcpy(blma[k], lm, m * sizeof(double));
		memcpy(lbia[k], lbi, m * sizeof(int));
		memcpy(ubia[k], ubi, m * sizeof(int));
		buba[k] = bub;
		blba[k] = blb;
	}

}



void Pmed::makedist()
{
	printf("Making distances\n");
	double tt = CPU_time();

	int i, j;

	Arc *col = new Arc[m];
	
	arcs = new Arc*[m];
	narcs = new int[m];

	for(i=0; i<m; i++)
	{
		if(i%1000 == 0)
			printf("%5.5d done\n", i);
		for(j=0; j<m; j++)
		{
			col[j].i = j;
			if(i==j)
				col[j].d = INFTY;
			else
				col[j].d = getdist(i, j);
		}
		qsort(col, m, sizeof(Arc), compareArcAZ);

		if(0)
		{
			printf("hh %d\n", i);
			for(j=0; j<m; j++)
				printf("%.2f\n", col[j].d);
		}

		narcs[i] = m1;

		arcs[i] = new Arc[narcs[i]];
		memcpy(arcs[i], col, m1 * sizeof(Arc));

		
	}

	if(0)
	{
		for(i=0; i<m; i++)
		{
			printf("ff = %d\n", i);
			for(j=0; j<narcs[i]; j++)
			{
				printf("%d %.2f\n", arcs[i][j].i, arcs[i][j].d);
			}
		}
	}

	delete[] col;
	makelarcs();
	printf("Done in %.2f seconds\n", CPU_time() - tt);
}

void Pmed::readdist()
{
	error("Soory, exit readdist");
	/*printf("Reading distances\n");
	char fn[200];
	sprintf(fn, "%s.dis", probname);
	FILE *ff = fopen(fn, "rb");
	if(!ff)
		error("Cannot open file ", fn);
	double tt = CPU_time();

	int i, j;

	arcs = new Arc*[m];
	narcs = new int[m];

	for(i=0; i<m; i++)
	{
		if(i%1000 == 0)
			printf("%5.5d done\n", i);
		narcs[i] = m1;
		arcs[i] = new Arc[narcs[i]];
		__int64 pos = (__int64) i * (__int64) m * (__int64) sizeof(Arc);
		_fseeki64(ff, pos, SEEK_SET);
		fread(arcs[i], sizeof(Arc), m1, ff);
		//printf("\n %d %d %f %f %f\n", m, m1, arcs[i][0].d, arcs[i][1].d, arcs[i][m1-1].d);
	}

	if(0)
	{
		for(i=0; i<m; i++)
		{
			printf("ff = %d\n", i);
			for(j=0; j<narcs[i]; j++)
			{
				printf("%d %.2f\n", arcs[i][j].i, arcs[i][j].d);
			}
		}
	}
	fclose(ff);
	//makelarcs();
	printf("Done in %.2f seconds\n", CPU_time() - tt);
	*/
}

void Pmed::writedist()
{

	printf("Writing distances\n");
	double tt = CPU_time();


	int i;

	arcs = new Arc*[m];
	narcs = new int[m];

	for(i=0; i<m; i++)
		arcs[i] = new Arc[m1];
	
	double ub = 0;
	
#pragma omp parallel for
	for(i=0; i<m; i++)
	{
		int j;
//		int id = omp_get_thread_num();
//		if(i%1000 == 0)
//			printf("%5.5d done by %d\n", i, id);
		Arc *col = new Arc[m];
		for(j=0; j<m; j++)
		{
			col[j].i = j;
			if(i==j)
				col[j].d = INFTY;
			else
				col[j].d = getdist(i, j);
		}
		qsort(col, m, sizeof(Arc), compareArcAZ);
		narcs[i] = m1;
		memcpy(arcs[i], col, m1 * sizeof(Arc));
		delete[] col;
	}
	
	printf("Done in %.2f seconds\n", CPU_time() - tt);

}


struct Arc1;

void Pmed:: readprob(int m_, int *narcs_, Arc1 **arcs_)
{
	int i, j;
	m = m_;
	n = 0;
	narcs1 = narcs_;
	arcs1 = (Arc**) arcs_;
	arcs = new Arc*[m];
	narcs = new int[m];
	memcpy(narcs, narcs1, m * sizeof(int));
	for(i=0; i<m; i++)
	{
		arcs[i] = new Arc[narcs1[i]];
		memcpy(arcs[i], arcs1[i], narcs1[i] * sizeof(Arc));
	}
	
	c = new double[m];
	int nfixed = 0;
	for(i=0; i<m; i++)
	{
		c[i] = 0;
	}
	sprintf(probname, "odm");

	int k = p1 - p0 + 1;

	nfix1a = new int[k];
	nfix0a = new int[k];
	brhoa = new Arcr*[k];
	blma = new double*[k];
	lbia = new int*[k];
	ubia = new int*[k];
	buba = new double[k];
	blba = new double[k];

	for(i=0; i<k; i++)
	{
		brhoa[i] = new Arcr[m];
		blma[i] = new double[m];
		lbia[i] = new int[m];
		ubia[i] = new int[m];
	}

}

void Pmed::par1()
{
	phi = 2.;
	stopphi = .001;
	phiFactor = 1.01;
	cfiter = 1;

	DUALGAP = 0.01;
	PRIMALFREQ = 1;
	FIXFREQ = 10;

	dispInterval = 0;
}

void Pmed::par2()
{
	phi = 2.;
	stopphi = .0001;
	phiFactor = 1.001;
	cfiter = 1;

	DUALGAP = 0.0001;
	PRIMALFREQ = 1;
	FIXFREQ = 10;

	dispInterval = 0;
}

void Pmed::par3()
{
	phi = 2.;
	stopphi = .00001;
	phiFactor = 1.0001;
	cfiter = 1;

	DUALGAP = 0.0001;
	PRIMALFREQ = 1;
	FIXFREQ = 10;

	dispInterval = 0;
}

void Pmed::par4()
{
	phi = .1;
	stopphi = .0001;
	phiFactor = 1.001;
	cfiter = 1;

	DUALGAP = 0.0001;
	PRIMALFREQ = 1;
	FIXFREQ = 10;

	dispInterval = 0;
}

int Pmed::optimize(int p_, double ub_, int *sol_)
{
	ttime =  CPU_time();
	p = p_;
	bub = ub_;
	bsol = sol_;
	blb = -1.e30;

	par1();

	lang();

	ttime =  CPU_time() - ttime;
	return 0;
}

int Pmed::reopt(int p_, int *sol_)
{
	int i;

	ttime =  CPU_time();
	p = p_;
	bsol = sol_;
	
	lang_set3();

	lang();

	ttime =  CPU_time() - ttime;

	return 0;
}

void CPXerror(int errorcode, char *msg)
{
	if(errorcode)
	{
		printf("\nCplex error with code %d\n%s\n", errorcode, msg);
		abort();
	}
}

int Pmed::solveCPX(int p_, int *sol_)
{
	int i;
	
	ttime =  CPU_time();
	p = p_;
	bsol = sol_;
	
	lang_set3();
	/*for(i=0; i<m; i++)
	{
		fixed[brho[i].i] = i;
	}
*/
	//if((bub-blb)/bub > 0.0001)
	if(0)
	{
		for(i=0; i<SATRY; i++)
			SA();

		//lang();
		fix(bub);
		memcpy(brho, rho, m * sizeof(Arcr));
		for(i=0; i<m; i++)
		{
			fixed[brho[i].i] = i;
		}
		char ss[500];
		sprintf(ss, "bub = %14.2f fixed1 = %d fixed0 = %d nfix = %d\n", bub, nfix1, nfix0, nfix1 + nfix0);
		print(ss);
	}
	
	char fnn[500];
	sprintf(fnn, "p%d.lp", omp_get_thread_num());
	writelp(fnn);
	// Open Cplex environment
	int status = 0;
	CPXENVptr env;
	env = CPXopenCPLEX (&status);
	CPXerror(status, "env = CPXopenCPLEX (&status);");

	// Create Cplex problem
	CPXLPptr  lp;
	lp = CPXcreateprob (env, &status, probname);
	CPXerror(status, "CPXcreateprob in main");
	// Read lp file
	CPXerror(CPXreadcopyprob(env, lp, fnn, 0), "Cannon fide lp file");
	CPXsetintparam(env, CPX_PARAM_THREADS, 1);
	CPXmipopt(env, lp);
	
	//double bub = 1.e30;
	CPXerror(CPXgetmipobjval(env, lp, &bub), "CPXgetmipobjval(env, lp, &bub)");
	double *x = new double[m-nfix0-nfix1];
	CPXerror(CPXgetmipx(env, lp, x, 0, m-nfix0-nfix1-1), "CPXgetmipobjval(env, lp, &bub)");
	//CPXerror(CPXsolution(env, lp, 0, &bub, x, 0, 0, 0), "CPXsolution in cutting plane");
	blb = bub;
	int nnodes = CPXgetnodecnt(env, lp);
	if(CPXgetnodeleftcnt(env, lp))
		CPXerror(CPXgetbestobjval(env, lp, &blb), "CPXgetbestobjval(env, lp, &blb)");
	CPXerror(CPXfreeprob(env, &lp), "CPXfreeprob(env, &lp)");
	CPXerror(CPXcloseCPLEX (&env), "CPXcloseCPLEX (&env)");

	int *flag2 = new int[m];
	for(i=0; i<m; i++)
		flag2[i] = 0;
	for(i=0; i<nfix1; i++)
		flag2[brho[i].i] = 1;
	for(i=nfix1; i<m-nfix0; i++)
	{
		if(x[i-nfix1] < 0.5)
			continue;
		flag2[brho[i].i] = 1;
	}


	double dff = 0.;
	int j;
	for(i=0; i<m; i++)
	{
		if(flag2[i])
		{
			bsol[i] = i;
			continue;
		}
		for(j=0; j<narcs[i]; j++)
		{
			int k = arcs[i][j].i;
			if(flag2[k])
			{
				bsol[i] = k;
				dff += arcs[i][j].d;
				break;
			}
		}
	}
	if(fabs(bub - dff) > EPSILON)
	{
		printf("%12.5f  == %12.5f\n", bub, dff);
		writelp2("ss.lp");
		writefix("ss.fix");
		writepp("ss.txt");
		error("bub != dff");
	}

	delete[] x;
	return 0;
}

int indff = 10;

int Pmed::optBB(int p_, int *sol_)
{
	int i;
	char sss[500];
	ttime =  CPU_time();
	p = p_;
	bsol = sol_;
	
	lang_set3();
	for(i=0; i<m; i++)
	{
		fixed[brho[i].i] = i;
	}
	
	int *ss = new int[m];
	for(i=0; i<m; i++)
	{
		ss[i] = 0;
	}
	for(i=0; i<p; i++)
	{
		ss[brho[i].i] = 1;
	}

	char fn[500];
	for(i=0; i<SATRY; i++)
		SA();
	sprintf(sss, "bub = %14.2f fixed1 = %d fixed0 = %d nfix = %d\n", bub, nfix1, nfix0, nfix1 + nfix0);
	print(sss);
	memcpy(rho, brho, m * sizeof(Arcr));
	memcpy(lm, blm, m * sizeof(double));
	//globalfix();

	indff++;

	solveCPX(p_, sol_);

	delete[] ss;
	ttime =  CPU_time() - ttime;
	return 0;
}



void Pmed::globalfix()
{
	char sss[500];
	print("Global fixing:\n");
	double ttime = CPU_time();
	int i, j, h, l;
	int nnarcs = 0;
	for(j=0; j<m; j++)
		nnarcs+= narcs[j];
	sprintf(sss, "narcs = %d\n", nnarcs);
	printf(sss);
	sprintf(sss, "m = %6d nfix1 = %6d nfix0 = %6d unfixed = %6d\n", m, nfix1, nfix0, m - nfix0 - nfix1);
	print(sss);

	int **ad = new int*[m];
	for(i=0; i<m; i++)
	{
		ad[i] = new int[m];
		for(j=0; j<m; j++)
		{
			ad[i][j] = 0;
		}
	}

	for(j=0; j<m; j++)
	{
		for(i=0; i<narcs[j]; i++)
		{
			ad[arcs[j][i].i][j] = 1;
		}
	}
	for(i=0; i<m; i++)
	{
		for(j=0; j<m; j++)
		{
			if(i == j || !ad[i][j])
				continue;
			for(h=0; h<m; h++)
			{
				if(h==i || h==j || !ad[j][h])
					continue;
				ad[i][h] = 1;
			}
		}
	}
	int hh = 0;
	for(i=0; i<m; i++)
		for(j=0; j<m; j++)
			hh += ad[i][j];
	printf("hh = %d\n", hh);

	int ** fx = new int *[m];
	for(j=0; j<m; j++)
	{
		fx[j] = new int [narcs[j]];
		for(i=0; i<narcs[j]; i++)
		{
			fx[j][i] = -1;
		}
	}

	int nnf0 = 0;
	int nnf1 = 0;
	int naf0 = 0;
	int naf1 = 0;

	// node reduced cost fixing
	int nfx1 = 0;
	for(i=nfix1; i<p; i++)
	{
		if(blb + brho[p].d - brho[i].d < bub + EPSILON)
			break;
		brho[i].fix = 1;
		nfix1 ++;
		nnf1 ++;
	}
	if(nfix1 == p)
	{
		for(i=p; i<m-nfix0; i++)
		{
			brho[i].fix = 0;
			nfix0 ++;
			nnf0 ++;
		}
	}
	sprintf(sss, "RCF1: nnf1 = %6d nnf0 = %6d naf1= %6d naf0 = %6d\n", nnf1, nnf0, naf1, naf0);
	print(sss);
	for(i=m-nfix0-1; i>=p; i--)
	{
		if(blb + brho[i].d - brho[p-1].d < bub + EPSILON)
			break;
		brho[i].fix = 0;
		nfix0 ++;
		nnf0 ++;
	}
	sprintf(sss, "RCF2: nnf1 = %6d nnf0 = %6d naf1= %6d naf0 = %6d\n", nnf1, nnf0, naf1, naf0);
	print(sss);
	// arc reduced cost fixing
	for(j=0; j<nfix1; j++)
	{
		for(i=0; i<narcs[brho[j].i]; i++)
		{
			if(brho[fixed[arcs[brho[j].i][i].i]].fix !=1 || fx[brho[j].i][i] != -1)
				continue;
			if(blb - arcs[brho[j].i][i].d + blm[brho[j].i] > bub)
			{
				naf1 ++;
				fx[brho[j].i][i] = 1;
			}
		}
	}

	sprintf(sss, "RCF3: nnf1 = %6d nnf0 = %6d naf1= %6d naf0 = %6d\n", nnf1, nnf0, naf1, naf0);
	print(sss);

	for(j=0; j<p; j++)
	{
		for(i=0; i<narcs[brho[j].i]; i++)
		{
			if(fx[brho[j].i][i] != -1)
				continue;
			if(blb + arcs[brho[j].i][i].d - blm[brho[j].i] > bub)
			{
				naf0 ++;
				fx[brho[j].i][i] = 0;
			}
		}
	}

	sprintf(sss, "RCF4: nnf1 = %6d nnf0 = %6d naf1= %6d naf0 = %6d\n", nnf1, nnf0, naf1, naf0);
	print(sss);
	
	int *nodf = new int[m];
	for(i=0; i<m; i++)
	{
		nodf[brho[i].i] = brho[i].fix;
	}

	// logical fixing
	while(1)
	{
		//for(j=0; j<m; j++)
		//{
		//	if(nodf[j] != 0 || narcs[j] != 1)
		//		continue;
		//	if(nodf[arcs[j][0].i] != 1)
		//	{
		//		nodf[arcs[j][0].i] = 1;
		//		nnf1 ++;
		//	}
		//	if(fx[j][0] != 1)
		//	{
		//		fx[j][0] = 1;
		//		naf1 ++;
		//	}
		//}
		//sprintf(sss, "LIF0: nnf1 = %6d nnf0 = %6d naf1= %6d naf0 = %6d\n", nnf1, nnf0, naf1, naf0);
		//print(sss);

	
		for(j=0; j<m; j++)
		{
			for(i=0; i<narcs[j]; i++)
			{
				if(fx[j][i] != -1)
					continue;
				if(nodf[arcs[j][i].i] == 0)
				{
					fx[j][i] = 0;
					naf0 ++;
				}
			}
		}
		sprintf(sss, "LIF1: nnf1 = %6d nnf0 = %6d naf1= %6d naf0 = %6d\n", nnf1, nnf0, naf1, naf0);
		print(sss);

		for(j=0; j<m; j++)
		{
			if(nodf[j] != 1)
				continue;
			for(i=0; i<narcs[j]; i++)
			{
				if(fx[j][i] != -1)
					continue;
				fx[j][i] = 0;
				naf0 ++;
			}
		}
		sprintf(sss, "LIF2: nnf1 = %6d nnf0 = %6d naf1= %6d naf0 = %6d\n", nnf1, nnf0, naf1, naf0);
		print(sss);

		for(j=0; j<m; j++)
		{
			for(i=0; i<narcs[j]; i++)
			{
				if(fx[j][i] != 1)
					continue;
				if(nodf[arcs[j][i].i] != 1)
				{
					nodf[arcs[j][i].i] = 1;
					nnf1 ++;
				}
				if(nodf[j] != 0)
				{
					nodf[j] = 0;
					nnf0 ++;
				}
				for(h = 0; h<i; h++)
				{
					if(fx[j][h] == -1)
					{
						fx[j][h] = 0;
						naf0 ++;
					}
					if(nodf[arcs[j][h].i] == -1)
					{
						nodf[arcs[j][h].i] = 0;
						nnf0 ++;
					}
				}
				for(h = i+1; h<narcs[j]; h++)
				{
					if(fx[j][h] == -1)
					{
						fx[j][h] = 0;
						naf0 ++;
					}
				}
				break;
			}
		}

		sprintf(sss, "LIF3: nnf1 = %6d nnf0 = %6d naf1= %6d naf0 = %6d\n", nnf1, nnf0, naf1, naf0);
		print(sss);
	
		for(j=0; j<m; j++)
		{
			if(nodf[j] != -1)
				continue;
			h = 0;
			for(i=0; i<narcs[j]; i++)
				h += fx[j][i] != 0;
			if(h == 0)
			{
				nodf[j] = 1;
				nnf0 ++;
			}
		}
		sprintf(sss, "LIF4: nnf1 = %6d nnf0 = %6d naf1= %6d naf0 = %6d\n", nnf1, nnf0, naf1, naf0);
		print(sss);

		for(j=0; j<m; j++)
		{
			if(nodf[j] != 0)
				continue;
			h = 0;
			for(i=0; i<narcs[j]; i++)
			{
				h += fx[j][i] == -1;
			}
			if(h != 1)
				continue;
			h = 0;
			for(i=0; i<narcs[j]; i++)
			{
				h += fx[j][i] == 1;
			}
			if(h != 0)
				continue;
			for(i=0; i<narcs[j]; i++)
			{
				if(fx[j][i] != -1)
					continue;
				fx[j][i] = 1;
				naf1 ++;
				break;
			}
		}
		sprintf(sss, "LIF5: nnf1 = %6d nnf0 = %6d naf1= %6d naf0 = %6d\n", nnf1, nnf0, naf1, naf0);
		print(sss);

		for(j=0; j<m; j++)
		{
			if(nodf[j] == 1)
				continue;
			for(i=0; i<narcs[j]; i++)
			{
				if(nodf[arcs[j][i].i] != 1)
					continue;
				for(h=i+1; h<narcs[j]; h++)
				{
					if(fx[j][h] != 0)
					{
						fx[j][h] = 0;
						naf0 ++;
					}
				}
				break;
			}
		}

		sprintf(sss, "LIF6: nnf1 = %6d nnf0 = %6d naf1= %6d naf0 = %6d\n", nnf1, nnf0, naf1, naf0);
		print(sss);

		for(j=0; j<m; j++)
		{
			if(nodf[j] == 1)
				continue;
			for(i=m-p; i<narcs[j]; i++)
			{
				if(fx[j][i] != 0)
				{
					fx[j][i] = 0;
					naf0 ++;
				}
			}
		}
		sprintf(sss, "LIF7: nnf1 = %6d nnf0 = %6d naf1= %6d naf0 = %6d\n", nnf1, nnf0, naf1, naf0);
		print(sss);

		int k;
		for(j=0; j<m; j++)
		{
			for(i=0; i<narcs[j]; i++)
			{
				if(fx[j][i] != 1)
					continue;
				 l = arcs[j][i].i;
				 for(h=0; h<m; h++)
				 {
					 if(!ad[l][h] || !ad[h][j])
						 continue;
					 for(k = 0; k<narcs[h]; k++)
					 {
						 if(arcs[h][k].i != l || fx[h][k] == 1)
							 continue;
						 fx[h][k] = 1;
						 naf1 ++;
					 }
				 }
			}
		}

		for(j=0; j<m; j++)
		{
			for(i=0; i<narcs[j]; i++)
			{
				if(fx[j][i] != 0)
					continue;
				 l = arcs[j][i].i;
				 for(h=0; h<m; h++)
				 {
					 if(!ad[l][h] || !ad[j][h])
						 continue;
					 for(k = 0; k<narcs[h]; k++)
					 {
						 if(arcs[h][k].i != l || fx[h][k] == 0)
							 continue;
						 fx[h][k] = 0;
						 naf0 ++;
					 }
				 }
			}
		}

		sprintf(sss, "LIF8: nnf1 = %6d nnf0 = %6d naf1= %6d naf0 = %6d\n", nnf1, nnf0, naf1, naf0);
		print(sss);

		for(j=0; j<m; j++)
		{
			if(nodf[j] != 0)
				continue;;
			for(i=0; i<narcs[j]; i++)
			{
				if(fx[j][i] != 0)
					break;
				if(nodf[arcs[j][i].i] != 0)
				{
					nodf[arcs[j][i].i] = 0;
					nnf0 ++;
				}
			}
		}
		sprintf(sss, "LIF9: nnf1 = %6d nnf0 = %6d naf1= %6d naf0 = %6d\n", nnf1, nnf0, naf1, naf0);
		print(sss);

		for(j=0; j<m; j++)
		{
			if(nodf[j] != 0)
				break;
			for(i=0; i<narcs[j]; i++)
			{
				if(fx[j][i] != 0)
					break;
			}
			if(fx[j][i] == 1 || nodf[arcs[j][i].i] != 1)
				continue;
			fx[j][i] = 1;
			naf1 ++;
		}
		sprintf(sss, "LIF10: nnf1 = %6d nnf0 = %6d naf1= %6d naf0 = %6d\n", nnf1, nnf0, naf1, naf0);
		print(sss);

		if(nnf1 + nnf0 + naf1 + naf0 == 0)
			break;
		nnf1 = nnf0 = naf1 = naf0 = 0;
		print("--\n");
	}

	for(j=0; j<m; j++)
	{
		for(i=narcs[j] - 1; i>= 0; i--)
		{
			if(fx[j][i] != 0)
				continue;
			shift(arcs[j] + i, narcs[j] - i);
			narcs[j] --;
		}
	}
	nfixA1 = 0;
	constf = 0.;
	for(i=0; i<m; i++)
	{
		if(nodf[brho[i].i] == 0 && narcs[brho[i].i] == 1)
		{
			nfixA1 ++;
			brho[i].d = INFTY;
			constf += arcs[brho[i].i][0].d;
		}

		if(nodf[brho[i].i] == -1 || nodf[brho[i].i] == brho[i].fix)
			continue;
		if(nodf[brho[i].i] == 1)
		{
			brho[i].fix = 1;
			nfix1 ++;
			continue;
		}
		if(nodf[brho[i].i] == 0)
		{
			brho[i].fix = 0;
			nfix0 ++;
			continue;
		}
	}

	qsort(brho, m, sizeof(Arcr), compareArcrAZ2);

	memcpy(rho, brho, m * sizeof(Arcr));

	nnarcs = 0;
	for(j=0; j<m; j++)
		nnarcs+= narcs[j];
	sprintf(sss, "narcs = %d\n", nnarcs);
	printf(sss);
	sprintf(sss, "m = %6d nfix1 = %6d nfix0 = %6d unfixed = %6d\n", m, nfix1, nfix0, m - nfix0 - nfix1);
	print(sss);

	delete[] nodf;
	for(i=0; i<m; i++)
	{
		delete[] fx[i];
		delete[] ad[i];
	}
	delete[] ad;
	delete[] fx;
	printf("time = %6.2f\n", CPU_time() - ttime);
	writelp1("ll.lp");
}


void Pmed::bb(double ll, int ii, int pp, int *ss)
{
	int i;
	//printf("BB ii = %5d, pp = %5d, ll = %10.3f\n", ii, pp, ll);
	/*if(pp > m - nfix0)
	{
		printf("**\n");
		return;
	}*/
	prim(ss);
	if(pp+1 > m - nfix0)
	{
		//printf("**\n");
		return;
	}
	for(i=ii; i<pp; i++)
	{
		ll += (brho[pp].d - brho[i].d);
		if((bub - ll)/bub <= OPTGAP)
		{
			//printf("!! i = %5d ll = %10.3f\n", i, ll);
			continue;
		}
		ss[brho[i].i] = 0;
		ss[brho[pp].i] = 1;
		bb(ll, i+1, pp+1, ss);
		ss[brho[i].i] = 1;
		ss[brho[pp].i] = 0;
	}
}

void Pmed::bbLANG(int ii, int pp, int *ss)
{
	int i;
	/*printf("BB ii = %5d, pp = %5d, blb = %10.3f\n", ii, pp, blb);
	if(pp > m - nfix0)
	{
		printf("**\n");
		return;
	}*/
	prim(ss);
	if(pp+1 > m - nfix0)
	{
		//printf("**\n");
		return;
	}
	for(i=ii; i<pp; i++)
	{
		brho[i].fix = 0;
		memcpy(rho, brho, m*sizeof(Arcr));
		nfix0 ++;
		rho[i] = brho[m-nfix0];
		rho[m-nfix0] = brho[i];
		memcpy(lm, blm, m*sizeof(double));
		blb = -1.e30;
		par4();
		/*printf("brho:\n");
		int j;
		for(j=0; j<m; j++)
		{
			printf("%d", brho[j].fix);
		}
		printf("\n");
		printf("rho:\n");
		for(j=0; j<m; j++)
		{
			printf("%d", rho[j].fix);
		}
		printf("\n");*/
		langBB();
		if((bub - blb)/bub <= OPTGAP)
		{
			//printf("!! i = %5d ll = %10.3f\n", i, blb);
			brho[i].fix = -1;
			nfix0 --;
			continue;
		}
		ss[brho[i].i] = 0;
		ss[brho[pp].i] = 1;
		bbLANG(i+1, pp+1, ss);
		brho[i].fix = -1;
		nfix0 --;
		ss[brho[i].i] = 1;
		ss[brho[pp].i] = 0;
	}
}

void Pmed::readprob(int b_, const char *fn)
{
	strcpy(probname, fn);
	double tt = CPU_time();
	FILE * ff = NULL;
	ff = fopen(fn, "rb");
	printf("Reading file %s\n", fn);
	if(!ff)
		error("Cannot open file", fn);
	fread(&m, sizeof(int), 1, ff);
	n = 0;
	
	int i, j;
	double **dd = new double*[m];
	for(i=0; i<m; i++)
	{
		dd[i] = new double[m];
		fread(dd[i], sizeof(double), m, ff);
	}
	narcs = new int[m];
	arcs = new Arc*[m];
	c = new double[m];
	for(i=0; i<m; i++)
	{
		arcs[i] = new Arc[m];
		for(j=0; j<m; j++)
		{
			if(i!=j)
			{
				arcs[i][j].d = dd[j][i];
				arcs[i][j].i = j;
				continue;
			}
		}
		arcs[i][i].d = INFTY;
		arcs[i][i].i = i;
		qsort(arcs[i], m, sizeof(Arc), compareArcAZ);
		if(0)
		{
			printf("\n");
			for(j=0; j<100; j++)
			{
				printf("%d %f\n", arcs[i][j].i, arcs[i][j].d);
			}
		}
		narcs[i] = m;
		c[i] = dd[i][i];
	}
	for(i=0; i<m; i++)
		delete[] dd[i];
	delete[] dd;
	isfull = true;
	m1 = m-1;
	m2 = 1;

	bsol = new int[m];
	for(i=0; i<m; i++)
	{
		bsol[i] = 0;
		//c[i] = 0.;
	}
	
	printf("Done in %.2f seconds.\n", CPU_time() - tt);
	fclose(ff);
}

void Pmed::readprob(int b_, double **dd)
{
	double tt = CPU_time();
	
	m = b_;

	n = 0;
	
	int i, j;
	
	narcs = new int[m];
	arcs = new Arc*[m];
	c = new double[m];
	for(i=0; i<m; i++)
	{
		arcs[i] = new Arc[m];
		for(j=0; j<m; j++)
		{
			if(i!=j)
			{
				arcs[i][j].d = dd[j][i];
				arcs[i][j].i = j;
				continue;
			}
		}
		arcs[i][i].d = INFTY;
		arcs[i][i].i = i;
		qsort(arcs[i], m, sizeof(Arc), compareArcAZ);
		if(0)
		{
			printf("\n");
			for(j=0; j<100; j++)
			{
				printf("%d %f\n", arcs[i][j].i, arcs[i][j].d);
			}
		}
		narcs[i] = m;
		c[i] = dd[i][i];
	}
	
	//isfull = true;
	m1 = m-1;
	m2 = 1;

	bsol = new int[m];
	for(i=0; i<m; i++)
	{
		bsol[i] = 0;
		//c[i] = 0.;
	}
	
	printf("Done in %.2f seconds.\n", CPU_time() - tt);
}

void Pmed::readprob(const char *fn, bool firstnumber)
{
	strcpy(probname, fn);
	double tt = CPU_time();
	FILE * ff = NULL;
	ff = fopen(fn, "r");
	printf("Reading file %s\n", fn);
	if(!ff)
		error("Cannot open file", fn);
	m = 0;
	n = 0;
	int nn;
	int mm;
	fscanf(ff, "%ld", &mm);
	fscanf(ff, "%ld", &nn);
	printf("mm = %10d\nnn = %10d\n", mm, nn);
	m = mm;
	n = nn;
	printf("m = %10d\nn = %10d\n", m, n);
	coor = new double*[m];
	int i, j, k;
	for(i=0; i<m; i++)
	{
		coor[i] = new double[n];
	}
	for(i=0; i<m; i++)
	{
		k = i;
		if(firstnumber)
		{
			fscanf(ff, "%ld", &k);
			k--;
		}
		for(j=0; j<n; j++)
			fscanf(ff, "%lf", &coor[k][j]);
	}
	printf("Done in %.2f seconds.\n", CPU_time() - tt);
	fclose(ff);
	if(0)
	{
		for(i=0; i<m; i++)
		{
			printf("%5d", i+1);
			for(j=0; j<n; j++)
				printf(" %10.2f", coor[i][j]);
			printf("\n");
		}
	}
	bsol = new int[m];
	c = new double[m];
	for(i=0; i<m; i++)
	{
		bsol[i] = 0;
		c[i] = 0.;
	}
}

Pmed::Pmed(void):coor(NULL), m(0), n(0), arcs(NULL), narcs(NULL), numarcs(0), dispInterval(100),
				nlarcs(NULL), larcs(NULL), bsol(NULL), c(NULL), isfull(false), fixed(NULL), nfixA1(0), constf(0.)
{
}

Pmed::~Pmed(void)
{
	int i;
	if(coor)
	{
		for(i=0; i<m; i++)
			delete[] coor[i];
		delete[] coor;
	}
	if(larcs)
	{
		for(i=0; i<m; i++)
			delete[] larcs[i];
		delete[] larcs;
	}
	if(nlarcs)
		delete[] nlarcs;

	int k = p1 - p0 - 1;
	for(i=0; i<k; i++)
	{
		delete[] brhoa[i];
		delete[] blma[i];
		delete[] lbia[i];
		delete[] ubia[i];
	}
	
	delete[] nfix1a;
	delete[] nfix0a;
	delete[] brhoa;
	delete[] blma;
	delete[] lbia;
	delete[] ubia;
	delete[] buba;
	delete[] blba;

}
