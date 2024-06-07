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

#ifdef CPX
#include <cplex.h>
#endif

#ifdef XPRS
#include <xprs.h>
#endif

struct Pair
{
	int i,j;
};

/* int r = 3; */

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
	}

	char fn[200];
	sprintf(fn, "%s-%d.sol", probname, p);
	FILE *ff = fopen(fn, "w");

	fprintf(ff, "p = %d\n", p);
	fprintf(ff, "bub = %lf\n", bub);
	fprintf(ff, "Check = %d\n", bb);
	fprintf(ff, "Medians:\n");
	for(i=0; i<p; i++)
		fprintf(ff, "%d\n", medians[i]);
	fclose(ff);

	delete[] medians;
	
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

//void Pmed::writecore()
//{
//
//
//	printf("Writing core problem..\n");
//	double tt = CPU_time();
//
//	int nnodes = (int)(MAXNODE * (double)p);
//	if(nnodes > m)
//		nnodes = m;
//	int nnarcs = MAXCORE / nnodes;
//	if(nnarcs > m - 1)
//		nnarcs = m - 1;
//	printf("Core nodes = %d\nCore leaving arcs = %d\n", nnodes, nnarcs);
//
//	int *flag = new int[nnodes];
//	Arc ** core = new Arc*[nnodes];
//	int i, j, k;
//	
//
//	for(i=0; i<nnodes; i++)
//	{
//		flag[i] = brho[i].i;
//		core[i] = new Arc[nnarcs];
//		if(nnarcs <= narcs[flag[i]])
//		{
//			memcpy(core[i], arcs[flag[i]], nnarcs * sizeof(Arc));
//		}
//		else
//		{
//			long long pos = (long long) flag[i] * (long long) m * (long long) sizeof(Arc);
//			_fseeki64(fd, pos, SEEK_SET);
//			fread(core[i], sizeof(Arc), nnarcs, fd);
//		}
//	}
//
//	char fn[200] = "core.lp";
//	FILE *ff = fopen(fn, "w");
//
//	fprintf(ff, "min");
//	int nn = 0;
//	for(i=0; i<nnodes; i++)
//	{
//		for(j=0; j<nnarcs; j++)
//		{
//			if(nn)
//				fprintf(ff, " + ");
//			if(nn%5 == 0)
//				fprintf(ff, "\n");
//			fprintf(ff, "%lfx%d,%d", core[i][j].d, flag[i], core[i][j].i);
//			nn ++;
//		}
//	}
//	fprintf(ff, "\ns.t.");
//
//	for(k=0; k<m; k++)
//	{
//		nn = 0;
//		for(i=0; i<nnodes; i++)
//		{
//			for(j=0; j<nnarcs; j++)
//			{
//				if(core[i][j].i != k)
//					continue;
//				if(nn)
//					fprintf(ff, " + ");
//				if(nn%5 == 0)
//					fprintf(ff, "\n");
//				fprintf(ff, "x%d,%d", flag[i], core[i][j].i);
//				nn ++;
//			}
//		}
//		for(i=0; i<nnodes; i++)
//		{
//			if(flag[i] == k)
//				fprintf(ff, " + y%d", flag[i]);
//		}
//		fprintf(ff, " = 1\n");
//	}
//
//	for(i=0; i<nnodes; i++)
//	{
//		for(j=0; j<nnarcs; j++)
//		{
//			fprintf(ff, "\nx%d,%d - y%d <= 0", flag[i], core[i][j].i, flag[i]);
//		}
//	}
//
//	nn = 0;
//	for(i=0; i<nnodes; i++)
//	{
//		if(nn)
//			fprintf(ff, " + ");
//		if(nn%5 == 0)
//			fprintf(ff, "\n");
//		fprintf(ff, "y%d", flag[i]);
//		nn++;
//	}
//	fprintf(ff, " = %d\n", p);
//
//	nn = 0;
//	fprintf(ff, "binaries");
//	for(i=0; i<nnodes; i++)
//	{
//		if(nn%5 == 0)
//			fprintf(ff, "\n");
//		fprintf(ff, "y%d ", flag[i]);
//		nn++;
//	}
//	fprintf(ff, "\nend\n");
//
//	fclose(ff);
//	delete[] flag;
//	for(i=0; i<nnodes; i++)
//		delete[] core[i];
//	delete[] core;
//	printf("Done in %.2f seconds\n", CPU_time() - tt);
//
//	printf("Solving core (wait %.2f seconds)\n", CORETIME);
//	CPXENVptr env;
//	CPXLPptr  lp;
//	int status;
//	env = CPXopenCPLEX (&status);
//	if(status)
//	{
//		error("Cannot open Cplex\n");
//	}
//	lp = CPXcreateprob (env, &status, "FGL");
//	if(status)
//	{
//		error("Cannot create problem\n");
//	}
//
//	CPXsetdblparam(env, CPX_PARAM_TILIM, CORETIME);
//	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
//	if(status = CPXreadcopyprob(env, lp, "core.lp", 0))
//		error("CPXreadcopyprob");
//
//	
//
//	if(status = CPXmipopt(env, lp))
//		error("CPXmipopt(env, lp)");
//
//	double objval = 0.;
//
//	if(status = CPXgetmipobjval(env, lp, &objval))
//		error("CPXgetmipobjval(env, lp, &objval)");
//
//	printf("objval = %.2lf\n", objval);
//
//	if(objval < bub)
//		bub = objval;
//
//	if(status = CPXfreeprob (env, &lp))
//	{
//		error("Cannot free problem\n");
//	}
//	if(status = CPXcloseCPLEX (&env))
//	{
//		error("Cannot close Cplex\n");
//	}
//}
//
//void Pmed::writecore1()
//{
//
//
//	printf("Writing core problem..\n");
//	double tt = CPU_time();
//
//	int nnodes = (int)(MAXNODE * p);
//	if(nnodes > m)
//		nnodes = m;
//	int nnarcs = m-1;
//
//	if(MAXCORE > nnodes*nnarcs)
//		MAXCORE = nnodes*nnarcs;
//	printf("Core nodes = %d\nCore leaving arcs = %d\n", nnodes, MAXCORE);
//
//	int *flag = new int[nnodes];
//	Arc ** core = new Arc*[nnodes];
//	int i, j, k;
//	
//	double * dd = new double[nnodes*nnarcs];
//
//	for(i=0; i<nnodes; i++)
//	{
//		flag[i] = brho[i].i;
//		core[i] = new Arc[nnarcs];
//		if(nnarcs <= narcs[flag[i]])
//		{
//			memcpy(core[i], arcs[flag[i]], nnarcs * sizeof(Arc));
//		}
//		else
//		{
//			long long pos = (long long) flag[i] * (long long) m * (long long) sizeof(Arc);
//			_fseeki64(fd, pos, SEEK_SET);
//			fread(core[i], sizeof(Arc), nnarcs, fd);
//		}
//		for(j=0; j<nnarcs; j++)
//		{
//			dd[i * nnarcs + j] = core[i][j].d - blm[core[i][j].i];
//		}
//	}
//
//	qsort(dd, nnodes*nnarcs, sizeof(double), compareDoubleAZ);
//	double maxe = dd[MAXCORE-1];
//	printf("d0 = %f\ndm = %f\ndl = %f\n", dd[0], maxe, dd[nnodes*nnarcs-1]);
//	delete[] dd;
//
//
//
//	int totarcs = 0;
//	Arc *cc = new Arc[m-1];
//	int *ncore = new int[nnodes];
//	for(i=0;i<nnodes; i++)
//	{
//		memcpy(cc, core[i], nnarcs*sizeof(Arc));
//		ncore[i] = 0;
//		for(j=0; j<nnarcs; j++)
//		{
//			if(cc[j].d - blm[cc[j].i] <= maxe)
//				ncore[i] ++;
//		}
//		totarcs += ncore[i];
//		delete[] core[i];
//		core[i] = new Arc[ncore[i]];
//		ncore[i] = 0;
//		for(j=0; j<nnarcs; j++)
//		{
//			if(cc[j].d - blm[cc[j].i] <= maxe)
//			{
//				core[i][ncore[i]] = cc[j];
//				ncore[i] ++;
//			}
//		}
//	}
//	delete[] cc;
//	
//	printf("Core nodes = %d\nCore leaving arcs = %d\n", nnodes, totarcs);
//
//	char fn[200] = "core.lp";
//	FILE *ff = fopen(fn, "w");
//
//	fprintf(ff, "min");
//	int nn = 0;
//	for(i=0; i<nnodes; i++)
//	{
//		for(j=0; j<ncore[i]; j++)
//		{
//			if(nn)
//				fprintf(ff, " + ");
//			if(nn%5 == 0)
//				fprintf(ff, "\n");
//			fprintf(ff, "%lfx%d,%d", core[i][j].d, flag[i], core[i][j].i);
//			nn ++;
//		}
//	}
//	fprintf(ff, "\ns.t.\n");
//
//	for(k=0; k<m; k++)
//	{
//		nn = 0;
//		for(i=0; i<nnodes; i++)
//		{
//			for(j=0; j<ncore[i]; j++)
//			{
//				if(core[i][j].i != k)
//					continue;
//				if(nn)
//					fprintf(ff, " + ");
//				if(nn%5 == 0)
//					fprintf(ff, "\n");
//				fprintf(ff, "x%d,%d", flag[i], core[i][j].i);
//				nn ++;
//			}
//		}
//		for(i=0; i<nnodes; i++)
//		{
//			if(flag[i] == k)
//				fprintf(ff, " + y%d", flag[i]);
//		}
//		fprintf(ff, " = 1\n");
//	}
//
//	for(i=0; i<nnodes; i++)
//	{
//		for(j=0; j<ncore[i]; j++)
//		{
//			fprintf(ff, "\nx%d,%d - y%d <= 0", flag[i], core[i][j].i, flag[i]);
//		}
//	}
//
//	nn = 0;
//	for(i=0; i<nnodes; i++)
//	{
//		if(nn)
//			fprintf(ff, " + ");
//		if(nn%5 == 0)
//			fprintf(ff, "\n");
//		fprintf(ff, "y%d", flag[i]);
//		nn++;
//	}
//	fprintf(ff, " = %d\n", p);
//
//	nn = 0;
//	fprintf(ff, "binaries");
//	for(i=0; i<nnodes; i++)
//	{
//		if(nn%5 == 0)
//			fprintf(ff, "\n");
//		fprintf(ff, "y%d ", flag[i]);
//		nn++;
//	}
//	fprintf(ff, "\nend\n");
//
//	fclose(ff);
//	delete[] flag;
//	delete[] ncore;
//	for(i=0; i<nnodes; i++)
//		delete[] core[i];
//	delete[] core;
//	printf("Done in %.2f seconds\n", CPU_time() - tt);
//
//	printf("Solving core (wait %.2f seconds)\n", CORETIME);
//	CPXENVptr env;
//	CPXLPptr  lp;
//	int status;
//	env = CPXopenCPLEX (&status);
//	if(status)
//	{
//		error("Cannot open Cplex\n");
//	}
//	lp = CPXcreateprob (env, &status, "FGL");
//	if(status)
//	{
//		error("Cannot create problem\n");
//	}
//
//	CPXsetdblparam(env, CPX_PARAM_TILIM, CORETIME);
//	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
//	if(status = CPXreadcopyprob(env, lp, "core.lp", 0))
//		error("CPXreadcopyprob");
//
//	
//
//	if(status = CPXmipopt(env, lp))
//		error("CPXmipopt(env, lp)");
//
//	double objval = 0.;
//
//	if(status = CPXgetmipobjval(env, lp, &objval))
//		error("CPXgetmipobjval(env, lp, &objval)");
//
//	printf("objval = %.2lf\n", objval);
//
//	if(objval < bub)
//		bub = objval;
//
//	if(status = CPXfreeprob (env, &lp))
//	{
//		error("Cannot free problem\n");
//	}
//	if(status = CPXcloseCPLEX (&env))
//	{
//		error("Cannot close Cplex\n");
//	}
//}
//
//
void Pmed::writecore3()
{


	printf("Writing core problem..\n");
	double tt = CPU_time();

	int nnodes = (int)(MAXNODE * (double)p);
	if(nnodes > m)
		nnodes = m;
	int nnarcs = m;
	int i, j;

	if(MAXCORE > nnodes*nnarcs)
		MAXCORE = nnodes*nnarcs;
	printf("Core nodes = %d\nCore leaving arcs = %d\n", nnodes, MAXCORE);

	

	int *flag1 = new int[nnodes];
	int *flag2 = new int[m];
	Arc ** core = new Arc*[m];
	for(i=0; i<m; i++)
	{
		try
		{
			core[i] = new Arc[nnodes];
		}
		catch(...)
		{
			error("Out of memory 1");
		}
	}
	int nn;
	for(i=0; i<m; i++)
	{
		flag2[i] = 0;
	}
	for(i=0; i<nnodes; i++)
	{
		flag1[i] = brho[i].i;
		flag2[flag1[i]] = 1;
	}
	for(i=0; i<m; i++)
	{
		nn = 0;
		for(j=0; j<m; j++)
		{
			if(!flag2[arcs[i][j].i])
				continue;
			core[i][nn]=arcs[i][j];
			nn++;
		}
	}
	double *dd = NULL;
	try
	{
		dd = new double[nnodes*nnarcs];
	}
	catch(...)
	{
		error("Out of memory.");
	}
	for(i=0; i<m; i++)
	{
		for(j=0; j<nnodes; j++)
		{
			dd[i*nnodes + j] = core[i][j].d - blm[i];
		}
	}
	
	qsort(dd, nnodes*nnarcs, sizeof(double), compareDoubleAZ);

	double maxe = dd[MAXCORE-1];
	//printf("minrd = %f\nmaxe = %f\n, maxrd = %f\n", dd[0], maxe, dd[nnodes*nnarcs-1]);
	delete[] dd;

	int km = MAXCORE / (m-1);
	if(km > nnodes)
		km = nnodes;
	//printf("\n6\n");

	int *ncore = new int[m];
	int nc = 0;
	for(i=0; i<m; i++)
	{
		ncore[i] = 0;
		for(j=0; j<nnodes; j++)
		{
			if(core[i][j].d - blm[i] <= maxe)
				ncore[i] ++;
		}
		if(ncore[i] < 2)
			ncore[i] = 2;
		nc += ncore[i];
	}

	

	printf("Core nodes = %d\nCore leaving arcs = %d\n", nnodes, nc);

	char fn[200] = "core.lp";
	FILE *ff = fopen(fn, "w");

	nn = 0;	
	fprintf(ff, "min");
	for(i=0; i<nnodes; i++)
	{
		if(nn)
			fprintf(ff, " + ");
		if(nn%5 == 0)
			fprintf(ff, "\n");
		fprintf(ff, "%lfy%d", c[flag1[i]], flag1[i]);
		nn++;
	}
	for(i=0; i<m; i++)
	{
		for(j=0; j<ncore[i]; j++)
		{
			if(nn)
				fprintf(ff, " + ");
			if(nn%5 == 0)
				fprintf(ff, "\n");
			fprintf(ff, "%lfx%d,%d", core[i][j].d, core[i][j].i, i);
			nn ++;
		}
	}
	fprintf(ff, "\ns.t.\n");

	for(i=0; i<m; i++)
	{
		nn = 0;
		for(j=0; j<ncore[i]; j++)
		{
			if(nn)
				fprintf(ff, " + ");
			if(nn%5 == 0)
				fprintf(ff, "\n");
			fprintf(ff, "x%d,%d", core[i][j].i, i);
			nn ++;
		}
		for(j=0; j<nnodes; j++)
		{
			if(flag1[j] == i)
				fprintf(ff, " + y%d", i);
		}
		fprintf(ff, " = 1\n");
	}

	for(i=0; i<m; i++)
	{
		for(j=0; j<ncore[i]; j++)
		{
			fprintf(ff, "\nx%d,%d - y%d <= 0", core[i][j].i, i, core[i][j].i);
		}
	}

	nn = 0;
	for(i=0; i<nnodes; i++)
	{
		if(nn)
			fprintf(ff, " + ");
		if(nn%5 == 0)
			fprintf(ff, "\n");
		fprintf(ff, "y%d", flag1[i]);
		nn++;
	}
	fprintf(ff, " = %d\n", p);

	nn = 0;
	fprintf(ff, "binaries");
	for(i=0; i<nnodes; i++)
	{
		if(nn%5 == 0)
			fprintf(ff, "\n");
		fprintf(ff, "y%d ", flag1[i]);
		nn++;
	}
	fprintf(ff, "\nend\n");

	fclose(ff);
	for(i=0; i<m; i++)
		delete[] core[i];
	delete[] core;
	delete[] ncore;
	printf("Done in %.2f seconds\n", CPU_time() - tt);

	printf("Solving core (wait %.2f seconds)\n", CORETIME);
#ifdef CPX

	CPXENVptr env;
	CPXLPptr  lp;
	int status;
	env = CPXopenCPLEX (&status);
	if(status)
	{
		error("Cannot open Cplex\n");
	}
	lp = CPXcreateprob (env, &status, "FGL");
	if(status)
	{
		error("Cannot create problem\n");
	}

	CPXsetdblparam(env, CPX_PARAM_TILIM, CORETIME);
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
	if(status = CPXreadcopyprob(env, lp, "core.lp", 0))
		error("CPXreadcopyprob");

	

	if(status = CPXmipopt(env, lp))
		error("CPXmipopt(env, lp)");

	double objval = 0.;

	status = CPXgetmipobjval(env, lp, &objval);
	//	error("CPXgetmipobjval(env, lp, &objval)");

	printf("objval = %.2lf\n", objval);

	if(status == 0 /*&& objval < bub*/)
	{
		//bub = objval;
		double *y = new double[nnodes];
		if(status = CPXgetmipx(env, lp, y, 0, nnodes-1))
			error("CPXgetmipx(env, lp, y, 0, nnodes-1)");

		for(i=0; i<m; i++)
		{
			flag[i] = 0;
		}
		for(i=0; i<nnodes; i++)
		{
			if(y[i] < 0.5)
				continue;
			flag[flag1[i]] = 1;
		}

		primal();
		//savesol();

		delete[] y;
	}

	if(status = CPXfreeprob (env, &lp))
	{
		error("Cannot free problem\n");
	}
	if(status = CPXcloseCPLEX (&env))
	{
		error("Cannot close Cplex\n");
	}
#endif

#ifdef XPRS

	int status;
	XPRSprob xprob;

	//if(status = XPRSinit(0))
	//	error("status = XPRSinit(0)");

	if(status = XPRScreateprob(&xprob))
		error("XPRScreateprob(&xprob)");

	if(status = XPRSreadprob(xprob, "core.lp", ""))
		error("XPRSreadprob(xprob, core.lp)");
	
	if(status = XPRSsetintcontrol(xprob, XPRS_MAXTIME, (int)-CORETIME))
		error("status = XPRSsetdblcontrol(xprob, XPRS_MAXTIME, CORETIME)");

	if(status = XPRSminim(xprob, "g"))
		error("XPRSminim");

	double objval = 0.;

	status = XPRSgetdblattrib(xprob, XPRS_MIPOBJVAL, &objval);

	printf("objval = %.2lf\n", objval);

	if(status == 0)
	{
		
		int ncols;
		XPRSgetintattrib(xprob, XPRS_ORIGINALCOLS, &ncols);

		double *y = new double[ncols];
		if(status = XPRSgetmipsol(xprob, y, NULL))
			error("XPRSgetmipsol(xprob, y, NULL)");

		for(i=0; i<m; i++)
		{
			flag[i] = 0;
		}
		for(i=0; i<nnodes; i++)
		{
			if(y[i] < 0.5)
				continue;
			flag[flag1[i]] = 1;
		}

		primal();
		//savesol();

		delete[] y;
	}
	if(status = XPRSdestroyprob(xprob))
		error("XPRSdestroyprob(xprob)");
	//if(status = XPRSfree())
	//	error("XPRSfree()");

#endif

	delete[] flag1;
	delete[] flag2;
}
//
//void Pmed::writecore2()
//{
//
//	if(isfull)
//	{
//		writecore3();
//		return;
//	}
//
//	printf("Writing core problem..\n");
//	double tt = CPU_time();
//
//	int nnodes = (int)(MAXNODE * (double)p);
//	if(nnodes > m)
//		nnodes = m;
//	int nnarcs = m;
//	int i, j;
//
//	if(MAXCORE > nnodes*nnarcs)
//		MAXCORE = nnodes*nnarcs;
//	printf("Core nodes = %d\nCore leaving arcs = %d\n", nnodes, MAXCORE);
//
//	
//
//	int *flag1 = new int[nnodes];
//	Arc ** core = new Arc*[m];
//	for(i=0; i<m; i++)
//	{
//		try
//		{
//			core[i] = new Arc[nnodes];
//		}
//		catch(...)
//		{
//			error("Out of memory 1");
//		}
//	}
//	Arc *col = new Arc[m];
//	
//	for(i=0; i<nnodes; i++)
//	{
//		flag1[i] = brho[i].i;
//		if(nnarcs <= narcs[flag1[i]])
//		{
//			memcpy(col, arcs[flag1[i]], nnarcs * sizeof(Arc));
//		}
//		else
//		{
//			long long pos = (long long) flag1[i] * (long long) m * (long long) sizeof(Arc);
//			_fseeki64(fd, pos, SEEK_SET);
//			fread(col, sizeof(Arc), nnarcs, fd);
//		}
//		for(j=0; j<nnarcs; j++)
//		{
//			core[col[j].i][i].d = col[j].d;
//			core[col[j].i][i].i = flag1[i];
//		}
//	}
//	delete[] col;
//	double *dd = NULL;
//	try
//	{
//		dd = new double[nnodes*nnarcs];
//	}
//	catch(...)
//	{
//		error("Out of memory.");
//	}
//	printf("_ %d _", nnodes*nnarcs);
//	for(i=0; i<m; i++)
//	{
//		qsort(core[i], nnodes, sizeof(Arc), compareArcAZ);
//		for(j=0; j<nnodes; j++)
//		{
//			dd[i*nnodes + j] = core[i][j].d - blm[i];
//		}
//	}
//	printf("1");
//	qsort(dd, nnodes*nnarcs, sizeof(double), compareDoubleAZ);
//	printf("2");
//
//	double maxe = dd[MAXCORE-1];
//	printf("3");
//	//printf("minrd = %f\nmaxe = %f\n, maxrd = %f\n", dd[0], maxe, dd[nnodes*nnarcs-1]);
//	delete[] dd;
//	printf("4");
//
//	int km = MAXCORE / (m-1);
//	if(km > nnodes)
//		km = nnodes;
//	//printf("\n6\n");
//
//	printf("5");
//	int *ncore = new int[m];
//	int nc = 0;
//	for(i=0; i<m; i++)
//	{
//		ncore[i] = 0;
//		for(j=0; j<nnodes; j++)
//		{
//			if(core[i][j].d - blm[i] <= maxe)
//				ncore[i] ++;
//		}
//		if(ncore[i] < 2)
//			ncore[i] = 2;
//		nc += ncore[i];
//	}
//	printf("6");
//
//	
//
//	printf("Core nodes = %d\nCore leaving arcs = %d\n", nnodes, nc);
//
//	char fn[200] = "core.lp";
//	FILE *ff = fopen(fn, "w");
//
//	fprintf(ff, "min");
//	int nn = 0;
//	for(i=0; i<nnodes; i++)
//	{
//		if(nn)
//			fprintf(ff, " + ");
//		if(nn%5 == 0)
//			fprintf(ff, "\n");
//		fprintf(ff, "%lfy%d", c[flag1[i]], flag1[i]);
//		nn++;
//	}
//	for(i=0; i<m; i++)
//	{
//		for(j=0; j<ncore[i]; j++)
//		{
//			if(nn)
//				fprintf(ff, " + ");
//			if(nn%5 == 0)
//				fprintf(ff, "\n");
//			fprintf(ff, "%lfx%d,%d", core[i][j].d, core[i][j].i, i);
//			nn ++;
//		}
//	}
//	fprintf(ff, "\ns.t.\n");
//
//	for(i=0; i<m; i++)
//	{
//		nn = 0;
//		for(j=0; j<ncore[i]; j++)
//		{
//			if(nn)
//				fprintf(ff, " + ");
//			if(nn%5 == 0)
//				fprintf(ff, "\n");
//			fprintf(ff, "x%d,%d", core[i][j].i, i);
//			nn ++;
//		}
//		for(j=0; j<nnodes; j++)
//		{
//			if(flag1[j] == i)
//				fprintf(ff, " + y%d", i);
//		}
//		fprintf(ff, " = 1\n");
//	}
//
//	for(i=0; i<m; i++)
//	{
//		for(j=0; j<ncore[i]; j++)
//		{
//			fprintf(ff, "\nx%d,%d - y%d <= 0", core[i][j].i, i, core[i][j].i);
//		}
//	}
//
//	nn = 0;
//	for(i=0; i<nnodes; i++)
//	{
//		if(nn)
//			fprintf(ff, " + ");
//		if(nn%5 == 0)
//			fprintf(ff, "\n");
//		fprintf(ff, "y%d", flag1[i]);
//		nn++;
//	}
//	fprintf(ff, " = %d\n", p);
//
//	nn = 0;
//	fprintf(ff, "binaries");
//	for(i=0; i<nnodes; i++)
//	{
//		if(nn%5 == 0)
//			fprintf(ff, "\n");
//		fprintf(ff, "y%d ", flag1[i]);
//		nn++;
//	}
//	fprintf(ff, "\nend\n");
//
//	fclose(ff);
//	for(i=0; i<m; i++)
//		delete[] core[i];
//	delete[] core;
//	delete[] ncore;
//	printf("Done in %.2f seconds\n", CPU_time() - tt);
//
//	printf("Solving core (wait %.2f seconds)\n", CORETIME);
//	CPXENVptr env;
//	CPXLPptr  lp;
//	int status;
//	env = CPXopenCPLEX (&status);
//	if(status)
//	{
//		error("Cannot open Cplex\n");
//	}
//	lp = CPXcreateprob (env, &status, "FGL");
//	if(status)
//	{
//		error("Cannot create problem\n");
//	}
//
//	CPXsetdblparam(env, CPX_PARAM_TILIM, CORETIME);
//	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
//	if(status = CPXreadcopyprob(env, lp, "core.lp", 0))
//		error("CPXreadcopyprob");
//
//	
//
//	if(status = CPXmipopt(env, lp))
//		error("CPXmipopt(env, lp)");
//
//	double objval = 0.;
//
//	status = CPXgetmipobjval(env, lp, &objval);
//	//	error("CPXgetmipobjval(env, lp, &objval)");
//
//	printf("objval = %.2lf\n", objval);
//
//	if(status == 0 && objval < bub)
//	{
//		//bub = objval;
//		double *y = new double[nnodes];
//		if(status = CPXgetmipx(env, lp, y, 0, nnodes-1))
//			error("CPXgetmipx(env, lp, y, 0, nnodes-1)");
//
//		for(i=0; i<m; i++)
//		{
//			flag[i] = 0;
//		}
//		for(i=0; i<nnodes; i++)
//		{
//			if(y[i] < 0.5)
//				continue;
//			flag[flag1[i]] = 1;
//		}
//
//		primal();
//		//savesol();
//
//		delete[] y;
//	}
//
//	if(status = CPXfreeprob (env, &lp))
//	{
//		error("Cannot free problem\n");
//	}
//	if(status = CPXcloseCPLEX (&env))
//	{
//		error("Cannot close Cplex\n");
//	}
//
//	delete[] flag1;
//}

void Pmed::writecore4()
{

	if(isfull)
	{
		writecore3();
		return;
	}

	printf("Writing core problem..\n");
	double tt = CPU_time();

	int nnodes = (int)(MAXNODE * (double)p);
	if(nnodes > m)
		nnodes = m;
	int nnarcs = m;
	int i, j;

	if(MAXCORE > nnodes*nnarcs)
		MAXCORE = nnodes*nnarcs;
	printf("Core nodes = %d\nCore leaving arcs = %d\n", nnodes, MAXCORE);

	int *flag1 = new int[nnodes];
	int *flag2 = new int[m];
	int *ncore = new int[m];

	for(i=0; i<m; i++)
	{
		flag2[i] = 0;
		ncore[i] = 0;
	}
	for(i=0; i<nnodes; i++)
	{
		flag1[i] = brho[i].i;
		flag2[flag1[i]] = 1;
	}

	int nnc = 0;

	for(i=0; i<m; i++)
	{
		for(j=0; j<narcs[i]; j++)
		{
			int h = arcs[i][j].i;
			if(flag2[h] != 1)
				continue;
			if(h==i)
				continue;
			ncore[i] ++;
			if(ncore[i] == nnodes)
				break;
		}
		//if(ncore[i] <2)
		//	ncore[i] = 2;
		nnc += ncore[i];
	}

	printf("Core size = %d\n", nnc);

	double *dd = NULL;
	try
	{
		dd = new double[nnc];
	}
	catch(...)
	{
		error("Out of memory.");
	}
	int cnt = 0;
	for(i=0; i<m; i++)
	{
		ncore[i] = 0;
		for(j=0; j<narcs[i]; j++)
		{
			int h = arcs[i][j].i;
			if(flag2[h] != 1)
				continue;
			if(h==i)
				continue;
			dd[cnt] = arcs[i][j].d - blm[i];
			//printf(" III %d %d %f\n", arcs[i][j].i, i, arcs[i][j].d);
			cnt++;
			ncore[i]++;
			if(ncore[i] == nnodes)
				break;
		}
	}
	qsort(dd, cnt, sizeof(double), compareDoubleAZ);
	

	if(MAXCORE > cnt)
		MAXCORE = cnt;
	double maxe = dd[MAXCORE-1];
	printf("MAXCORE = %d\nmine = %f\nmaxe = %f\nmax = %f\n", MAXCORE, dd[0], maxe, dd[cnt-1]);
	delete[] dd;

	nnc = 0;
	for(i=0; i<m; i++)
	{
		ncore[i] = 0;
		for(j=0; j<narcs[i]; j++)
		{
			if(arcs[i][j].d - blm[i] > maxe)
				break;
			if(!flag2[arcs[i][j].i])
				continue;
			if(arcs[i][j].i==i)
				continue;
			ncore[i] ++;
		}
		nnc += ncore[i];
	}

	
	Arc ** core = new Arc*[m];
	for(i=0; i<m; i++)
	{
		try
		{
			if(ncore[i] < MINARCS)
				core[i] = new Arc[MINARCS];
			else
				core[i] = new Arc[ncore[i]];
		}
		catch(...)
		{
			error("Out of memory 1");
		}
	}
	nnc = 0;
	Arc *col = new Arc[m];
	int *order = new int[m];
	for(i=0; i<m; i++)
	{
		int nn = 0;
		if(ncore[i] < MINARCS)
		{
			int kk = narcs[i];
			for(j=0; j<narcs[i]; j++)
			{
				if(arcs[i][j].d - blm[i] > maxe)
				{
					kk = j;
					break;
				}
				if(!flag2[arcs[i][j].i])
					continue;
				if(arcs[i][j].i==i)
					continue;
				core[i][nn] = arcs[i][j];
				nn ++;
			}
			for(j=kk; j<narcs[i]; j++)
			{
				if(!flag2[arcs[i][j].i])
					continue;
				if(arcs[i][j].i==i)
					continue;
				core[i][nn] = arcs[i][j];
				nn ++;
				if(nn == MINARCS)
					break;
			}
			if(nn < MINARCS)
			{
				long long pos = (long long) i * (long long) m * (long long) sizeof(int);
				/* _fseeki64(fd, pos, SEEK_SET); */
				fread(order, sizeof(int), m, fd);
				for(j=narcs[i]; j<m; j++)
				{
					if(!flag2[order[j]])
						continue;
					if(order[j]==i)
						continue;
					core[i][nn].i = order[j];
					core[i][nn].d = getdist(i, order[j]);
					nn ++;
					if(nn == MINARCS)
						break;
				}
			}
		}
		else
		{
			for(j=0; j<narcs[i]; j++)
			{
				if(arcs[i][j].d - blm[i] > maxe)
					break;
				if(!flag2[arcs[i][j].i])
					continue;
				if(arcs[i][j].i==i)
					continue;
				core[i][nn] = arcs[i][j];
				nn ++;
			}
		}
		ncore[i] = nn;
		nnc += ncore[i];
	}
	delete[] col;
	delete[] order;
	//printf("Core\n");

	/*for(i=0; i<m; i++)
	{
		for(j=0; j<ncore[i]; j++)
		{
			printf("%d %d %f\n", core[i][j].i, i, core[i][j].d);
		}
	}*/

	printf("Core nodes = %d\nCore leaving arcs = %d\n", nnodes, nnc);

	char fn[200] = "core.lp";
	FILE *ff = fopen(fn, "w");

	fprintf(ff, "min");
	int nn = 0;
	for(i=0; i<nnodes; i++)
	{
		if(nn)
			fprintf(ff, " + ");
		if(nn%5 == 0)
			fprintf(ff, "\n");
		fprintf(ff, "%lfy%d", c[flag1[i]], flag1[i]);
		nn++;
	}
	for(i=0; i<m; i++)
	{
		for(j=0; j<ncore[i]; j++)
		{
			if(nn)
				fprintf(ff, " + ");
			if(nn%5 == 0)
				fprintf(ff, "\n");
			fprintf(ff, "%lfx%d,%d", core[i][j].d, core[i][j].i, i);
			nn ++;
		}
	}
	fprintf(ff, "\ns.t.\n");

	for(i=0; i<m; i++)
	{
		nn = 0;
		for(j=0; j<ncore[i]; j++)
		{
			if(nn)
				fprintf(ff, " + ");
			if(nn%5 == 0)
				fprintf(ff, "\n");
			fprintf(ff, "x%d,%d", core[i][j].i, i);
			nn ++;
		}
		for(j=0; j<nnodes; j++)
		{
			if(flag1[j] == i)
				fprintf(ff, " + y%d", i);
		}
		fprintf(ff, " = %d\n", r);
	}

	for(i=0; i<m; i++)
	{
		for(j=0; j<ncore[i]; j++)
		{
			fprintf(ff, "\nx%d,%d - y%d <= 0", core[i][j].i, i, core[i][j].i);
		}
	}

	nn = 0;
	for(i=0; i<nnodes; i++)
	{
		if(nn)
			fprintf(ff, " + ");
		if(nn%5 == 0)
			fprintf(ff, "\n");
		fprintf(ff, "y%d", flag1[i]);
		nn++;
	}
	fprintf(ff, " = %d\n", p);

	nn = 0;
	fprintf(ff, "binaries");
	for(i=0; i<nnodes; i++)
	{
		if(nn%5 == 0)
			fprintf(ff, "\n");
		fprintf(ff, "y%d ", flag1[i]);
		nn++;
	}
	fprintf(ff, "\nend\n");

	fclose(ff);
	for(i=0; i<m; i++)
		delete[] core[i];
	delete[] core;
	delete[] ncore;
	printf("Done in %.2f seconds\n", CPU_time() - tt);

	printf("Solving core (wait %.2f seconds)\n", CORETIME);

#ifdef CPX

	CPXENVptr env;
	CPXLPptr  lp;
	int status;
	env = CPXopenCPLEX (&status);
	if(status)
	{
		error("Cannot open Cplex\n");
	}
	lp = CPXcreateprob (env, &status, "FGL");
	if(status)
	{
		error("Cannot create problem\n");
	}

	CPXsetdblparam(env, CPX_PARAM_TILIM, CORETIME);
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);
	//CPXsetintparam(env, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_FEASIBILITY);
	if(status = CPXreadcopyprob(env, lp, "core.lp", 0))
		error("CPXreadcopyprob");

	

	if(status = CPXmipopt(env, lp))
		error("CPXmipopt(env, lp)");

	double objval = 0.;

	status = CPXgetmipobjval(env, lp, &objval);
	//	error("CPXgetmipobjval(env, lp, &objval)");

	printf("objval = %.2lf\n", objval);

	if(status == 0 /*&& objval < bub*/)
	{
		//bub = objval;
		double *y = new double[nnodes];
		if(status = CPXgetmipx(env, lp, y, 0, nnodes-1))
			error("CPXgetmipx(env, lp, y, 0, nnodes-1)");

		for(i=0; i<m; i++)
		{
			flag[i] = 0;
		}
		for(i=0; i<nnodes; i++)
		{
			if(y[i] < 0.5)
				continue;
			flag[flag1[i]] = 1;
		}

		primal();
		//savesol();

		delete[] y;
	}

	if(status = CPXfreeprob (env, &lp))
	{
		error("Cannot free problem\n");
	}
	if(status = CPXcloseCPLEX (&env))
	{
		error("Cannot close Cplex\n");
	}
#endif

#ifdef XPRS

	int status;
	XPRSprob xprob;

	//if(status = XPRSinit(0))
	//	error("status = XPRSinit(0)");

	if(status = XPRScreateprob(&xprob))
		error("XPRScreateprob(&xprob)");

	if(status = XPRSreadprob(xprob, "core.lp", ""))
		error("XPRSreadprob(xprob, core.lp)");
	
	if(status = XPRSsetintcontrol(xprob, XPRS_MAXTIME, (int)-CORETIME))
		error("status = XPRSsetdblcontrol(xprob, XPRS_MAXTIME, CORETIME)");

	if(status = XPRSminim(xprob, "g"))
		error("XPRSminim");

	double objval = 0.;

	status = XPRSgetdblattrib(xprob, XPRS_MIPOBJVAL, &objval);

	printf("objval = %.2lf\n", objval);
	int ncols;
		XPRSgetintattrib(xprob, XPRS_ORIGINALCOLS, &ncols);

	double *y = new double[ncols];
	status = XPRSgetmipsol(xprob, y, NULL);
	//error("XPRSgetmipsol(xprob, y, NULL)");
	if(status == 0)
	{
		
		for(i=0; i<m; i++)
		{
			flag[i] = 0;
		}
		for(i=0; i<nnodes; i++)
		{
			if(y[i] < 0.5)
				continue;
			flag[flag1[i]] = 1;
		}

		primal();
		//savesol();

	}
	delete[] y;
	if(status = XPRSdestroyprob(xprob))
		error("XPRSdestroyprob(xprob)");
	//if(status = XPRSfree())
	//	error("XPRSfree()");

#endif

	delete[] flag1;
	delete[] flag2;
}
void Pmed::writelp()
{
	char fn[200];
	double tt = CPU_time();
	printf("Writing LP..\n");
	sprintf(fn, "%s.lp", probname);
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
		fprintf(ff, "+ y%d = 1\n", i, r);
	}
	for(i=0; i<m; i++)
	{
		for(j=0; j<narcs[i]; j++)
		{
			fprintf(ff, "y%d - x%d,%d >= 0\n", arcs[i][j].i, arcs[i][j].i, i);
		}
	}
	nn = 0;
	for(i=0; i<m; i++)
	{
		if(nn)
			fprintf(ff, " + ");
		fprintf(ff, "y%d", i);
		nn ++;
		if(nn%5 == 0)
			fprintf(ff, "\n");
	}
	fprintf(ff, " = %d\n", p);
	fprintf(ff, "\nbinaries\n");
	nn = 0;
	for(i=0; i<m; i++)
	{
		fprintf(ff, " y%d", i);
		nn ++;
		if(nn%5 == 0) 
			fprintf(ff, "\n");
	}
	fprintf(ff, "\nend\n");
	fclose(ff);
	printf("Done in %.2f seconds\n", CPU_time() - tt);
}

void Pmed::primal()
{
	Arc *col = new Arc[m-1];
	int *sol = new int[m];
	int *order = new int[m];
	int k, h, i;
	double ub = 0;
	int nas = 0;
	for(k=0 ; k<m ; k++)
	{
		int kn = 0;
		if (flag[k] == 1) 
		{
			sol[k] = k;
			ub += c[k];
			kn ++;
			//continue;
		}

		if(kn == r)
			continue;

		bool bb;
		bb = false;
		for(h=0; h<narcs[k]; h++)
		{
			i = arcs[k][h].i;
			if(flag[i])
			{
				ub += arcs[k][h].d;
				sol[k] = i;
				nas ++;
				kn++;
				if(kn == r)
				{
					bb = true;
					break;
				}
								
			}
		}
		if(!bb)
		{
			long long pos = (long long) k * (long long) m * (long long) sizeof(int);
			/* _fseeki64(fd, pos, SEEK_SET); */
			fread(order, sizeof(int), m, fd);
			for(h=narcs[k]; h<m; h++)
			{
				i = order[h];
				if(flag[i])
				{
					ub += getdist(k, i);
					sol[k] = i;
					nas ++;
					break;
				}
			}
		}
	}

	/*double sd = 0.;
	for(k=0; k<m; k++)
	{
		sd += getdist(k, sol[k]);
		printf(" %d %d %lf\n", k, sol[k], getdist(k, sol[k]));
	}

	printf("sd = %f\n", sd);*/


	
	if(nas == m * r - p && ub < bub)
	{
		bub = ub;
		//for(i=0; i<p; i++)
		//	bestMedians[i] = rho[i].i;
		memcpy(bsol, sol, m*sizeof(int));
	}
	delete[] col;
	delete[] sol;
	delete[] order;
}

void Pmed::dualvalue()
{
	int i, j, h;
	for(i=0; i<m; i++)
	{
		flag[i] = false;
		rho[i].i = i;
		rho[i].d = c[i] - lm[i];
	}
	lb = .0;

	for(j=0; j<m ; j++)
	{
		for(h = 0; h<narcs[j]; h++)
		{
			i = arcs[j][h].i;
			if(arcs[j][h].d - lm[j] >= -EPSILON)
				break;
			rho[i].d += (arcs[j][h].d - lm[j]);
		}
		lb += r * lm[j];
	}

	qsort(rho, m, sizeof(Arc), compareArcAZ);

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
		memcpy(brho, rho, m*sizeof(Arc));
	}
}

void Pmed::subgradient()
{
	int i, j, k, h;
	for(k=0; k<m; k++)
	{
		j = rho[k].i;
		sg[j] = r - flag[j];
		for(h=0; h<narcs[j]; h++)
		{
			i = arcs[j][h].i;
			if(arcs[j][h].d - lm[j] >= -EPSILON) 
				break;
			if(flag[i] == 1) 
				sg[j]--;
		}
	}
	normSG = 0.0;
	for(i=0; i<m; i++)
		normSG += sg[i] * sg[i];
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
		if( ubi[i] > m - 1)
			continue;
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


void Pmed::greedy()
{
	printf("\nGreedy algorithm.\n");
	double tt = CPU_time();

	int *first = new int[m];
	int *second = new int[m];
	int i, k, j;
	for(i=0; i<m; i++)
	{
		first[i] = -1;
		second[i] = -1;
	}

	double ub = 0.;

	double dub = 0;
	double bestdub = 0.;

	int besti = -1;
	int fi;

	for(k=m-1; k>=p; k--)
	{
		besti = -1;
		bestdub = 1.e30;
		for(i=0; i<m; i++)
		{
			if(first[i] != -1)
				continue;
			fi = -1;
			for(j=0; j<narcs[i]; j++)
			{
				if(first[arcs[i][j].i] == -1)
				{
					fi = j;
					break;
				}
			}
			if(fi == -1)
				continue;
			dub = arcs[i][fi].d;
			if(nlarcs[i])
			{
				for(int l=0; l<nlarcs[i]; l++)
				{
					j = larcs[i][l].i;
					if(first[j] == -1)
						continue;
					if(arcs[j][first[j]].i == i)
					{
						if(second[j] == -1)
						{
							dub = bestdub +1.;
							break;
						}
						dub += arcs[j][second[j]].d - arcs[j][first[j]].d;
					}
				}
			}
			else
			{
				dub = bestdub +1.;
			}
			if(dub < bestdub)
			{
				bestdub = dub;
				besti = i;
			}

		}
		if(besti == -1)
		{
			//break;
			error("ku");
		}
		ub += bestdub;
		for(j=0; j<narcs[besti]; j++)
		{
			if(first[arcs[besti][j].i] == -1)
			{
				first[besti] = j;
				break;
			}
		}
		second[besti] = -1;
		for(j=j+1; j<narcs[besti]; j++)
		{
			if(first[arcs[besti][j].i] == -1)
			{
				second[besti] = j;
				break;
			}
		}
		//for(int l=0; l<narcs[besti]; l++)
		for(i=0; i<m; i++)
		{
			if(besti == i)
				continue;
			//i = arcs[besti][l].i;
			if(first[i] == -1)
				continue;
			if(arcs[i][first[i]].i != besti)
				continue;
			first[i] = second[i];
			if(second[i] == -1)
				continue;
			second[i] = -1;
			for(j=first[i]+1; j<narcs[i]; j++)
			{
				if(first[arcs[i][j].i] == -1)
				{
					second[i] = j;
					break;
				}
			}
			
		}
		//for(int l=0; l<narcs[besti]; l++)
		for(i=0; i<m; i++)
		{
			if(besti == i)
				continue;
			//i = arcs[besti][l].i;
			if(second[i] == -1)
				continue;
			if(besti != arcs[i][second[i]].i)
				continue;
			fi = second[i];
			second[i] = -1;
			for(j=fi+1; j<narcs[i]; j++)
			{
				if(first[arcs[i][j].i] == -1)
				{
					second[i] = j;
					break;
				}
			}
		}
		if(k%100 == 0)
			printf("it = %6d ub = %15.2f\n", k, ub);
	}

	printf("Cheking...\n");
	double cub = 0;
	int np = 0;
	for(i=0; i<m; i++)
	{
		if(first[i] == -1)
			continue;
		cub += arcs[i][first[i]].d;
		np ++;
	}
	printf("cub = %.2f np = %d n-p = %d\n", cub, np, m-p);

	printf("Cheking...\n");
	cub = 0;
	np = 0;
	for(i=0; i<m; i++)
	{
		if(first[i] == -1)
			continue;
		for(j=0; j<narcs[i]; j++)
		{
			if(first[arcs[i][j].i] == -1)
			{
				cub += arcs[i][j].d;
				np ++;
				if(first[i] != j)
				{
					printf("\ni = %d first = %d j = %d\n", i, first[i], j);
					printf("%f %f\n", arcs[i][first[i]].d, arcs[i][j].d);
				}

				break;
			}
		}
	}
	printf("cub = %.2f np = %d n-p = %d\n", cub, np, m-p);

	delete[] first;
	delete[] second;
	bub = cub;
	printf("Done in %.2f seconds\n", CPU_time() - tt);
}

void Pmed::greedy1()
{
	printf("\nGreedy algorithm.\n");
	double tt = CPU_time();

	int *first = new int[m];
	int *second = new int[m];
	int i, k, j;
	for(i=0; i<m; i++)
	{
		first[i] = -1;
		second[i] = -1;
	}

	double ub = 0.;

	double dub = 0;
	double bestdub = 0.;

	int besti = -1;
	int fi;

	for(k=m-1; k>=p; k--)
	{
		if(k == 691)
			printf("");
		besti = -1;
		bestdub = 1.e30;
		for(i=0; i<m; i++)
		{
			if(first[i] != -1)
				continue;
			fi = -1;
			for(j=0; j<narcs[i]; j++)
			{
				if(first[arcs[i][j].i] == -1)
				{
					fi = j;
					break;
				}
			}
			if(fi == -1)
				continue;
			dub = arcs[i][fi].d;
			for(j=0; j<m; j++)
			{
				if(first[j] == -1)
					continue;
				if(arcs[j][first[j]].i == i)
				{
					if(second[j] == -1)
					{
						dub = bestdub +1.;
						break;
					}
					dub += arcs[j][second[j]].d - arcs[j][first[j]].d;
				}
			}
			if(dub < bestdub)
			{
				bestdub = dub;
				besti = i;
			}

		}
		if(besti == -1)
			break;
		ub += bestdub;
		for(j=0; j<narcs[besti]; j++)
		{
			if(first[arcs[besti][j].i] == -1)
			{
				first[besti] = j;
				break;
			}
		}
		second[besti] = -1;
		for(j=j+1; j<narcs[besti]; j++)
		{
			if(first[arcs[besti][j].i] == -1)
			{
				second[besti] = j;
				break;
			}
		}
		for(i=0; i<m; i++)
		{
			if(first[i] == -1)
				continue;
			if(arcs[i][first[i]].i != besti)
				continue;
			first[i] = second[i];
			if(second[i] == -1)
				continue;
			second[i] = -1;
			for(j=first[i]+1; j<narcs[i]; j++)
			{
				if(first[arcs[i][j].i] == -1)
				{
					second[i] = j;
					break;
				}
			}
			
		}
		for(i=0; i<m; i++)
		{
			if(second[i] == -1)
				continue;
			if(besti != arcs[i][second[i]].i)
				continue;
			fi = second[i];
			second[i] = -1;
			for(j=fi+1; j<narcs[i]; j++)
			{
				if(first[arcs[i][j].i] == -1)
				{
					second[i] = j;
					break;
				}
			}
		}
		//printf("it = %6d ub = %15.2f\n", k, ub);
	}

	printf("Cheking...\n");
	double cub = 0;
	int np = 0;
	for(i=0; i<m; i++)
	{
		if(first[i] == -1)
			continue;
		cub += arcs[i][first[i]].d;
		np ++;
	}
	printf("cub = %.2f np = %d n-p = %d\n", cub, np, m-p);

	printf("Cheking...\n");
	cub = 0;
	np = 0;
	for(i=0; i<m; i++)
	{
		if(first[i] == -1)
			continue;
		for(j=0; j<narcs[i]; j++)
		{
			if(first[arcs[i][j].i] == -1)
			{
				cub += arcs[i][j].d;
				np ++;
				if(first[i] != j)
				{
					printf("\ni = %d first = %d j = %d\n", i, first[i], j);
					printf("%f %f\n", arcs[i][first[i]].d, arcs[i][j].d);
				}

				break;
			}
		}
	}
	printf("cub = %.2f np = %d n-p = %d\n", cub, np, m-p);

	delete[] first;
	delete[] second;
	printf("Done in %.2f seconds\n", CPU_time() - tt);
}

void Pmed::addarcs(int j)
{
	//printf("1");
	//printf("\nAdd arcs = %d was %d\n", j, narcs[j]);
	Arc *col = new Arc[narcs[j]];

	int ml = m2;
	if(ml + narcs[j] > m)
		ml = m - narcs[j];
	//printf("ml = %d\n", ml);


	memcpy(col, arcs[j], narcs[j] * sizeof(Arc));
	delete[] arcs[j];
	arcs[j] = new Arc[narcs[j] + ml];
	memcpy(arcs[j], col, narcs[j] * sizeof(Arc));

	int *order = new int[ml];
	long long pos = ((long long) j * (long long) m  + (long long) narcs[j]) * (long long) sizeof(int);
	/* _fseeki64(fd, pos, SEEK_SET); */
	fread(order, sizeof(int), ml, fd);
	int i;
	for(i=0; i<ml; i++)
	{
		arcs[j][narcs[j]+i].i = order[i];
		arcs[j][narcs[j]+i].d = getdist(order[i], j);
	}
	narcs[j] += ml;
	
	delete[] col;
	delete[] order;
	//printf("2\n");
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
		printf("%-20s%5d\n", "Number of medians per node:", r);
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
	rho = new Arc[m];
	brho = new Arc[m];
	
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

void Pmed::lang_init()
{
	if(!isfull && n)
	{
		printf("Open distance file\n");
		char fn[200];
		//sprintf(fn, "%s.dis", probname);
		sprintf(fn, "dis.dis");
		fd = fopen(fn, "rb");
		if(!fd)
			error("Cannot open file ", fn);
	}

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
		blm[i] = arcs[i][0].d;
		lmub[i] = arcs[i][narcs[i]-1].d;
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
	rho = new Arc[m];
	brho = new Arc[m];
}

void Pmed::lang_final()
{
	if(!isfull && n)
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

void Pmed::lang_set1()
{
	// set algorithm parameters

	if(!isfull)
	{
		m1 = (int)(MEMSIZE / (long long)sizeof(Arc) / (long long)m);
		if(m1 >= m)
			m1 = m - 1;
		m2 = 100;
	}

	phi = 2;
	stopphi = .005;
	phiFactor = 2;
	cfiter = 30;

	/*phi = 1.5;
	stopphi = .001;
	phiFactor = 1.01;
	cfiter = 1;*/

	EPSCORE = 1.e30;
	
	/*MAXNODE = 8;
	MAXCORE = 12*m;
	MINARCS = 3;*/
	//CORETIME = 60.;

	PRIMALFREQ = 2000000000;
}

void Pmed::lang_set2()
{
	// set algorithm parameters

	/*for(int i=0; i<m; i++)
	{
		lbi[i] = -1;
		ubi[i] = 1;
		lm[i] = arcs[i][0].d;
		blm[i] = arcs[i][0].d;
		//lmub[i] = arcs[i][narcs[i]-1].d;
	}*/

	/*phi = 2;
	stopphi = .005;
	phiFactor = 2;
	cfiter = 30;*/

	phi = 0.1;
	stopphi = .001;
	phiFactor = 1.01;
	cfiter = 1;

	EPSCORE = 1.e30;

	PRIMALFREQ = 0;
	
	
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

void Pmed::lang_set3()
{
	// set algorithm parameters

	phi = .005;
	stopphi = .001;
	phiFactor = 1.01;
	cfiter = 2;

	EPSCORE = 1.e30;

	PRIMALFREQ = 0;
	
	
}


void Pmed::lang()
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
		printf("%-20s%5d\n", "Number of medians per node:", r);
	}
	
	// start time
	double btime = CPU_time();

	// iteration counter
	int it = 0;
	fiter = 0;

	//double stopphi = 0.001;
	memcpy(lm, blm, m*sizeof(double));

	// main loop
	
	while(true)
	{	
		// evaluate lagrangean function
		dualvalue();
		if(it==0)
		{
			Eta = (bub - blb);
		}
		if (PRIMALFREQ && (it % PRIMALFREQ == 0))
			primal();
		if(dispInterval && it%dispInterval == 0)
		{
			printf("%6d %15.2f %12.2f %12.2f %12.2f %g\n", it, blb, bub, 
				(bub - blb) / bub * 100., CPU_time() - btime, phi);
		}
		// check optimality
		if( blb/bub >= 1. - DUALGAP
			)
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

	//if(PRIMALFREQ)
	//	primal();
	
	if(dispInterval)
	{
		printf("%6d %15.2f %12.2f %12.2f %12.2f\n", it, blb, bub, 
				(bub - blb) / bub * 100., CPU_time() - btime);
	}

	
}

//void Pmed::lagrangean1()
//{
//	if(dispInterval)
//	{
//		printf("*************************************************************\n");
//		printf("**                                                         **\n");
//		printf("**       Lagrangian heuristic for p-Median promlem         **\n");
//		printf("**       authors: Igor Vasiliev (vil@icc.ru)               **\n");
//		printf("**                                                         **\n");
//		printf("*************************************************************\n");
//		printf("%-20s%5d\n", "Number of nodes:", m);
//		printf("%-20s%5d\n", "Number of medians:", p);
//	}
//	printf("Open distance file\n");
//	char fn[200];
//	sprintf(fn, "%s.dis", probname);
//	fd = fopen(fn, "rb");
//	if(!fd)
//		error("Cannot open file ", fn);
//
//	int i;
//	// allocate and initialize lagrangean multipliers
//	lm = new double[m];
//	blm = new double[m];
//	lmub = new double[m];
//	lbi = new int[m];
//	ubi = new int[m];
//	for(i=0; i<m; i++)
//	{
//		lbi[i] = -1;
//		ubi[i] = 1;
//		lm[i] = arcs[i][0].d;
//		lmub[i] = arcs[i][narcs[i]-1].d;
//	}
//	// allocate subgradiend
//	sg = new double[m];
//
//	// start time
//	double btime = CPU_time();
//
//	// initialize bounds
//	blb = -1.e30;
//	// initialize best medians
//	flag = new int[m];
//	// initialize reduced costs
//	rho = new Arc[m];
//	brho = new Arc[m];
//	
//	
//
//	// set algorithm parameters
//	phi = 2.;
//	phiFactor = 1.1;
//	cfiter = 1;
//
//	// iteration counter
//	int it = 0;
//	fiter = 0;
//
//	double stopphi = 0.001;
//
//	// main loop
//	for(int g=0; g<2; g++)
//	{
//		switch(g)
//		{
//		case 0:
//			{
//				//printf("*");
//				phi = 2.;
//				stopphi = .001;
//				phiFactor = 1.01;
//				cfiter = 1;
//				break;
//			}
//		case 1:
//			{
//				phi = .001;
//				stopphi = .0001;
//				phiFactor = 1.01;
//				cfiter = 2;
//				break;
//			}
//		}
//		while(true)
//		{	
//			// evaluate lagrangean function
//			dualvalue();
//			if(it==0)
//			{
//				Eta = (bub - blb);
//			}
//			//primal();
//			if(dispInterval && it%dispInterval == 0)
//			{
//				printf("%6d %15.2f %12.2f %12.2f %12.2f %g\n", it, blb, bub, 
//					(bub - blb) / bub * 100., CPU_time() - btime, phi);
//			}
//			// check optimality
//			if( blb/bub >= 1. - DUALGAP
//				)
//			{
//				if(dispInterval)
//					printf("UB=LB\n");
//				break;
//			}
//			subgradient();
//			// check optimality
//			if(normSG < EPSILON)
//			{
//				if(dispInterval)
//						printf("nornSG = 0\n");
//				//bestUB = ub = lb;
//				break;
//			}
//			// go ahead subgradient
//			step1();
//			// check stop condition
//			if(phi < stopphi)
//			//if(theta < stopphi)
//			{
//				if(dispInterval)
//					printf("Stop by theta.\n");
//				break;
//			}
//			if(Eta < THETA_MIN)
//			{
//				if(dispInterval)
//					printf("Stop by theta.\n");
//				break;
//			}
//			it++;
//		}
//	}
//	if(dispInterval)
//	{
//		printf("%6d %15.2f %12.2f %12.2f %12.2f\n", it, blb, bub, 
//				(bub - blb) / bub * 100., CPU_time() - btime);
//	}
//
//	writecore();
//
//	fclose(fd);
//	//free memory
//	delete[] lm;
//	delete[] blm;
//	delete[] lmub;
//	delete[] sg;
//	delete[] rho;
//	delete[] brho;
//	delete[] flag;
//	delete[] lbi;
//	delete[] ubi;
//}

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
	printf("Reading distances\n");
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
		long long pos = (long long) i * (long long) m * (long long) sizeof(Arc);
		/* _fseeki64(ff, pos, SEEK_SET); */
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
}

void Pmed::writedist()
{
	if(!n)
		return;
	printf("Writing distances\n");
	char fn[200];
	//sprintf(fn, "%s.dis", probname);
	sprintf(fn, "dis.dis");
	FILE *ff = fopen(fn, "wb");
	double tt = CPU_time();


	int i, j;

	Arc *col = new Arc[m];
	int *order = new int[m];
	arcs = new Arc*[m];
	narcs = new int[m];
	
	double ub = 0;

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
		for(j=0; j<m; j++)
		{
			order[j] = col[j].i;
		}
		fwrite(order, sizeof(int), m, ff);
		narcs[i] = m1;
		arcs[i] = new Arc[m1];
		memcpy(arcs[i], col, m1 * sizeof(Arc));

		//printf("\n %d %d %f %f %f\n", m, m1, col[0].d, col[1].d, col[m1-1].d);
		
	}
	
	delete[] col;
	delete[] order;
	fclose(ff);
	printf("Done in %.2f seconds\n", CPU_time() - tt);
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

void Pmed::readprobmatr(char *fn)
{
	strcpy(probname, fn);
	double tt = CPU_time();
	
	FILE *ff = fopen(fn, "r");

	fscanf(ff, "%d", &m);
  printf("m=%d",m);
	n = 0;

	
	int i, j;
	double ** dd = new double*[m];
	for(i=0; i<m; i++)
	{
		dd[i] = new double[m];
		for(j=0; j<m; j++)
		{
			fscanf(ff, "%lf", &dd[i][j]);
      //printf("%d\n",dd[i][j]);
		}
		//dd[i][i] = INFTY;
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
	
	//isfull = true;
	m1 = m-1;
	m2 = 1;

	bsol = new int[m];
	for(i=0; i<m; i++)
	{
		bsol[i] = 0;
		//c[i] = 0.;
	}
	
	for(i=0; i<m; i++)
		delete[] dd[i];
	delete[] dd;
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
	fscanf(ff, "%d%d", &m, &n);
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
				nlarcs(NULL), larcs(NULL), bsol(NULL), c(NULL), isfull(false)
{
#ifdef XPRS
	int status;
	if(status = XPRSinit(0))
		error("status = XPRSinit(0)");
#endif
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
	if(arcs)
	{
		for(i=0; i<m; i++)
			delete[] arcs[i];
		delete[] arcs;
	}
	if(narcs)
		delete[] narcs;
	if(larcs)
	{
		for(i=0; i<m; i++)
			delete[] larcs[i];
		delete[] larcs;
	}
	if(nlarcs)
		delete[] nlarcs;
	if(bsol)
		delete[] bsol;
	if(isfull && n)
	{
		FILE *ff = fopen("dis.dis", "wb");
		fclose(ff);
	}
#ifdef XPRS
	int status;
	if(status = XPRSfree())
		error("status = XPRSfree(0)");
#endif
}
