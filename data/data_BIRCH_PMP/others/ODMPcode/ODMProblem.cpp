// ODMProblem.cpp: implementation of the ODMProblem class.
//
//////////////////////////////////////////////////////////////////////

#include "ODMProblem.h"
#include "ommom.h"

#ifdef _MSC_VER
using namespace std;
#endif
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////


extern int nprob1;
extern int nprob2;
extern int nprob3;

int compareGraphComponentsZA( const void *arg1, const void *arg2 )
{
   if( (* ( GraphComponent* ) arg1).nmodels < (* ( GraphComponent* ) arg2 ).nmodels )
	   return 1;
   else
	   return -1;
}

struct IT
{
	 int i, p;
	 int z;
};


int compareITZA( const void *arg1, const void *arg2 )
{
	if(((IT*) arg1)->z < ((IT*) arg2)->z)
		return 1;
	if(((IT*) arg1)->z > ((IT*) arg2)->z)
		return -1;
	if(((IT*) arg1)->p < ((IT*) arg2)->p)
		return 1;
	else
		return -1;
}


ODMProblem::ODMProblem(char* probname_, int nmodels_, __int64* modeloptions_, double* modelcosts_, 
					   double* modeldemands_, char** modelnames_, int noptions_, 
					   int* hardoptions_, 
					   int dp_, int** solutions_, double* values_, double constructioncost_):
					   probname(probname_), nmodels(nmodels_), modeloptions(modeloptions_), modelcosts(modelcosts_), 
					   modeldemands(modeldemands_), modelnames(modelnames_), noptions(noptions_), 
					   hardoptions(hardoptions_), 
					   dp(dp_), solutions(solutions_), values(values_),
					   constructioncost(constructioncost_)
						
{
	if(noptions>32)
		error("Too many options.");
	if(dp>100)
		error("dp is too big.");

	timer=new CPUTIMER;
	soltimer=new CPUTIMER;
	String fn;
	sprintf(fn, "%sODM.log", probname);
	setlogfile(fn);


	timer->start();

	print("\n-----\n");
	sprintf(fn, "%-30s", "Graph splitting..");
	print(fn);
	divide();
	timer->printTimer();
	print("\n");
	printst();
	solved=false;

	lb = new double[p1-p0+1];
	for(int p=p0; p<=p1; p++)
	{
		lb[p-p0] = 0.;
	}



	
}

void dive(int *comp, int u, int ind, Arc1** earcs, int* nearcs, Arc1 **larcs, int *nlarcs)
{
	int i;
	for(i=0; i<nearcs[u]; i++)
	{
		if(comp[earcs[u][i].u] != -1)
			continue;
		comp[earcs[u][i].u]= ind;
		dive(comp, earcs[u][i].u, ind, earcs, nearcs, larcs, nlarcs);
	}
	for(i=0; i<nlarcs[u]; i++)
	{
		if(comp[larcs[u][i].u] != -1)
			continue;
		comp[larcs[u][i].u]= ind;
		dive(comp, larcs[u][i].u, ind, earcs, nearcs, larcs, nlarcs);
	}
}

ODMProblem::ODMProblem(char* probname_, int p1_):
					   probname(probname_), nmodels(0), modeloptions(0), modelcosts(0), 
					   modeldemands(0), modelnames(0), noptions(0), 
					   hardoptions(0), 
					   constructioncost(0.), p0(0), p1(p1_),
					   modelhardoptions(0), modelcomponents(0)
						
{
	timer=new CPUTIMER;
	soltimer=new CPUTIMER;
	String fn;
	sprintf(fn, "%sODM.log", probname);
	setlogfile(fn);
	timer->start();

	print("\n-----\n");
	sprintf(fn, "%-30s", "Graph reading..");
	print(fn);
	FILE *ff = fopen(probname, "r");
	fscanf(ff, "%ld%ld", &nmodels, &narcs);
	Arc2 *arcs = new Arc2[narcs];

	int i, j;

	for(i=0; i<narcs; i++)
	{
		fscanf(ff, "%ld%ld%lf", &arcs[i].u, &arcs[i].v, &arcs[i].cost);
		arcs[i].u--;
		arcs[i].v--;
	}
	fclose(ff);
	timer->printTimer();
	print("\n");

	print("\n-----\n");
	sprintf(fn, "%-30s", "Graph splitting..");
	print(fn);
	nearcs = new int[nmodels];
	int *nlarcs = new int[nmodels];
	for(i=0; i<nmodels; i++)
	{
		nearcs[i] = nlarcs[i] = 0;
	}
	for(j=0; j<narcs; j++)
	{
		//if(arcs[j].cost > 0.)
		{
			nearcs[arcs[j].v] ++;
			nlarcs[arcs[j].u] ++;
		}
	}
	
	earcs = new Arc1*[nmodels];
	Arc1 **larcs = new Arc1*[nmodels];
	for(i=0; i<nmodels; i++)
	{
		earcs[i] = new Arc1[nearcs[i]];
		larcs[i] = new Arc1[nlarcs[i]];
		nearcs[i] = nlarcs[i] = 0;
	}
	for(j=0; j<narcs; j++)
	{
		//if(arcs[j].cost > 0.)
		{
			earcs[arcs[j].v][nearcs[arcs[j].v]].u = arcs[j].u;
			earcs[arcs[j].v][nearcs[arcs[j].v]].cost = arcs[j].cost;
			nearcs[arcs[j].v]++;
			larcs[arcs[j].u][nlarcs[arcs[j].u]].u = arcs[j].v;
			larcs[arcs[j].u][nlarcs[arcs[j].u]].cost = arcs[j].cost;
			nlarcs[arcs[j].u]++;
		}
	}
	for(i=0; i<nmodels; i++)
	{
		//qsort(earcs[i], nearcs[i], sizeof(Arc1), compareArcCostArrayAZ);
	}
	
	ngraphcomponents = 0;
	modelcomponents = new int[nmodels];
	for(i=0; i<nmodels; i++)
	{
		modelcomponents[i] = -1;
	}

	for(i=0; i<nmodels; i++)
	{
		if(modelcomponents[i] != -1)
			continue;
		modelcomponents[i] = ngraphcomponents;
		dive(modelcomponents, i, ngraphcomponents, earcs, nearcs, larcs, nlarcs);
		ngraphcomponents++;
	}

	//printf("ngraphcomponents = %d\n", ngraphcomponents);

	subgraphs=new GraphComponent*[ngraphcomponents];
	p0=0;
	narcs=0;
	toobig=false;
	for(i=0; i<ngraphcomponents; i++)
	{
		//printf("gf = %d\n", i);
		subgraphs[i]=new GraphComponent(this, i, earcs, nearcs);
		p0+= subgraphs[i]->p0;
		narcs+= subgraphs[i]->narcs;
		toobig = toobig || subgraphs[i]->narcs>MANYARCS;
	}

	if (p0 > p1)
	{
		char ss[250];
		sprintf(ss, "Number of roots %d more than the max munber of medians %d.", p0, p1);
		error(ss);
	}
	
	dp = p1 - p0;

	for(i=0; i<ngraphcomponents; i++)
	{
		subgraphs[i]->getp1();
	}
	//qsort(subgraphs, ngraphcomponents, sizeof(GraphComponent), compareGraphComponentsZA);
	timer->printTimer();
	print("\n");
	printst();
	solved=false;

	lb = new double[dp+1];
	values = new double[dp+1];
	solutions = new int*[dp+1];
	for(int p=p0; p<=p1; p++)
	{
		solutions[p-p0] = new int[nmodels];
		lb[p-p0] = 0.;
	}
	delete arcs;
	for(i=0; i<nmodels; i++)
	{
		//delete[] earcs[i];
		delete[] larcs[i];
	}
	//delete[] earcs;
	delete[] larcs;
	//delete[] nearcs;
	delete[] nlarcs;
}

ODMProblem::~ODMProblem()
{

	int i;
	print("-----\n");
	String s;
	sprintf(s, "%-30s", "Cleaning up..");
	printf(s);

	delete[] modelhardoptions;
	delete[] modelcomponents;

	for(i=0; i<ngraphcomponents; i++)
	{
		delete subgraphs[i];
	}
	delete[] subgraphs;

	delete[] lb;
	delete[] values;
	for(i=0; i<dp+1; i++)
	{
		delete[] solutions[i];
	}
	delete[] solutions;

	for(i=0; i<nmodels; i++)
	{
		delete[] earcs[i];
	}
	delete[] earcs;

	timer->printTimer();
	print("\n");
	
	print("-----\n");
	delete timer;
	delete soltimer;
}

int ODMProblem::divide()
{
	__int64 hopt=0;
	int i;
	for( i=0; i<noptions; i++)
	{
		if(hardoptions[i])
		{
			hopt+= (__int64) 1 << i;
		}
	}

	modelcomponents=new int[nmodels];
	modelhardoptions=new __int64[nmodels];
	for(i=0; i<nmodels; i++)
	{
		modelcomponents[i]=-1;
		modelhardoptions[i]=modeloptions[i] & hopt;
	}

	ngraphcomponents=0;

#ifndef _STUPID_
	int begmod=0;
	while(1)
	{
		bool empty=true;
		int i;
		for( i=begmod; i<nmodels; i++)
		{
			if(modelcomponents[i]==-1)
			{
				begmod=i;
				modelcomponents[i]=ngraphcomponents++;
				empty=false;
				break;
			}
		}
		if(empty)
			break;
		for(i=begmod+1; i<nmodels; i++)		
		{
			if(modelhardoptions[begmod]==modelhardoptions[i])
			{
				modelcomponents[i]=ngraphcomponents-1;
			}
		}
	}
#else
	for(i=0; i<nmodels; i++)
	{
		modelcomponents[i] = 0;
	}
	ngraphcomponents = 1;
#endif

	subgraphs=new GraphComponent*[ngraphcomponents];
	p0=0;
	narcs=0;
	toobig=false;
	for(i=0; i<ngraphcomponents; i++)
	{
		subgraphs[i]=new GraphComponent(this, i);
		p0+= subgraphs[i]->p0;
		narcs+= subgraphs[i]->narcs;
		toobig = toobig || subgraphs[i]->narcs>MANYARCS;
	}

	p1=MIN(nmodels, p0+dp);

	for(i=0; i<ngraphcomponents; i++)
	{
		subgraphs[i]->getp1();
	}

	

	return 0;

}

int ODMProblem::printst()
{

	String s;
	sprintf(s, "%-30s%15d\n", "Number of models", nmodels);
	print(s);
	sprintf(s, "%-30s%15d\n", "Number of arcs", narcs);
	print(s);
	sprintf(s, "%-30s%15d\n", "Number of roots", p0);
	print(s);
	sprintf(s, "%-30s%15d\n", "Max number of models", p1);
	print(s);
	sprintf(s, "%-30s%15d\n", "Number of components", ngraphcomponents);
	print(s);
	for(int i=0; i<ngraphcomponents; i++)
	{
		sprintf(s, "C%4.4d: %6d %8d %4d %4d\n", i, subgraphs[i]->nmodels,
			     subgraphs[i]->narcs, subgraphs[i]->p0, subgraphs[i]->p1);
		print(s);
	}

	return 0;
}

int ODMProblem::optCplex()
{

	/*soltimer->start(1);
	String s;
	sprintf(s, "%-30s\n", "Looking for optimal solution using decomposition..");
	print(s);
	if(toobig)
	{
		print("The problem is to big. Can do it.\n");
		return 0;
	}
	for(int i=0; i<ngraphcomponents; i++)
	{
		subgraphs[i]->optCplex();
	}

	sprintf(s, "Knapsack %4d-%4d: init", p0, p1);
	print(s);
	int status=0;
	env = CPXopenCPLEX (&status);
	if(!env)
		error("env = CPXopenCPLEX (&status);");
	// create problem
	lp=CPXcreateprob(env,&status,"ll");
	if(!lp)
		error("lp==CPXcreateprob(env,&status,k);");

	status = CPXchgprobtype (env, lp, 1);
	if(status)
		error("status = CPXchgprobtype (env, lp, CPXPROB_MIP);");

	//CPXsetintparam (env, CPX_PARAM_SCRIND, 1);


	
	
	createKnapsackLP();
	double* x=new double[nknapvars];

	for(int p=p0; p<=p1; p++)
	{
		sprintf(s, " f(%d)", p);
		print(s);
		int cnt=1;
		int indices[1]={0};
		double rhs[1]={p};


		status = CPXchgrhs (env, lp, cnt, indices, rhs);
		if(status)
		{
			error("status = CPXchgrhs (env, lp, cnt, indices, rhs);");
		}

		status = CPXmipopt(env, lp);
		if(status)
			error("status = CPXmipopt(env, lp);");

		status = CPXgetmipobjval (env, lp, values+p-p0);
		if(status)
			error("status = CPXgetmipobjval (env, lp, values+p-p0);");
		sprintf(s, "=%.2f", p*constructioncost+values[p-p0]);
		print(s);

		status = CPXgetmipx(env, lp, x, 0, nknapvars-1);
		if(status)
			error("status = CPXgetmipx(env, lp, x, 0, nknapvars-1);");


		for(int i=0; i<nmodels; i++)
		{
			solutions[p-p0][i]=i;
		}

		int cur=0;
		for(int k=0; k<ngraphcomponents; k++)
		{
			for(int pp=subgraphs[k]->p0; pp<=subgraphs[k]->p1; pp++)
			{
				if(x[cur++]<0.5)
					continue;
				for(int i=0; i<subgraphs[k]->nmodels; i++)
				{
					solutions[p-p0][subgraphs[k]->models[i]] = 
						subgraphs[k]->models[subgraphs[k]->solutions[pp-subgraphs[k]->p0][i]];
				}
			}
		}
	}


	status=0;
	status=CPXfreeprob(env,&lp);
	if(status)
		error("status=CPXfreeprob(env,&lp);");
	status = CPXcloseCPLEX (&env);
	if(status)
		error("status = CPXcloseCPLEX (&env);");
	
	

	print("\n");
	soltimer->stop();
	sprintf(s, "%-30s%15d\n", "Solution time (secs)", soltimer->seconds());
	print(s);
	sprintf(s, "%-30s", "Solution time ");
	print(s);
	soltimer->printTimer();
	print("\n");
	solved=true;
	delete[] x;
	*/
	return 0;
}

int ODMProblem::optCplexFull()
{
	/*soltimer->start(true);
	String s;
	sprintf(s, "%-30s\n", "Looking for optimal solution..");
	print(s);
	if(toobig)
	{
		print("The problem is to big. Can do it.\n");
		return 0;
	}
	sprintf(s, "%4d-%4d: init", p0, p1);
	print(s);
	if(narcs>50000)
	{
		print(" it is big. be patient.");
	}

	int status=0;
	env = CPXopenCPLEX (&status);
	if(!env)
		error("env = CPXopenCPLEX (&status);");
	// create problem
	lp=CPXcreateprob(env,&status,"ll");
	if(!lp)
		error("lp==CPXcreateprob(env,&status,k);");

	status = CPXchgprobtype (env, lp, 1);
	if(status)
		error("status = CPXchgprobtype (env, lp, CPXPROB_MIP);");

	//CPXsetintparam (env, CPX_PARAM_SCRIND, 1);

	createLP();
	double* x=new double[nmodels+narcs];

	for(int p=p0; p<=p1; p++)
	{
		sprintf(s, " f(%d)", p);
		print(s);
		int cnt=1;
		int indices[1]={0};
		double rhs[1]={p};


		status = CPXchgrhs (env, lp, cnt, indices, rhs);
		if(status)
		{
			error("status = CPXchgrhs (env, lp, cnt, indices, rhs);");
		}

		status = CPXmipopt(env, lp);
		if(status)
			error("status = CPXmipopt(env, lp);");

		status = CPXgetmipobjval (env, lp, values+p-p0);
		if(status)
			error("status = CPXgetmipobjval (env, lp, values+p-p0);");
		sprintf(s, "=%.2f", p*constructioncost+values[p-p0]);
		print(s);

		status = CPXgetmipx(env, lp, x, 0, nmodels+narcs-1);
		if(status)
			error("status = CPXgetmipx(env, lp, x, 0, nmodels+narcs-1);");


		for(int i=0; i<nmodels; i++)
		{
			solutions[p-p0][i]=i;
		}

		int cur=nmodels;
		for(int k=0; k<ngraphcomponents; k++)
		{
			for(int i=0; i<subgraphs[k]->nmodels; i++)
			{
				for(int j=0; j<subgraphs[k]->nearcs[i]; j++)
				{
					if(x[cur++]<0.5)
						continue;
					solutions[p-p0][subgraphs[k]->models[i]] = 
						subgraphs[k]->models[subgraphs[k]->earcs[i][j].u];
				}
			}
		}


		//CPXwriteprob(env, lp, "mm.lp", 0);

	}

	status=0;
	status=CPXfreeprob(env,&lp);
	if(status)
		error("status=CPXfreeprob(env,&lp);");
	status = CPXcloseCPLEX (&env);
	if(status)
		error("status = CPXcloseCPLEX (&env);");
	print("\n");

	soltimer->stop();
	sprintf(s, "%-30s%15d\n", "Solution time (secs)", soltimer->seconds());
	print(s);
	sprintf(s, "%-30s", "Solution time ");
	print(s);
	soltimer->printTimer();
	print("\n");
	solved=true;
	delete[] x;
	*/
	return 0;
}

int ODMProblem::createLP()
{
	/*int nvar=nmodels+narcs;

	int status=0;
	double* obj=new double[nvar];;
	double* lb=new double[nvar];;
	double* ub=new double[nvar];
	char** colname=new char*[nvar];
	int* indices=new int[nvar];
	char* ctype=new char[nvar];

	int cur=0;
	for(int k=0; k<ngraphcomponents; k++)
	{
		for(int i=0;i<subgraphs[k]->nmodels;i++)
		{
			obj[cur]=0.;
			lb[cur]=0.;
			ub[cur]=1.;
			colname[cur]=new char[50];
			sprintf(colname[cur],"y%d", subgraphs[k]->models[i]);
			//cout << subgraphs[k]->models[i] << endl;
			indices[cur]=subgraphs[k]->models[i];
			ctype[cur]='B';
			obj[cur]=0.;
			cur++;
		}
	}
	cur=nmodels;
	for(k=0; k<ngraphcomponents; k++)
	{
		for(int i=0;i<subgraphs[k]->nmodels;i++)
		{
			for(int j=0; j<subgraphs[k]->nearcs[i]; j++)
			{
				obj[cur]=0.;
				lb[cur]=0.;
				ub[cur]=1.;
				colname[cur]=new char[50];
				sprintf(colname[cur],"x%d,%d", 
					subgraphs[k]->models[subgraphs[k]->earcs[i][j].u], subgraphs[k]->models[i]);
				indices[cur]=cur;
				ctype[cur]='C';
				obj[cur]=subgraphs[k]->earcs[i][j].cost;
				cur++;
			}
		}
	}
	status=CPXaddcols (env, lp, cur, 0, obj, 0, 0, 0, lb, ub,  colname);
	if(status)
		error("status=CPXaddcols (env, lp, n, 0, obj, 0, 0, 0, lb, ub, colname);");

	status=CPXchgctype (env, lp, cur, indices, ctype);
	if(status)
		error("status=CPXchgctype (env, lp, nvar, indices, ctype);");


	int maxr=nmodels+narcs+1;
	int maxe=2*nmodels+3*narcs;
	int rcnt=0;
	int nzcnt=0;
	double *rhs=new double[maxr];
	char *sense=new char[maxr];
	int *rmatbeg=new int[maxr+1];
	int *rmatind=new int[maxe];
	double *rmatval=new double[maxe];

	rmatbeg[rcnt]=nzcnt;
	sense[rcnt]='E';
	rhs[rcnt]=p0;
	for(int i=0; i<nmodels; i++)
	{
		rmatind[nzcnt]=i;
		rmatval[nzcnt]=1.;
		nzcnt++;
	}
	rcnt++;

	cur=nmodels;
	for(k=0; k<ngraphcomponents; k++)
	{
		for(int i=0; i<subgraphs[k]->nmodels; i++)
		{
			rmatbeg[rcnt]=nzcnt;
			sense[rcnt]='E';
			rhs[rcnt]=1.;
			for(int j=0; j<subgraphs[k]->nearcs[i]; j++)
			{
				rmatind[nzcnt]=cur++;
				rmatval[nzcnt]=1.;
				nzcnt++;
			}
			rmatind[nzcnt]=subgraphs[k]->models[i];
			rmatval[nzcnt]=1.;
			nzcnt++;
			rcnt++;
		}
	}
	cur=nmodels;
	for(k=0; k<ngraphcomponents; k++)
	{
		for(i=0; i<subgraphs[k]->nmodels; i++)
		{
			for(int j=0; j<subgraphs[k]->nearcs[i]; j++)
			{
				rmatbeg[rcnt]=nzcnt;
				sense[rcnt]='L';
				rhs[rcnt]=0.;
				rmatind[nzcnt]=cur++;
				rmatval[nzcnt]=1.;
				nzcnt++;
				rmatind[nzcnt]=subgraphs[k]->models[subgraphs[k]->earcs[i][j].u];
				rmatval[nzcnt]=-1.;
				nzcnt++;
				rcnt++;
			}
		}
	}



	status = CPXaddrows (env, lp, 0, rcnt, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, 0, 0);
	if(status)
		error("status = CPXaddrows (env, lp, 0, rcnt..");
		


	String s;
	sprintf(s, "%s.lp", probname);
	//CPXwriteprob(env, lp, s, 0);


	delete[] rhs;
	delete[] sense;
	delete[] rmatbeg;
	delete[] rmatind;
	delete[] rmatval;
	delete[] obj;
	delete[] lb;
	delete[] ub;
	for( i=0; i<nvar; i++)
	{
		delete[] colname[i];
	}
	delete[] colname;
	delete[] indices;
	delete[] ctype;
*/
	return 0;

}

int ODMProblem::writesol(char* name, int classe)
{
	String s;
	sprintf(s, "%-30s", "Saving solution..");
	print(s);
	if(!solved)
	{
		print("no solution.\n");
		return 0;
	}
	sprintf(s, "%stc.%d", name, classe);
	ofstream to(s);
	int p;
	for(p=p1; p>=p0; p--)
	{
		sprintf(s, "%20f\n", values[p-p0]+constructioncost*p );
		to << s;
	}

	to.close();

	sprintf(s, "%stcFULL.%d", name, classe);
	to.open(s);

	
	for( p=p1; p>=p0; p--)
	{
		double ll = lb[p-p0]+constructioncost*p;
		double bb = values[p-p0]+constructioncost*p;

		sprintf(s, "%4d%20f%20f%%\n", p, bb, 100.*(ll/bb) );
		to << s;
	}

	to.close();

	for(p=p0; p<=p1; p++)
	{
		sprintf(s, "%sout%d.%d", name, classe, p);
		to.open(s);
		for(int i=0; i<nmodels; i++)
		{
			sprintf(s, "%s,%s\n", modelnames[i], modelnames[solutions[p-p0][i]]);
			to << s;
		}
		to.close();
	}

	timer->printTimer();
	print("\n");

	sprintf(s, "%s.sol", name);
	to.open(s);

	
	for( p=p1; p>=p0; p--)
	{
		double ll = lb[p-p0];
		double bb = values[p-p0];

		sprintf(s, "%4d%20f%20f\n", p, bb, ll );
		to << s;
	}

	to.close();
	return 0;
}

int ODMProblem::writesol1(char* name)
{
	String s;
	sprintf(s, "%-30s", "Saving solution..");
	print(s);
	if(!solved)
	{
		print("no solution.\n");
		return 0;
	}
	int p, i, j, k;

	for(p=p0; p<p1; p++)
	{
		double vf = 0.;
		for(j=0; j<nmodels; j++)
		{
			if(j == solutions[p-p0][j])
				continue;
			bool bb = false;
			for(k=0; k<nearcs[j]; k++)
			{
				if(earcs[j][k].u == solutions[p-p0][j])
				{
					vf += earcs[j][k].cost;
					bb = true;
					break;
				}
			}
			//if(!bb)
			//	error("Error in writesol");
		}
		if(fabs(vf - values[p-p0])>1.e-5)
		{
			printf("\nMismatch p=%d vf = %12.2f  values = %12.2f\n", p, vf, values[p-p0]);
		}
	}

	sprintf(s, "%s", name);
	ofstream to(s);
	for(p=p0; p<=p1; p++)
	{
		//sprintf(s, "%20f\n", values[p-p0]+constructioncost*p );
		sprintf(s, "%6d%20f\n", p, values[p-p0]);
		to << s;
	}

	for(p=p0; p<=p1; p++)
	{
		to << "p=" << p << "\n";
		for(int i=0; i<nmodels; i++)
		{
			sprintf(s, "%d,%d\n", i, solutions[p-p0][i]);
			to << s;
		}
	}

	timer->printTimer();
	print("\n");
	to.close();
	return 0;
}


void ODMProblem::reopt()
{
	printf("\nReopt\n");
	int i, j, p, q;
	if(ngraphcomponents > 1)
	{
		double** vv = new double*[ngraphcomponents];
		int** ps = new int*[ngraphcomponents];
		int **fl = new int*[ngraphcomponents];
		int nfl = 0;
		for(i=0; i<ngraphcomponents; i++)
		{
			nfl += subgraphs[i]->p1 - subgraphs[i]->p0 + 1;
			fl[i] = new int[subgraphs[i]->p1 - subgraphs[i]->p0 + 1];
			for(j=subgraphs[i]->p0; j<= subgraphs[i]->p1; j++)
			{
				fl[i][j - subgraphs[i]->p0] = 0;
			}
		}
		for(i=0; i<ngraphcomponents; i++)
		{
			ps[i] = new int[2];
		}
		for(j=0; j<ngraphcomponents; j++)
		{
			//printf("Component %d:\n", j);
			for(i=0; i<ngraphcomponents; i++)
			{
				if(i<j)
				{
					ps[i][0] = subgraphs[i]->p0;
					ps[i][1] = subgraphs[i]->p1;
					if(ps[i][1] > subgraphs[i]->nmodels)
						ps[i][1] = subgraphs[i]->nmodels;
					vv[i] = subgraphs[i]->lb;
				}
				if(i>j)
				{
					ps[i-1][0] = subgraphs[i]->p0;
					ps[i-1][1] = subgraphs[i]->p1;
					if(ps[i-1][1] > subgraphs[i]->nmodels)
						ps[i-1][1] = subgraphs[i]->nmodels;
					vv[i-1] = subgraphs[i]->lb;
				}
			}
			Commom* knapsack = new Commom(ngraphcomponents-1, vv, ps, p0-subgraphs[j]->p0, p1-subgraphs[j]->p0);
			knapsack->DP();

			double *nlb = new double[subgraphs[j]->p1 - subgraphs[j]->p0 + 1];

			for(p = p0; p <= p1; p++)
			{
				//printf(" p=%4d\n",p);
				int q1 = p - p0 + subgraphs[j]->p0;
				if(q1 > subgraphs[j]->nmodels)
					q1 = subgraphs[j]->nmodels;
				for(q=subgraphs[j]->p0; q<=q1; q++)
				{
					double lbv = subgraphs[j]->lb[q-subgraphs[j]->p0] + knapsack->lb[p-q-p0+subgraphs[j]->p0];
					if((values[p-p0] - lbv)/values[p-p0] < OPTGAP)
						continue;
					if(fabs(subgraphs[j]->values[q-subgraphs[j]->p0]-subgraphs[j]->lb[q-subgraphs[j]->p0]) < OPTGAP)
						continue;
					if((subgraphs[j]->values[q-subgraphs[j]->p0]-subgraphs[j]->lb[q-subgraphs[j]->p0])/
						subgraphs[j]->values[q-subgraphs[j]->p0] < OPTGAP)
						continue;
					fl[j][q-subgraphs[j]->p0] = 1;
				}

			}

			delete[] nlb;
			delete knapsack;
		}

		print("Unsolved:\n");
		int unsl = 0;
		char ss[500];
		for(i=0; i<ngraphcomponents; i++)
		{
			for(p=subgraphs[i]->p0; p<=subgraphs[i]->p1; p++)
			{
				if(!fl[i][p-subgraphs[i]->p0])
					continue;
				unsl ++;
			}
		}

		sprintf(ss, "unsl = %d %.2f\n", unsl, (double)unsl/(double) nfl *100.);
		print(ss);
		nprob2 = unsl;
		int k;

#pragma omp parallel for schedule(dynamic,1) private(p)
		for(i=0; i<ngraphcomponents; i++)
		{
			for(p=subgraphs[i]->p0; p<=subgraphs[i]->p1; p++)
			{
				if(!fl[i][p-subgraphs[i]->p0])
					continue;
				sprintf(ss, "C%4.4d p = %d\n", i, p);
				print(ss);
				subgraphs[i]->reopt(p);
			}
		}

		for(i=0; i<ngraphcomponents; i++)
		{
			ps[i][0] = subgraphs[i]->p0;
			ps[i][1] = subgraphs[i]->p1;
			vv[i] = subgraphs[i]->values;
		}
		Commom* knapsack = new Commom(ngraphcomponents, vv, ps, p0, p1);

		knapsack->DP();
		//knapsack->optimize();

		for(p=p0; p<=p1; p++)
		{
			values[p-p0] = knapsack->values[p-p0];

			for( i=0; i<nmodels; i++)
			{
				solutions[p-p0][i]=i;
			}

			int cur=0;
			for(int k=0; k<ngraphcomponents; k++)
			{
				for(int i=0; i<subgraphs[k]->nmodels; i++)
				{
					solutions[p-p0][subgraphs[k]->models[i]] = 
						subgraphs[k]->models[subgraphs[k]->solutions[knapsack->solutions[p-p0][k]
						-subgraphs[k]->p0][i]];
				}
			}
		}

		delete knapsack;

		for( i=0; i<ngraphcomponents; i++)
		{
			ps[i][0] = subgraphs[i]->p0;
			ps[i][1] = subgraphs[i]->p1;
			vv[i] = subgraphs[i]->lb;
		}
	
		knapsack = new Commom(ngraphcomponents, vv, ps, p0, p1);
		knapsack->DP();

		for(p=p0; p<=p1; p++)
		{
			lb[p-p0] = knapsack->lb[p-p0];
		}

		for(i=0; i<ngraphcomponents-1; i++)
		{
			//delete[] vv[i];
			delete[] ps[i];
		}
		delete[] vv;
		delete[] ps;
	}
	else
	{
		int unsl = 0;

		for(p=p0; p<=p1; p++)
		{
			if(p<nmodels)
				unsl += (values[p-p0]-lb[p-p0])/values[p-p0] > OPTGAP;
		}

		char ss[500];
		sprintf(ss, "unsl = %d %.2f\n", unsl, (double)unsl/(double) (p1-p0) *100.);
		print(ss);

		for(p=p0+2; p<=p1; p++)
		{
			if((values[p-p0]-lb[p-p0])/values[p-p0] <= OPTGAP || p >=nmodels)
				continue;
			sprintf(ss, "C%4.4d p = %d\n", 0, p);
			print(ss);
			subgraphs[0]->reopt(p);

			values[p-p0] = subgraphs[0]->values[p-p0];
			lb[p-p0] = subgraphs[0]->lb[p-p0];

			int i;

			for( i=0; i<nmodels; i++)
			{
				solutions[p-p0][i]=i;
			}

			for(int i=0; i<nmodels; i++)
			{
				solutions[p-p0][i] = subgraphs[0]->solutions[p-p0][i];
			}
		}

	}

	print("--------------------------------------------\n");
	print("Results\n");
	char s[500];
	for(p=p0; p<=p1;p++)
	{
		double l1 = lb[p-p0];//+p*constructioncost;
		double l2 = values[p-p0];//+p*constructioncost;
		sprintf(s, "%5d%12.2f%12.2f%12.2f%\%\n", p, l1, l2, 100.*l1/l2);
		print(s);
	}
	print("--------------------------------------------\n");
}

void ODMProblem::optBB()
{
	printf("\nReoptBB\n");
	int i, j, p, q;
	if(ngraphcomponents >1)
	{
		double** vv = new double*[ngraphcomponents];
		int** ps = new int*[ngraphcomponents];
		int **fl = new int*[ngraphcomponents];
		int nfl = 0;
		for(i=0; i<ngraphcomponents; i++)
		{
			nfl += subgraphs[i]->p1 - subgraphs[i]->p0 + 1;
			fl[i] = new int[subgraphs[i]->p1 - subgraphs[i]->p0 + 1];
			for(j=subgraphs[i]->p0; j<= subgraphs[i]->p1; j++)
			{
				fl[i][j - subgraphs[i]->p0] = 0;
			}
		}
		for(i=0; i<ngraphcomponents; i++)
		{
			ps[i] = new int[2];
		}
		for(j=0; j<ngraphcomponents; j++)
		{
			//printf("Component %d:\n", j);
			for(i=0; i<ngraphcomponents; i++)
			{
				if(i<j)
				{
					ps[i][0] = subgraphs[i]->p0;
					ps[i][1] = subgraphs[i]->p1;
					if(ps[i][1] > subgraphs[i]->nmodels)
						ps[i][1] = subgraphs[i]->nmodels;
					vv[i] = subgraphs[i]->lb;
				}
				if(i>j)
				{
					ps[i-1][0] = subgraphs[i]->p0;
					ps[i-1][1] = subgraphs[i]->p1;
					if(ps[i-1][1] > subgraphs[i]->nmodels)
						ps[i-1][1] = subgraphs[i]->nmodels;
					vv[i-1] = subgraphs[i]->lb;
				}
			}
			Commom* knapsack = new Commom(ngraphcomponents-1, vv, ps, p0-subgraphs[j]->p0, p1-subgraphs[j]->p0);
			knapsack->DP();

			double *nlb = new double[subgraphs[j]->p1 - subgraphs[j]->p0 + 1];

			for(p = p0; p <= p1; p++)
			{
				//printf(" p=%4d\n",p);
				int q1 = p - p0 + subgraphs[j]->p0;
				if(q1 > subgraphs[j]->nmodels)
					q1 = subgraphs[j]->nmodels;
				for(q=subgraphs[j]->p0; q<=q1; q++)
				{
					double lbv = subgraphs[j]->lb[q-subgraphs[j]->p0] + knapsack->lb[p-q-p0+subgraphs[j]->p0];
					if((values[p-p0] - lbv)/values[p-p0] < OPTGAP)
						continue;
					if(fabs(subgraphs[j]->values[q-subgraphs[j]->p0]-subgraphs[j]->lb[q-subgraphs[j]->p0]) < OPTGAP)
						continue;
					if((subgraphs[j]->values[q-subgraphs[j]->p0]-subgraphs[j]->lb[q-subgraphs[j]->p0])/
						subgraphs[j]->values[q-subgraphs[j]->p0] < OPTGAP)
						continue;
					fl[j][q-subgraphs[j]->p0] = 1;
				}

			}

			delete[] nlb;
			delete knapsack;
		}

		print("Unsolved:\n");
		int unsl = 0;
		char ss[500];
		for(i=0; i<ngraphcomponents; i++)
		{
			for(p=subgraphs[i]->p0; p<=subgraphs[i]->p1; p++)
			{
				if(!fl[i][p-subgraphs[i]->p0])
					continue;
				unsl ++;
			}
		}
		sprintf(ss, "unsl = %d %.2f\n", unsl, (double)unsl/(double) nfl *100.);
		print(ss);
		nprob3 = unsl;

#pragma omp parallel for schedule(dynamic,1) private(p)
		for(i=0; i<ngraphcomponents; i++)
		{
			for(p=subgraphs[i]->p0; p<=subgraphs[i]->p1; p++)
			{
				if(!fl[i][p-subgraphs[i]->p0])
					continue;
				sprintf(ss, "C%4.4d p = %d\n", i, p);
				print(ss);
				subgraphs[i]->optBB(p);
			}
		}

		for(i=0; i<ngraphcomponents; i++)
		{
			ps[i][0] = subgraphs[i]->p0;
			ps[i][1] = subgraphs[i]->p1;
			vv[i] = subgraphs[i]->values;
		}
		Commom* knapsack = new Commom(ngraphcomponents, vv, ps, p0, p1);

		knapsack->DP();
		//knapsack->optimize();

		for(p=p0; p<=p1; p++)
		{
			values[p-p0] = knapsack->values[p-p0];

			for( i=0; i<nmodels; i++)
			{
				solutions[p-p0][i]=i;
			}

			int cur=0;
			for(int k=0; k<ngraphcomponents; k++)
			{
				for(int i=0; i<subgraphs[k]->nmodels; i++)
				{
					solutions[p-p0][subgraphs[k]->models[i]] = 
						subgraphs[k]->models[subgraphs[k]->solutions[knapsack->solutions[p-p0][k]
						-subgraphs[k]->p0][i]];
				}
			}
		}

		delete knapsack;

		for( i=0; i<ngraphcomponents; i++)
		{
			ps[i][0] = subgraphs[i]->p0;
			ps[i][1] = subgraphs[i]->p1;
			vv[i] = subgraphs[i]->lb;
		}
	
		knapsack = new Commom(ngraphcomponents, vv, ps, p0, p1);
		knapsack->DP();

		for( p=p0; p<=p1; p++)
		{
			lb[p-p0] = knapsack->lb[p-p0];
		}

		for(i=0; i<ngraphcomponents-1; i++)
		{
			//delete[] vv[i];
			delete[] ps[i];
		}
		delete[] vv;
		delete[] ps;
	}
	else
	{
		int unsl = 0;

		for(p=p0; p<=p1; p++)
		{
			if(p<nmodels)
				unsl += (values[p-p0]-lb[p-p0])/values[p-p0] > OPTGAP;
		}

		char ss[500];
		sprintf(ss, "unsl = %d %.2f\n", unsl, (double)unsl/(double) (p1-p0) *100.);
		print(ss);

		for(p=p0+2; p<=p1; p++)
		{
			if((values[p-p0]-lb[p-p0])/values[p-p0] <= OPTGAP || p >=nmodels)
				continue;
			sprintf(ss, "C%4.4d p = %d\n", 0, p);
			print(ss);
			subgraphs[0]->optBB(p);

			values[p-p0] = subgraphs[0]->values[p-p0];
			lb[p-p0] = subgraphs[0]->lb[p-p0];

			int i;

			for( i=0; i<nmodels; i++)
			{
				solutions[p-p0][i]=i;
			}

			for(int i=0; i<nmodels; i++)
			{
				solutions[p-p0][i] = subgraphs[0]->solutions[p-p0][i];
			}
		}
	}

	print("--------------------------------------------\n");
	print("Results\n");
	char s[500];
	for(p=p0; p<=p1;p++)
	{
		double l1 = lb[p-p0];//+p*constructioncost;
		double l2 = values[p-p0];//+p*constructioncost;
		sprintf(s, "%5d%12.2f%12.2f%12.2f%\%\n", p, l1, l2, 100.*l1/l2);
		print(s);
	}
	print("--------------------------------------------\n");
}

int ODMProblem::writegraph1()
{
	String s;
	sprintf(s, "%-30s", "Saving graph..");
	print(s);
	sprintf(s, "%s.graph", probname);
	ofstream to(s);

	to << "digraph G{\noverlap=false\n";

	for(int k=0; k<ngraphcomponents; k++)
	{
		for(int i=0; i<subgraphs[k]->nmodels; i++)
		{
			//sprintf(s, "%d [label=\"%d:%s\"]\n", 
			//	subgraphs[k]->models[i], subgraphs[k]->models[i], 
			//	modelnames[subgraphs[k]->models[i]]);
			//to << s;
			for(int j=0; j<subgraphs[k]->nearcs[i]; j++)
			{
				sprintf(s, "%d -> %d [label=\"%.2f\"]\n", 
					subgraphs[k]->models[subgraphs[k]->earcs[i][j].u], 
					subgraphs[k]->models[i], subgraphs[k]->earcs[i][j].cost);
				to << s;
			}
		}
	}

	to << "}\n";
	to.close();
	timer->printTimer();
	print("\n");
	return 0;
}

int ODMProblem::writegraph()
{
	String s;
	sprintf(s, "%-30s", "Saving graph..");
	print(s);
	sprintf(s, "%s.txt", probname);
	for(int i=0; i<strlen(s); i++)
		if(s[i] == '/')
			s[i] = '-';

	ofstream to(s);

	to << nmodels << " " << narcs << "\n";

	for(int k=0; k<ngraphcomponents; k++)
	{
		for(int i=0; i<subgraphs[k]->nmodels; i++)
		{
			//sprintf(s, "%d [label=\"%d:%s\"]\n", 
			//	subgraphs[k]->models[i], subgraphs[k]->models[i], 
			//	modelnames[subgraphs[k]->models[i]]);
			//to << s;
			for(int j=0; j<subgraphs[k]->nearcs[i]; j++)
			{
				sprintf(s, "%d %d %.4f\n", 
					subgraphs[k]->models[subgraphs[k]->earcs[i][j].u] + 1, 
					subgraphs[k]->models[i] + 1, subgraphs[k]->earcs[i][j].cost);
				to << s;
			}
		}
	}

	to.close();
	timer->printTimer();
	print("\n");
	return 0;
}

int ODMProblem::createKnapsackLP(bool llb)
{
/*	nknapvars=0;
	for(int k=0; k<ngraphcomponents; k++)
	{
		nknapvars+= subgraphs[k]->p1 - subgraphs[k]->p0 +1;
	}

	int status=0;
	double* obj=new double[nknapvars];;
	double* lb=new double[nknapvars];;
	double* ub=new double[nknapvars];
	char** colname=new char*[nknapvars];
	int* indices=new int[nknapvars];
	char* ctype=new char[nknapvars];

	int cur=0;
	for( k=0; k<ngraphcomponents; k++)
	{
		for(int p=subgraphs[k]->p0;p<=subgraphs[k]->p1;p++)
		{
			obj[cur]=subgraphs[k]->values[p-subgraphs[k]->p0];
			if(llb)
			{
				obj[cur]=subgraphs[k]->lb[p-subgraphs[k]->p0];
			}
			lb[cur]=0.;
			ub[cur]=1.;
			colname[cur]=new char[50];
			sprintf(colname[cur],"z%d,%d", k, p);
			//cout << subgraphs[k]->models[i] << endl;
			indices[cur]=cur;
			ctype[cur]='B';
			cur++;
		}
	}
	
	status=CPXaddcols (env, lp, cur, 0, obj, 0, 0, 0, lb, ub,  colname);
	if(status)
		error("status=CPXaddcols (env, lp, n, 0, obj, 0, 0, 0, lb, ub, colname);");

	status=CPXchgctype (env, lp, cur, indices, ctype);
	if(status)
		error("status=CPXchgctype (env, lp, nvar, indices, ctype);");


	int maxr=ngraphcomponents+1;
	int maxe=2*nknapvars;
	int rcnt=0;
	int nzcnt=0;
	double *rhs=new double[maxr];
	char *sense=new char[maxr];
	int *rmatbeg=new int[maxr+1];
	int *rmatind=new int[maxe];
	double *rmatval=new double[maxe];

	

	cur=0;
	rmatbeg[rcnt]=nzcnt;
	sense[rcnt]='E';
	rhs[rcnt]=p0;
	for( k=0; k<ngraphcomponents; k++)
	{
		for(int p=subgraphs[k]->p0;p<=subgraphs[k]->p1;p++)
		{
			rmatind[nzcnt]=cur++;
			rmatval[nzcnt]=p;
			nzcnt++;
		}
	}
	rcnt++;

	cur=0;
	rmatbeg[rcnt]=nzcnt;
	sense[rcnt]='E';
	rhs[rcnt]=p0;
	for( k=0; k<ngraphcomponents; k++)
	{
		rmatbeg[rcnt]=nzcnt;
		sense[rcnt]='E';
		rhs[rcnt]=1;
		for(int p=subgraphs[k]->p0;p<=subgraphs[k]->p1;p++)
		{
			rmatind[nzcnt]=cur++;
			rmatval[nzcnt]=1.;
			nzcnt++;
		}
		rcnt++;
	}
	

	status = CPXaddrows (env, lp, 0, rcnt, nzcnt, rhs, sense, rmatbeg, rmatind, rmatval, 0, 0);
	if(status)
		error("status = CPXaddrows (env, lp, 0, rcnt..");
		


	String s;
	sprintf(s, "%skn.lp", probname);
	//CPXwriteprob(env, lp, s, 0);


	delete[] rhs;
	delete[] sense;
	delete[] rmatbeg;
	delete[] rmatind;
	delete[] rmatval;
	delete[] obj;
	delete[] lb;
	delete[] ub;
	for(int i=0; i<nknapvars; i++)
	{
		delete[] colname[i];
	}
	delete[] colname;
	delete[] indices;
	delete[] ctype;
*/
	return 0;

}

int ODMProblem::greedyFull()
{

	soltimer->start(true);
	String s;
	sprintf(s, "%-30s\n", "Looking for solution by greedy..");
	print(s);
	sprintf(s, "%-30s", "Arcs sorting..");
	print(s);
	arcsorting();
	timer->printTimer();
	print("\n");
	
	sprintf(s, "%4d-%4d: init", p0, p1);
	print(s);
	if(nmodels>2000)
	{
		print(" it is big. be patient.");
	}

	int size=sizeof(int)*8;
	int ncols=nmodels/size+(nmodels%size>0);
	unsigned int** adj=new unsigned int*[nmodels];
	int i;
	for( i=0; i<nmodels; i++)
	{
		adj[i]=new unsigned int[ncols];
		for(int j=0; j<ncols; j++)
		{
			adj[i][j]=0;
		}
	}


	int* nearcs=new int[nmodels];
	Arc1** earcs=new Arc1*[nmodels];
	int k;
	for( k=0; k<ngraphcomponents; k++)
	{
		for(int i=0; i<subgraphs[k]->nmodels; i++)
		{
			nearcs[subgraphs[k]->models[i]]=subgraphs[k]->nearcs[i];
			earcs[subgraphs[k]->models[i]]=new Arc1[subgraphs[k]->nearcs[i]];
		}
	}

	for( k=0; k<ngraphcomponents; k++)
	{
		for(int i=0; i<subgraphs[k]->nmodels; i++)
		{
			for(int j=0; j<subgraphs[k]->nearcs[i]; j++)
			{
				earcs[subgraphs[k]->models[i]][j].u=
					subgraphs[k]->models[subgraphs[k]->earcs[i][j].u];
				//earcs[subgraphs[k]->models[i]][j].v=
				//	subgraphs[k]->models[subgraphs[k]->earcs[i][j].v];
				earcs[subgraphs[k]->models[i]][j].cost=subgraphs[k]->earcs[i][j].cost;
			}
		}
	}
	for( i=0; i<nmodels; i++)
	{
		for(int j=0; j<nearcs[i]; j++)
		{
			adj[earcs[i][j].u][i/size]+= (1<<(i%size));
		}
	}
	
	/*{
		String s;
		print("\nadj\n");
		for(int i=0; i<nmodels; i++)
		{
			for(int j=0; j<nmodels; j++)
			{
				sprintf(s, "%2d", (adj[i][j/size]&(1<<j%size))>0);
				print(s);
			}
			print("\n");
		}
	}*/


	sprintf(s, " f(%d)", p0);
	print(s);

	int* cursol=new int[nmodels];
	
	for( i=0; i<nmodels; i++)
	{
		cursol[i]=-1;
	}

	double objval=0.;
	for(i=0; i<nmodels; i++)
	{
		if(nearcs[i]==0)
			continue;
		for(int j=0; j<nearcs[i]; j++)
		{
			if(nearcs[earcs[i][j].u]==0)
			{
				cursol[i]=j;
				objval+= earcs[i][j].cost;
				break;
			}
		}
	}

	values[0]=objval;
	sprintf(s, "=%.2f", values[0]+p0*constructioncost);
	print(s);
	for(i=0; i<nmodels; i++)
	{
		if(cursol[i]==-1)
		{
			solutions[0][i]=i;
		}
		else
		{
			solutions[0][i]=earcs[i][cursol[i]].u;
		}
	}

	for(int p=p0+1; p<=p1; p++)
	{
		double min=1;
		double curval=0;
		int best=-1;
		sprintf(s, " f(%d)", p);
		print(s);

		int i;
		for( i=0; i<nmodels; i++)
		{
			if(cursol[i]==-1)
				continue;
			curval = -earcs[i][cursol[i]].cost;
			for(int j=0; j<nmodels; j++)
			{
				if(i==j || cursol[j]==-1)
					continue;
				if((adj[i][j/size]&(1<<j%size))==0)
					continue;
				int kl=-1;
				//for(int k=0; k<cursol[j]; k++)
				for(int k=cursol[j]-1;k>=0; k--)
				{
					if(earcs[j][k].u==i)
					{
						kl=k;
						break;
					}
				}
				if(kl!=-1)
				{
					curval+= earcs[j][kl].cost - earcs[j][cursol[j]].cost;
				}
			}
			if(min>curval)
			{
				min=curval;
				best=i;
			}
		}

		values[p-p0]=values[p-p0-1]+min;
		cursol[best]=-1;
		
		for(int j=0; j<nmodels; j++)
		{
			if(i==j || cursol[j]==-1)
				continue;
			for(int k=0; k<cursol[j]; k++)
			{
				if(earcs[j][k].u==best)
				{
					cursol[j]=k;
					break;
				}
			}
		}

		sprintf(s, "=%.2f", values[p-p0]+p*constructioncost);
		print(s);
		for(i=0; i<nmodels; i++)
		{
			if(cursol[i]==-1)
			{
				solutions[p-p0][i]=i;
			}
			else
			{
				solutions[p-p0][i]=earcs[i][cursol[i]].u;
			}
		}

		
	}

	for( i=0; i<nmodels; i++)
	{
		delete[] earcs[i];
		delete[] adj[i];
	}
	delete[] earcs;
	delete[] adj;
	delete[] nearcs;
	delete[] cursol;
	print("\n");
	soltimer->stop();
	sprintf(s, "%-30s%15d\n", "Solution time (secs)", soltimer->seconds());
	print(s);
	sprintf(s, "%-30s", "Solution time ");
	print(s);
	soltimer->printTimer();
	print("\n");
	solved=true;
	return 0;
	
}

int ODMProblem::arcsorting()
{
	for(int k=0; k<ngraphcomponents; k++)
	{
		subgraphs[k]->arcsorting();
	}

	return 0;
}

int ODMProblem::greedy()
{
	soltimer->start(true);
	String s;
	sprintf(s, "%-30s\n", "Looking for solution by greedy with decomposition..");
	print(s);

	int k;
	for( k=0; k<ngraphcomponents; k++)
	{
		subgraphs[k]->greedy();
	}

	sprintf(s, "Knapsack %4d-%4d: init", p0, p1);
	print(s);
	int* comp=new int[ngraphcomponents];
	for( k=0; k<ngraphcomponents; k++)
	{
		comp[k]=subgraphs[k]->p0;
	}
	sprintf(s, " f(%d)", p0);
	print(s);
	int p;
	for( p=p0; p<=p1; p++)
	{
		values[p-p0]=0;
		int k;
		for(k=0; k<ngraphcomponents; k++)
		{
			values[p-p0]+= subgraphs[k]->values[comp[k]-subgraphs[k]->p0];
		}
		sprintf(s, "=%.2f", p*constructioncost+values[p-p0]);
		print(s);

		for( k=0; k<ngraphcomponents; k++)
		{
			for(int i=0; i<subgraphs[k]->nmodels; i++)
			{
				solutions[p-p0][subgraphs[k]->models[i]] = 
					subgraphs[k]->models[subgraphs[k]->solutions[comp[k]-subgraphs[k]->p0][i]];
			}
		}

		if(p==p1)
			break;
		sprintf(s, " f(%d)", p);
		print(s);

		double min=1.;
		int best=-1;
		for(k=0; k<ngraphcomponents; k++)
		{
			if(comp[k]!=subgraphs[k]->p1 && 
				min> subgraphs[k]->values[comp[k]+1-subgraphs[k]->p0] - subgraphs[k]->values[comp[k]-subgraphs[k]->p0])
			{
				min=subgraphs[k]->values[comp[k]+1-subgraphs[k]->p0] - subgraphs[k]->values[comp[k]-subgraphs[k]->p0];
				best=k;
			}
		}
		comp[best]++;
	}

	delete[] comp;
	print("\n");
	soltimer->stop();
	sprintf(s, "%-30s%15d\n", "Solution time (secs)", soltimer->seconds());
	print(s);
	sprintf(s, "%-30s", "Solution time ");
	print(s);
	soltimer->printTimer();
	print("\n");
	solved=true;
	return 0;
}



int ODMProblem::LH(bool high)
{
/*	soltimer->start(true);
	String s;
	sprintf(s, "%-30s\n", "Looking for solution by Lagrangean heuristic..");
	print(s);

	for(int k=0; k<ngraphcomponents; k++)
	{
		subgraphs[k]->greedy(true, high);
	}

	sprintf(s, "Knapsack %4d-%4d: init", p0, p1);
	print(s);
	int status=0;
	env = CPXopenCPLEX (&status);
	if(!env)
		error("env = CPXopenCPLEX (&status);");
	// create problem
	lp=CPXcreateprob(env,&status,"ll");
	if(!lp)
		error("lp==CPXcreateprob(env,&status,k);");

	status = CPXchgprobtype (env, lp, 1);
	if(status)
		error("status = CPXchgprobtype (env, lp, CPXPROB_MIP);");

	//CPXsetintparam (env, CPX_PARAM_SCRIND, 1);


	
	
	createKnapsackLP();
	double* x=new double[nknapvars];

	for(int p=p0; p<=p1; p++)
	{
		sprintf(s, " f(%d)", p);
		print(s);
		int cnt=1;
		int indices[1]={0};
		double rhs[1]={p};


		status = CPXchgrhs (env, lp, cnt, indices, rhs);
		if(status)
		{
			error("status = CPXchgrhs (env, lp, cnt, indices, rhs);");
		}

		status = CPXmipopt(env, lp);
		if(status)
			error("status = CPXmipopt(env, lp);");

		status = CPXgetmipobjval (env, lp, values+p-p0);
		if(status)
			error("status = CPXgetmipobjval (env, lp, values+p-p0);");
		sprintf(s, "=%.2f", p*constructioncost+values[p-p0]);
/*		print(s);

		status = CPXgetmipx(env, lp, x, 0, nknapvars-1);
		if(status)
			error("status = CPXgetmipx(env, lp, x, 0, nknapvars-1);");


		for(int i=0; i<nmodels; i++)
		{
			solutions[p-p0][i]=i;
		}

		int cur=0;
		for(int k=0; k<ngraphcomponents; k++)
		{
			for(int pp=subgraphs[k]->p0; pp<=subgraphs[k]->p1; pp++)
			{
				if(x[cur++]<0.5)
					continue;
				for(int i=0; i<subgraphs[k]->nmodels; i++)
				{
					solutions[p-p0][subgraphs[k]->models[i]] = 
						subgraphs[k]->models[subgraphs[k]->solutions[pp-subgraphs[k]->p0][i]];
				}
			}
		}
	}


	status=0;
	status=CPXfreeprob(env,&lp);
	if(status)
		error("status=CPXfreeprob(env,&lp);");
	status = CPXcloseCPLEX (&env);
	if(status)
		error("status = CPXcloseCPLEX (&env);");


	sprintf(s, "\nLB Knapsack %4d-%4d: init", p0, p1);
	print(s);
	status=0;
	env = CPXopenCPLEX (&status);
	if(!env)
		error("env = CPXopenCPLEX (&status);");
	// create problem
	lp=CPXcreateprob(env,&status,"ll");
	if(!lp)
		error("lp==CPXcreateprob(env,&status,k);");

	status = CPXchgprobtype (env, lp, 1);
	if(status)
		error("status = CPXchgprobtype (env, lp, CPXPROB_MIP);");

	//CPXsetintparam (env, CPX_PARAM_SCRIND, 1);


	
	
	createKnapsackLP(true);
	//double* x=new double[nknapvars];

	for( p=p0; p<=p1; p++)
	{
		sprintf(s, " LB(%d)", p);
		print(s);
		int cnt=1;
		int indices[1]={0};
		double rhs[1]={p};


		status = CPXchgrhs (env, lp, cnt, indices, rhs);
		if(status)
		{
			error("status = CPXchgrhs (env, lp, cnt, indices, rhs);");
		}

		status = CPXmipopt(env, lp);
		if(status)
			error("status = CPXmipopt(env, lp);");

		status = CPXgetmipobjval (env, lp, lb+p-p0);
		if(status)
			error("status = CPXgetmipobjval (env, lp, lb+p-p0);");
		sprintf(s, "=%.2f", p*constructioncost+values[p-p0]);
/*		print(s);

	}


	status=0;
	status=CPXfreeprob(env,&lp);
	if(status)
		error("status=CPXfreeprob(env,&lp);");
	status = CPXcloseCPLEX (&env);
	if(status)
		error("status = CPXcloseCPLEX (&env);");
	
	

	print("\n");
	soltimer->stop();
	sprintf(s, "%-30s%15d\n", "Solution time (secs)", soltimer->seconds());
	print(s);
	sprintf(s, "%-30s", "Solution time ");
	print(s);
	soltimer->printTimer();
	print("\n");
	solved=true;
	delete[] x;


	// testing lagrangean for knapsack

	double** vv = new double*[ngraphcomponents];
	int** ps = new int*[ngraphcomponents];
	for(int i=0; i<ngraphcomponents; i++)
	{
		ps[i] = new int[2];
		ps[i][0] = subgraphs[i]->p0;
		ps[i][1] = subgraphs[i]->p1;
		vv[i] = subgraphs[i]->values;
	}
	Commom* knapsack = new Commom(ngraphcomponents, vv, ps, p0, p1);
	knapsack->printdata();
	knapsack->optimize();

	for(i=0; i<ngraphcomponents; i++)
	{
		delete[] ps[i];
	}
	delete[] ps;
	delete[] vv;

	delete knapsack;
*/
	return 0;
}

int ODMProblem::optimize(int lh, int level)
{

	soltimer->start(true);
	String s;
	sprintf(s, "%-30s\n", "Optimizing..");
	print(s);
	sprintf(s, "Method: gready");
	print(s);
	if(lh)
	{
		print("+lagrangean heuristic ");
		switch(level)
		{
		case 0:
			print("fast and low quality\n");
		case 1:
			print("medium velocity and quality\n");
		case 2:
			print("slow and high quality\n");
		}
	}
	else
	{
		print("\n");
	}
//#pragma omp parallel for
	for(int k=0; k<ngraphcomponents; k++)
	{
		//subgraphs[k]->printg();
		subgraphs[k]->greedy(lh, level);
	}
	//error("Stop");
	// testing lagrangean for knapsack

	double** vv = new double*[ngraphcomponents];
	int** ps = new int*[ngraphcomponents];
	int i;
	for(i=0; i<ngraphcomponents; i++)
	{
		ps[i] = new int[2];
		ps[i][0] = subgraphs[i]->p0;
		ps[i][1] = subgraphs[i]->p1;
		vv[i] = subgraphs[i]->values;
	}
	Commom* knapsack = new Commom(ngraphcomponents, vv, ps, p0, p1);

	//knapsack->printdata();
	//knapsack->optimize();
	knapsack->DP();


	for(int p=p0; p<=p1; p++)
	{
		values[p-p0] = knapsack->values[p-p0];

		for( i=0; i<nmodels; i++)
		{
			solutions[p-p0][i]=i;
		}

		int cur=0;
		for(int k=0; k<ngraphcomponents; k++)
		{
			for(int i=0; i<subgraphs[k]->nmodels; i++)
			{
				solutions[p-p0][subgraphs[k]->models[i]] = 
					subgraphs[k]->models[subgraphs[k]->solutions[knapsack->solutions[p-p0][k]
					-subgraphs[k]->p0][i]];
			}
		}
	}

	for(i=0; i<ngraphcomponents; i++)
	{
		delete[] ps[i];
	}
	delete[] ps;
	delete[] vv;

	delete knapsack;

	
	if(lh)
	{
		double** vv = new double*[ngraphcomponents];
		int** ps = new int*[ngraphcomponents];
		int i;
		for( i=0; i<ngraphcomponents; i++)
		{
			ps[i] = new int[2];
			ps[i][0] = subgraphs[i]->p0;
			ps[i][1] = subgraphs[i]->p1;
			vv[i] = subgraphs[i]->lb;
		}
		Commom* knapsack = new Commom(ngraphcomponents, vv, ps, p0, p1);
		//knapsack->printdata();
		//knapsack->optimize();
		knapsack->DP();

		int p;
		for( p=p0; p<=p1; p++)
		{
			lb[p-p0] = knapsack->lb[p-p0];
		}

		for(i=0; i<ngraphcomponents; i++)
		{
			delete[] ps[i];
		}
		delete[] ps;
		delete[] vv;
		delete knapsack;

		print("--------------------------------------------\n");
		print("Guarantee\n");
		for(p=p0; p<=p1;p++)
		{
			double l1 = lb[p-p0];//+p*constructioncost;
			double l2 = values[p-p0];//+p*constructioncost;
			sprintf(s, "%4d  %12.2f%%\n", p, 100.*l1/l2);
			print(s);
		}
		print("--------------------------------------------\n");
	}
	sprintf(s, "%-30s", "Solution time");
	print(s);
	soltimer->printTimer();
	print("\n");
	solved = true;
	return 0;
}

int ODMProblem::optimize()
{

	soltimer->start(true);
	String s;
	sprintf(s, "%-30s\n", "Optimizing by new version..");
	print(s);
	
	int k;
	IT * it = new IT[ngraphcomponents];
	for(k=0; k<ngraphcomponents; k++)
	{
		it[k].i = k;
		it[k].z = subgraphs[k]->narcs;
	}
 
	qsort(it, ngraphcomponents, sizeof(IT), compareITZA);
	int h;
#pragma omp parallel for reduction(+:nprob1) schedule(dynamic,1) private(h)
	for(k=0; k<ngraphcomponents; k++)
	{
		h = it[k].i;
		nprob1 += subgraphs[h]->p1 - subgraphs[h]->p0 + 1;
		subgraphs[h]->optimize();
		//subgraphs[k]->greedy(1,4);
	}


	delete[] it;

	if(ngraphcomponents > 1)
	{
	
		double** vv = new double*[ngraphcomponents];
		int** ps = new int*[ngraphcomponents];
		int i;
		for(i=0; i<ngraphcomponents; i++)
		{
			ps[i] = new int[2];
			ps[i][0] = subgraphs[i]->p0;
			ps[i][1] = subgraphs[i]->p1;
			vv[i] = subgraphs[i]->values;
		}
		Commom* knapsack = new Commom(ngraphcomponents, vv, ps, p0, p1);

		knapsack->DP();
		//knapsack->optimize();

		for(int p=p0; p<=p1; p++)
		{
			values[p-p0] = knapsack->values[p-p0];

			for( i=0; i<nmodels; i++)
			{
				solutions[p-p0][i]=i;
			}

			int cur=0;
			for(int k=0; k<ngraphcomponents; k++)
			{
				for(int i=0; i<subgraphs[k]->nmodels; i++)
				{
					solutions[p-p0][subgraphs[k]->models[i]] = 
						subgraphs[k]->models[subgraphs[k]->solutions[knapsack->solutions[p-p0][k]
						-subgraphs[k]->p0][i]];
				}
			}
		}

		delete knapsack;
	
		for( i=0; i<ngraphcomponents; i++)
		{
			ps[i] = new int[2];
			ps[i][0] = subgraphs[i]->p0;
			ps[i][1] = subgraphs[i]->p1;
			vv[i] = subgraphs[i]->lb;
		}
	
		knapsack = new Commom(ngraphcomponents, vv, ps, p0, p1);
		knapsack->DP();
		//knapsack->optimize();

		int p;
		for( p=p0; p<=p1; p++)
		{
			lb[p-p0] = knapsack->lb[p-p0];
		}

		for(i=0; i<ngraphcomponents; i++)
		{
			delete[] ps[i];
		}
		delete[] ps;
		delete[] vv;
		delete knapsack;
	}
	else
	{
		for(int p=p0; p<=p1; p++)
		{
			values[p-p0] = subgraphs[0]->values[p-p0];
			lb[p-p0] = subgraphs[0]->lb[p-p0];

			int i;

			for( i=0; i<nmodels; i++)
			{
				solutions[p-p0][i]=i;
			}

			for(int i=0; i<nmodels; i++)
			{
				solutions[p-p0][i] = subgraphs[0]->solutions[p-p0][i];
			}
		}
	}

	print("--------------------------------------------\n");
	print("Results\n");
	for(int p=p0; p<=p1;p++)
	{
		double l1 = lb[p-p0];//+p*constructioncost;
		double l2 = values[p-p0];//+p*constructioncost;
		sprintf(s, "%5d%12.2f%12.2f%12.2f%\%\n", p, l1, l2, 100.*l1/l2);
		print(s);
	}
	print("--------------------------------------------\n");
	sprintf(s, "%-30s", "Solution time");
	print(s);
	soltimer->printTimer();
	print("\n");
	solved = true;
	return 0;
}

int ODMProblem::optimizeParInit(int lh, int level)
{

	String s;
	sprintf(s, "%-30s\n", "Optimizing..");
	print(s);
	sprintf(s, "Method: gready");
	print(s);
	if(lh)
	{
		print("+lagrangean heuristic ");
		switch(level)
		{
		case 0:
			print("fast and low quality\n");
		case 1:
			print("medium velocity and quality\n");
		case 2:
			print("slowt and high quality\n");
		}
	}
	else
	{
		print("\n");
	}

	return 0;
}

int ODMProblem::optimizeParFinal(int lh, int level)
{

	String s;
	// testing lagrangean for knapsack

	double** vv = new double*[ngraphcomponents];
	int** ps = new int*[ngraphcomponents];
	int i;
	for(i=0; i<ngraphcomponents; i++)
	{
		ps[i] = new int[2];
		ps[i][0] = subgraphs[i]->p0;
		ps[i][1] = subgraphs[i]->p1;
		vv[i] = subgraphs[i]->values;
	}
	Commom* knapsack = new Commom(ngraphcomponents, vv, ps, p0, p1);
	//knapsack->printdata();
	knapsack->optimize();

	print("optimized\n");

	for(int p=p0; p<=p1; p++)
	{

		values[p-p0] = knapsack->values[p-p0];
		sprintf(s, "values[%d-%d]=%.2f\n", p, p0, values[p-p0]);
		print(s);
		for( i=0; i<nmodels; i++)
		{
			solutions[p-p0][i]=i;
		}

		int cur=0;
		for(int k=0; k<ngraphcomponents; k++)
		{
			for(int i=0; i<subgraphs[k]->nmodels; i++)
			{
				solutions[p-p0][subgraphs[k]->models[i]] = 
					subgraphs[k]->models[subgraphs[k]->solutions[knapsack->solutions[p-p0][k]
					-subgraphs[k]->p0][i]];
			}
		}
		print("end\n");
	}


	for(i=0; i<ngraphcomponents; i++)
	{
		delete[] ps[i];
	}
	delete[] ps;
	delete[] vv;

	delete knapsack;

	
	print("solution got\n");

	if(lh)
	{
		double** vv = new double*[ngraphcomponents];
		int** ps = new int*[ngraphcomponents];
		int i;
		for( i=0; i<ngraphcomponents; i++)
		{
			ps[i] = new int[2];
			ps[i][0] = subgraphs[i]->p0;
			ps[i][1] = subgraphs[i]->p1;
			vv[i] = subgraphs[i]->lb;
		}
		Commom* knapsack = new Commom(ngraphcomponents, vv, ps, p0, p1);
		//knapsack->printdata();
		knapsack->optimize();

		int p;
		for( p=p0; p<=p1; p++)
		{
			lb[p-p0] = knapsack->lb[p-p0];
		}

		for(i=0; i<ngraphcomponents; i++)
		{
			delete[] ps[i];
		}
		delete[] ps;
		delete[] vv;
		delete knapsack;

		print("--------------------------------------------\n");
		print("Guarantee\n");
		for(p=p0; p<=p1;p++)
		{
			double l1 = lb[p-p0]+p*constructioncost;
			double l2 = values[p-p0]+p*constructioncost;
			sprintf(s, "%4d  %12.2f%%\n", p, 100.*l1/l2);
			print(s);
		}
		print("--------------------------------------------\n");
	}
	sprintf(s, "%-30s", "Solution time");
	print(s);
	soltimer->printTimer();
	print("\n");
	solved = true;

	return 0;
}
