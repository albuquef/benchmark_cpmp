// GraphComponent.cpp: implementation of the GraphComponent class.
//
//////////////////////////////////////////////////////////////////////

#include "GraphComponent.h"
#include "ODMProblem.h"
#include "Convert.h"
#include "Pmed.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

int compareRoAZ( const void *arg1, const void *arg2 )
{

	if(((Ro*) arg1)->fix==1)
		return -1;
	if(((Ro*) arg2)->fix==1)
		return 1;
	if(((Ro*) arg1)->value>((Ro*) arg2)->value) 
		return 1;
	else
		return -1;
	
}

GraphComponent::GraphComponent(ODMProblem* prob_, int cind_, Arc1 **earcs_, int *nearcs_): prob(prob_), cind(cind_)
{
	nmodels=0;
	int i, j;
	for( i=0; i<prob->nmodels; i++)
	{
		if(prob->modelcomponents[i]==cind)
			nmodels++;
	}
	models=new int[nmodels];
	nmodels=0;
	for( i=0; i<prob->nmodels; i++)
	{
		if(prob->modelcomponents[i]==cind)
			models[nmodels++]=i;
	}

	int *invmap = new int[prob->nmodels];
	for(i=0; i<nmodels; i++)
		invmap[models[i]] = i;


	narcs=0;
	nearcs=new int[nmodels];
	earcs=new Arc1*[nmodels];

	for(i=0; i<nmodels; i++)
	{
		nearcs[i]=nearcs_[models[i]];
		earcs[i] = new Arc1[nearcs[i]];
		memcpy(earcs[i], earcs_[models[i]], nearcs[i] * sizeof(Arc1));
		for(j=0; j<nearcs[i]; j++)
		{
			earcs[i][j].u = invmap[earcs[i][j].u];
		}
		narcs += nearcs[i];
	}

	delete[] invmap;

	p0=0;
	for(i=0; i<nmodels; i++)
	{
		p0+= nearcs[i]==0;
	}
	arcsorting();

	solutions = NULL;
	values = NULL;
	lb = NULL;
	
	pmed = new Pmed();
	pmed->p0 = p0;
	getp1();
	pmed->readprob(nmodels, nearcs, earcs);
	pmed->lang_init();
}

GraphComponent::GraphComponent(ODMProblem* prob_, int cind_): prob(prob_), cind(cind_)
{
	nmodels=0;
	int i;
	for( i=0; i<prob->nmodels; i++)
	{
		if(prob->modelcomponents[i]==cind)
			nmodels++;
	}
	models=new int[nmodels];
	nmodels=0;
	for( i=0; i<prob->nmodels; i++)
	{
		if(prob->modelcomponents[i]==cind)
			models[nmodels++]=i;
	}

	narcs=0;
	nearcs=new int[nmodels];
	for(i=0; i<nmodels; i++)
	{
		nearcs[i]=0;
	}

	for(i=0; i<nmodels; i++)
	{
		int j;
		for( j=0; j<nmodels; j++)
		{
			if(i==j)
				continue;
			if( prob->modelhardoptions[models[i]]==prob->modelhardoptions[models[j]] &&
				(prob->modeloptions[models[i]] & prob->modeloptions[models[j]]) 
				== prob->modeloptions[models[j]])
			{
				narcs++;
				nearcs[j]++;
			}
		}
	}
	earcs=new Arc1*[nmodels];
	for( i=0; i<nmodels; i++)
	{
		earcs[i]=new Arc1[nearcs[i]];
		nearcs[i]=0;
	}
	for(i=0; i<nmodels; i++)
	{
		int j;
		for( j=0; j<nmodels; j++)
		{
			if(i==j)
				continue;
			if( prob->modelhardoptions[models[i]]==prob->modelhardoptions[models[j]] &&
				(prob->modeloptions[models[i]] & prob->modeloptions[models[j]]) 
				== prob->modeloptions[models[j]])
			{
				earcs[j][nearcs[j]].u=i;
				//earcs[j][nearcs[j]].v=j;
				earcs[j][nearcs[j]].cost=(prob->modelcosts[models[i]]-prob->modelcosts[models[j]])*
					prob->modeldemands[models[j]];
				nearcs[j]++;
			}
		}
	}
	p0=0;
	for(i=0; i<nmodels; i++)
	{
		p0+= nearcs[i]==0;
	}

	solutions = NULL;
	values = NULL;
	lb = NULL;

	
}

GraphComponent::~GraphComponent()
{
	delete[] models;
	delete[] nearcs;
	delete[] values;
	delete[] lb;
	if(solutions)
	{
		for(int i=0; i<=p1-p0; i++)
		{
			delete[] solutions[i];
		}
		delete[] solutions;
	}
	for(int i=0; i<nmodels; i++)
		delete[] earcs[i];
	delete[] earcs;

	pmed->lang_final();
	delete pmed;
}

int GraphComponent::getp1()
{

	p1=prob->p1-prob->p0+p0;
	p1=MIN(nmodels, p1);
	pmed->p1 = p1;

	return 0;
}

int GraphComponent::optCplex()
{


	return 0;
	
}

int GraphComponent::createLP()
{
	
	return 0;

}

int GraphComponent::arcsorting()
{

	for(int i=0; i<nmodels; i++)
	{
		qsort(earcs[i], nearcs[i], sizeof(Arc1), compareArcCostArrayAZ);
		//std::sort(earcs[i], earcs[i] + nearcs[i], &compareArcCostArrayAZ00);
	}

	return 0;
}


int GraphComponent::greedy(int lh, int high)
{
	String s;
	sprintf(s, "%-30s\n", "Looking for solution by greedy..");
	sprintf(s, "C%4.4d %3d-%3d: init", cind, p0, p1);
	print(s);
	arcsorting();
	
	
	if(nmodels>2000)
	{
		print(" it is big. be patient.");
	}
	solutions=new int*[p1-p0+1];
	int p;
	for( p=p0; p<=p1; p++)
	{
		solutions[p-p0]=new int[nmodels];
	}
	values=new double[p1-p0+1];
	lb=new double[p1-p0+1];

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


	//sprintf(s, "%d=%d\n", models[561], models[2127]);
	//print(s);

	
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


	//sprintf(s, " f(%d)", p0);
	//print(s);

	int* cursol=new int[nmodels];
	//int* seconds=new int[nmodels];

	for( i=0; i<nmodels; i++)
	{
		cursol[i]=-1;
	}

	double objval=0.;
	int cnt =0;
	for(i=0; i<nmodels; i++)
	{
		//seconds[i]=-1;
		if(nearcs[i]==0)
		{
			cnt++;
			continue;
		}
		for(int j=0; j<nearcs[i]; j++)
		{
			if(nearcs[earcs[i][j].u]==0)
			{
				cursol[i]=j;
				objval+= earcs[i][j].cost;
				break;
			}
		}
		//sprintf(s, "i=%d\n", i);
		//print(s);
		/*for(j=+1; j<nearcs[i]; j++)
		{
			if(nearcs[earcs[i][j].u]==0)
			{
				//seconds[i]=j;
				break;
			}
		}*/
	}
	//printg();
	

	values[0]=objval;
	lb[0] = objval;
	sprintf(s, "=%.2f", values[0]);
	print(s);
	cnt=0;
	for(i=0; i<nmodels; i++)
	{
		if(cursol[i]==-1)
		{
			cnt++;
			solutions[0][i]=i;
		}
		else
		{
			solutions[0][i]=earcs[i][cursol[i]].u;
		}
	}
	
	//printf("p0=%d cnt=%d\n", p0, cnt);
	
	//error("hhh");
	

	for(p=p0+1; p<=p1; p++)
	{
		double min=1e30;
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
		int j;
		for( j=0; j<nmodels; j++)
		{
			if(best==j || cursol[j]==-1)
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

		sprintf(s, "=%.2f", values[p-p0]);
		print(s);
		

		if(p == p0+1)
		{
			lb[1] = values[1];
		}
		else
		{
			lb[p-p0] = 0.;
		}

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
		if( lh && p > p0+1)
		{
			sprintf(s, " LH(%d)", p);
			print(s);
			LH(p, high);

			sprintf(s, "=%.2f GAP(%d)=%.2f", lb[p-p0],p, 100.*(1. - lb[p-p0]/values[p-p0]));
			print(s);
			for(int j=0; j<nmodels; j++)
			{
				if(solutions[p-p0][j] == j)
				{
					cursol[j] = -1;
				}
				else
				{
					for(int h=0; h<nearcs[j]; h++)
					{
						int i=earcs[j][h].u;
						if(solutions[p-p0][j] == i)
						{
							cursol[j] = h;
							break;
						}
					}
				}
			}
		}

				
	}

	
	for( i=0; i<nmodels; i++)
	{
		delete[] adj[i];
	}
	delete[] adj;
	delete[] cursol;
	print("\n");
	return 0;
}

int GraphComponent::reopt(int p_)
{
	char s[500];
	sprintf(s, "%5d%12.2f%12.2f%12.2f%10d\n", p_, values[p_-p0], lb[p_-p0], lb[p_-p0]/values[p_-p0]*100, 
		pmed->nfix1a[p_-p0] + pmed->nfix0a[p_-p0]);
	print(s);
	pmed->par2();
	pmed->reopt(p_, solutions[p_-p0]);
	lb[p_-p0] = pmed->blb;
	values[p_-p0] = pmed->bub;
	sprintf(s, "%5d%12.2f%12.2f%12.2f%10d%10.2f\n", p_, values[p_-p0], lb[p_-p0], lb[p_-p0]/values[p_-p0]*100, 
		pmed->nfix0 + pmed->nfix1, pmed->ttime);
	print(s);
	return 0;
}


int GraphComponent::optBB(int p_)
{
	char s[500];
	sprintf(s, "%5d%12.2f%12.2f%12.2f%10d\n", p_, values[p_-p0], lb[p_-p0], lb[p_-p0]/values[p_-p0]*100, 
		pmed->nfix1a[p_-p0] + pmed->nfix0a[p_-p0]);
	print(s);
	pmed->par3();

	//pmed->optBB(p_, solutions[p_-p0]);
	pmed->solveCPX(p_, solutions[p_-p0]);
	
	lb[p_-p0] = pmed->blb;
	values[p_-p0] = pmed->bub;
	sprintf(s, "%5d%12.2f%12.2f%12.2f%10d%10.2f\n", p_, values[p_-p0], lb[p_-p0], lb[p_-p0]/values[p_-p0]*100, 
		pmed->nfix0 + pmed->nfix1, pmed->ttime);
	print(s);
	return 0;
}

int GraphComponent::optimize()
{
	String s;
	sprintf(s, "%-30s\n", "Looking for solution by new version..");
	sprintf(s, "C%4.4d %3d-%3d: init\n", cind, p0, p1);
	print(s);
	//arcsorting();
	
	
	if(nmodels>2000)
	{
		print(" It is big, be patient.\n");
	}
	solutions=new int*[p1-p0+1];
	int p;
	for( p=p0; p<=p1; p++)
	{
		solutions[p-p0]=new int[nmodels];
	}
	values=new double[p1-p0+1];
	lb=new double[p1-p0+1];

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


	for( i=0; i<nmodels; i++)
	{
		for(int j=0; j<nearcs[i]; j++)
		{
			adj[earcs[i][j].u][i/size]+= (1<<(i%size));
		}
	}
	
	int* cursol=new int[nmodels];
	
	for( i=0; i<nmodels; i++)
	{
		cursol[i]=-1;
	}

	double objval=0.;
	int cnt =0;
	for(i=0; i<nmodels; i++)
	{
		if(nearcs[i]==0)
		{
			cnt++;
			continue;
		}
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
	lb[0] = objval;
	sprintf(s, "%6d%12.2f\n", p0, values[0]);
	print(s);
	cnt=0;
	for(i=0; i<nmodels; i++)
	{
		if(cursol[i]==-1)
		{
			cnt++;
			solutions[0][i]=i;
		}
		else
		{
			solutions[0][i]=earcs[i][cursol[i]].u;
		}
	}
	
	//pmed->lang_init();
	pmed->lang_set1();
	
	for(p=p0+1; p<=p1; p++)
	{
		double min=1e30;
		double curval=0;
		int best=-1;
		sprintf(s, "%6d", p);
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
		int j;
		for( j=0; j<nmodels; j++)
		{
			if(best==j || cursol[j]==-1)
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

		//sprintf(s, "=%.2f", values[p-p0]);
		//print(s);
		

		if(p == p0+1)
		{
			lb[1] = values[1];
		}
		else
		{
			lb[p-p0] = 0.;
		}

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
		
		if(p<= p0+1)
		{
			sprintf(s, "%12.2f\n", values[p-p0]);
			print(s);
			continue;
		}
		if(p >= nmodels)
		{
			values[p-p0] = lb[p-p0] = 0.;
			for(i=0; i<nmodels; i++)
			{
				solutions[p-p0][i] = i;
			}
			sprintf(s, "%12.2f\n", values[p-p0]);
			print(s);
			continue;
		}
		
		if(p > p0 + 2)
			pmed->lang_set2();
		pmed->optimize(p, values[p-p0], solutions[p-p0]);
		//pmed.optimize(p, values[p-p0], solutions[p-p0]);
		lb[p-p0] = pmed->blb;
		values[p-p0] = pmed->bub;
		//error("");
		//LH(p, high);

		sprintf(s, "%12.2f%12.2f%12.2f%10d%10.2f\n", values[p-p0], lb[p-p0], lb[p-p0]/values[p-p0]*100, 
			pmed->nfix0 + pmed->nfix1, pmed->ttime);
		print(s);
		for(int j=0; j<nmodels; j++)
		{
			if(solutions[p-p0][j] == j)
			{
				cursol[j] = -1;
			}
			else
			{
				for(int h=0; h<nearcs[j]; h++)
				{
					int i=earcs[j][h].u;
					if(solutions[p-p0][j] == i)
					{
						cursol[j] = h;
						break;
					}
				}
			}
		}
	}

	//pmed->lang_final();
	
	for( i=0; i<nmodels; i++)
	{
		delete[] adj[i];
	}
	delete[] adj;
	delete[] cursol;
	return 0;
}

int GraphComponent::funLH(int p, double* lambda, double* grad, Ro* ro, 
		double& lb, double& ub, int* sol)
{
	// initialize the node redused costs
	int i;
	for( i=0; i<nmodels; i++)
	{
		ro[i].i = i;
		ro[i].fix = (nearcs[i] == 0);
		ro[i].value = -lambda[i];
		ro[i].sol = 0;
	}


	// compute the node reduced costs
	lb=.0;
	int j;
	for( j=0; j<nmodels; j++)
	{
		for(int h=0; h<nearcs[j]; h++)
		{
			int i=earcs[j][h].u;
			if(earcs[j][h].cost-lambda[j] >= 0.)
				break;
			ro[i].value += (earcs[j][h].cost-lambda[j]);
			//if(i==0 || i==2 || i==6)
			{
				//printf("i=%d j=%d ff=%f\n", i, j, earcs[j][h].cost-lambda[j]);
			}
		}
		lb+=lambda[j];
	}
	// sort the reduced cost
	qsort(ro,(size_t) nmodels, sizeof(Ro), compareRoAZ);

	/*printf("lambda: ");
	for(i=0; i<nmodels; i++)
	{
		printf(" %f", lambda[i]);
	}
	printf("\n");

	printf("ro: ");
	for(i=0; i<nmodels; i++)
	{
		printf(" %d-%d-%f", ro[i].i, ro[i].fix, ro[i].value);
	}
	printf("\n");
*/
	//compute the lb
	for(i=0; i<p; i++)
	{
		lb += ro[i].value;
		ro[ro[i].i].sol = 1;
	}

	// compute the subgradient
	int k;
	for( k=0; k<nmodels; k++)
	{
		int j=ro[k].i;
		grad[j]=1.0-ro[j].sol;
		for(int h=0;h<nearcs[j];h++)
		{
			int i=earcs[j][h].u;
			if(earcs[j][h].cost-lambda[j]>= 0.) 
				break;
			if(ro[i].sol == 1) 
				grad[j]--;

		}
	}

	// compute ub
	double ub1 = 0.;

	for( k=0; k<nmodels; k++)
	{
		int j = ro[k].i;
		if(ro[j].sol == 1) 
			continue;
		for(int h=0;h<nearcs[j];h++)
		{
			int i = earcs[j][h].u;
			if(ro[i].sol == 1)
			{
				ub1 += earcs[j][h].cost;
				break;
			}
		}
	}

	if(ub1<0)
		error("fdfsdf");
	if(ub1 < ub )
	{
		//printf("sol:");
		for(int i=0; i< p; i++)
		{
			sol[i] = ro[i].i;
		//	printf(" %d", sol[i]);
		}
		//printf("\n");
		ub = ub1;
	}

	//printf("lb = %f ub = %f ub1 = %f\n", lb, ub, ub1);

	//printf("sol:");
	//for(i=0; i<p; i++)
	//{
	//	printf(" %d", sol[i]);
	//}
	//printf("\n");

	return 0;
}

int GraphComponent::LH(int p, int high)
{
	// The lagrangean multipliers
	double* lambda = new double[nmodels];
	// The subgradient
	double* grad = new double[nmodels];
	// The node reduced costs
	Ro* ro = new Ro[nmodels];
	// The solution
	int* sol = new int[p];

	int cnt = 0;
	int i;
	for(i=0; i<nmodels; i++)
	{
		lambda[i] = 1.;
		if(solutions[p-p0][i] == i )
		{
			//if(cnt==p)
			//	error("cnt=pp");
			sol[cnt++] = solutions[p-p0][i];
		}
	}
	//printf("p=%d cnt=%d\n", p, cnt);
	if(p!=cnt)
		error("Be-be");

	// the lower bound
	double lb = 0.;

	// 
	//CPUTIMER* timer=new CPUTIMER;

	//timer->start();
	String sstop;
	double phi=2.;

	double phid = 1.04;

	switch(high)
	{
	case 1:
		phid = 1.01;
		break;
	case 2:
		phid = 1.005;
	}
	
	int istop;
	int maxitn=1000000;

	double epsg=0.005;
	double epsgap=0.01;
	double epsh=0.001;

	//values[p-p0] = 100000000;

	funLH(p, lambda, grad, ro, lb, values[p-p0], sol);

	//String s;
	//sprintf(s, "%5d%15.2f%15.6f\n", 0, lb, phi);
	//print(s);

	double h;

	double dg;
	//main circle
	int itn;
	{for( itn=1;itn<=maxitn;itn++)
	{
		

		phi/=phid;
		dg=0;
		{for(int i=0;i<nmodels;i++)
		{
			dg+=grad[i]*grad[i];
		}}
		istop=2;
		if(sqrt(dg)<epsg)
			break;
		//dg=sqrt(dg);

		h=phi*(1.05*values[p-p0]-lb)/dg;
		{for(int i=0;i<nmodels;i++)
		{
			lambda[i]+=h*grad[i];
		}}

		funLH(p, lambda, grad, ro, lb, values[p-p0], sol);

		if(itn%1==0)
		{
			//sprintf(s, "%5d%15.2f%15.6f%15.6f\n", itn, lb, values[p-p0], h);
			//print(s);
		}

		if(lb>values[p-p0]-epsgap)
		{
			lb = values[p-p0];
			istop=6;
			break;
		}

		//printf("h=%f dg=%f\n", h, dg);
		istop=3;
		if(h<epsh)
			break;
	}
	if(itn>maxitn)
	{
		istop=4;
		itn--;
	}

	switch(istop)
	{
	case 2:
		sprintf(sstop,"%s", "epsg");
		break;
	case 3:
		sprintf(sstop,"%s", "epsh");
		break;
	case 4:
		sprintf(sstop,"%s", "maxint");
		break;
	case 6:
		sprintf(sstop,"%s", "epsgap");
		break;
	default:
		sprintf(sstop,"%s", "unknown reason");
	}}

	//sprintf(s, "%5d %15.2f %15.6f\n", 0, lb, h);
	//print(s);

	//sprintf(s, "lb=%15.2f\n", lb);
	//print(s);
	

	//print("Stopped by ", sstop, "\n");

	//timer->stop();

	//print("time ");
	//timer->printTimer();
	//print("\n");

	// retrieve the solution

	//printf("sol:");
	//for(i=0; i<nmodels; i++)
	//{
	//	printf(" %d", solutions[p-p0][i]);
	//}
	//printf("\n");
	for(i=0; i<nmodels; i++)
	{
		solutions[p-p0][i] = -1;
	}
	for(i=0; i<p; i++)
	{
		solutions[p-p0][sol[i]] = sol[i];
	}

	for(int j=0; j<nmodels; j++)
	{
		if(solutions[p-p0][j] == j)
			continue;
		for(int h=0; h<nearcs[j]; h++)
		{
			int i = earcs[j][h].u;
			if(solutions[p-p0][i] == i)
			{
				solutions[p-p0][j] = i;
				break;
			}
		}
	}

	//printf("sol:");
	//for(i=0; i<nmodels; i++)
	//{
	//	printf(" %d", solutions[p-p0][i]);
	//}
	//printf("\n");

	GraphComponent::lb[p-p0] = lb;

	delete[] lambda;
	delete[] grad;
	delete[] ro;
	delete[] sol;
	//error("Break point in GraphComponent::LH()");
	return 0;
}

int GraphComponent::packprob(Buffer& buf)
{
	buf.size = 0;
	buf.size = 5*sizeof(int) + 2*nmodels*sizeof(nmodels) + narcs*sizeof(Arc1);
	buf.data = new char[buf.size];
	unsigned int cnt=0;

	memcpy(buf.data+cnt, &cind, sizeof(int));
	cnt +=sizeof(int);
	memcpy(buf.data+cnt, &nmodels, sizeof(int));
	cnt +=sizeof(int);
	memcpy(buf.data+cnt, &narcs, sizeof(int));
	cnt +=sizeof(int);
	memcpy(buf.data+cnt, &p0, sizeof(int));
	cnt +=sizeof(int);
	memcpy(buf.data+cnt, &p1, sizeof(int));
	cnt +=sizeof(int);

	memcpy(buf.data+cnt, models, nmodels*sizeof(int));
	cnt += nmodels*sizeof(int);
	memcpy(buf.data+cnt, nearcs, nmodels*sizeof(int));
	cnt += nmodels*sizeof(int);

	for(int i=0; i<nmodels; i++)
	{
		memcpy(buf.data+cnt, earcs[i], nearcs[i]*sizeof(Arc1));
		cnt += nearcs[i]*sizeof(Arc1);
	}

	solsize =(p1-p0+1)*2*sizeof(double) + (p1-p0+1)*(nmodels)*sizeof(int);

	return 0;
}

GraphComponent::GraphComponent(Buffer& buf)
{
	prob = NULL;
	unsigned int cnt =0;

	memcpy(&cind, buf.data+cnt, sizeof(int));
	cnt +=sizeof(int);
	memcpy(&nmodels, buf.data+cnt, sizeof(int));
	cnt +=sizeof(int);
	memcpy(&narcs, buf.data+cnt, sizeof(int));
	cnt +=sizeof(int);
	memcpy(&p0, buf.data+cnt, sizeof(int));
	cnt +=sizeof(int);
	memcpy(&p1, buf.data+cnt, sizeof(int));
	cnt +=sizeof(int);


	models = new int[nmodels];
	nearcs = new int[nmodels];

	memcpy(models, buf.data+cnt, nmodels*sizeof(int));
	cnt += nmodels*sizeof(int);
	memcpy(nearcs, buf.data+cnt, nmodels*sizeof(int));
	cnt += nmodels*sizeof(int);

	earcs = new Arc1*[nmodels];
	for(int i=0; i<nmodels; i++)
	{
		earcs[i] = new Arc1[nearcs[i]];
		memcpy(earcs[i], buf.data+cnt, nearcs[i]*sizeof(Arc1));
		cnt += nearcs[i]*sizeof(Arc1);
	}

	

	solutions=NULL;
	values=NULL;
	lb = NULL;


}

int GraphComponent::printg()
{
	print("\n$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	String s;
	sprintf(s, "cind=%d nmodels=%d narcs=%d\n", cind, nmodels, narcs);
	print(s);
	sprintf(s, "p0=%d p1=%d\n", p0, p1);
	print(s);
	for(int i=0; i<nmodels; i++)
	{
		sprintf(s, "%d(%d):", i, models[i]);
		print(s);
		for(int h=0; h<nearcs[i]; h++)
		{
			sprintf(s, " %d[%.2f]", earcs[i][h].u, earcs[i][h].cost);
			print(s);
		}
		print("\n");
	}
	print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
	return 0;
}

int GraphComponent::packsol(Buffer &buf)
{
	buf.size = (p1-p0+1)*2*sizeof(double) + (p1-p0+1)*(nmodels)*sizeof(int);
	buf.data = new char[buf.size];
	unsigned int cnt = 0;
	memcpy(buf.data+cnt, values, (p1-p0+1)*sizeof(double));
	cnt += (p1-p0+1)*sizeof(double);
	memcpy(buf.data+cnt, lb, (p1-p0+1)*sizeof(double));
	cnt += (p1-p0+1)*sizeof(double);
	for(int p=p0; p<=p1; p++)
	{
		memcpy(buf.data+cnt, solutions[p-p0], nmodels*sizeof(int));
		cnt += nmodels*sizeof(int);
	}

	/*String s;
	sprintf(s, "Solution of %d packed %d bytes of %d\n", cind, cnt, buf.size);
	print(s);
	for(p=p0; p<=p1; p++)
	{
		sprintf(s, "%d: f=%.2f lb=%.2f", p, values[p-p0], lb[p-p0]);
		print(s);
		for(int i=0; i<nmodels; i++)
		{
			sprintf(s, " %d", solutions[p-p0][i]);
			print(s);
		}
		print("\n");
	}*/
	return 0;
}

int GraphComponent::unpacksol(Buffer &buf)
{
	values = new double[p1-p0+1];
	lb = new double[p1-p0+1];
	solutions = new int*[p1-p0+1];
	unsigned int cnt = 0;
	memcpy(values, buf.data+cnt, (p1-p0+1)*sizeof(double));
	cnt += (p1-p0+1)*sizeof(double);
	memcpy(lb, buf.data+cnt, (p1-p0+1)*sizeof(double));
	cnt += (p1-p0+1)*sizeof(double);
	for(int p=p0; p<=p1; p++)
	{
		solutions[p-p0] = new int[nmodels];
		memcpy(solutions[p-p0], buf.data+cnt, nmodels*sizeof(int));
		cnt += nmodels*sizeof(int);
	}
	/*String s;
	sprintf(s, "Solution of %d size=%d\n", cind, cnt);
	print(s);
	for(p=p0; p<=p1; p++)
	{
		sprintf(s, "%d: f=%.2f lb=%.2f", p, values[p-p0], lb[p-p0]);
		print(s);
		for(int i=0; i<nmodels; i++)
		{
			sprintf(s, " %d", solutions[p-p0][i]);
			print(s);
		}
		print("\n");
	}*/	
	return 0;
}


