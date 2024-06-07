// ommom.cpp: implementation of the Commom class.
//
//////////////////////////////////////////////////////////////////////

#include "ommom.h"
#include "time.h"

double CPU_time();

int Commom::DP1()
{
	double tt = CPU_time();
	char ss[500];
	//print("\nKnapsack problem optimize...");
	int w=0, i=0, j=0;
	int k = p1 - n + 1;
	double INF = 1.e30;

	double ** c = new double*[k+1];
	for(w=0; w<=k; w++)
	{
		c[w] = new double[n];
		for(i=0; i<n ;i++)
			c[w][i] = INF;
	}
	for(i=0; i<n; i++)
	{
		for(w=ps[i][0]; w<=ps[i][1]; w++)
		{
			c[w][i]= costs[i][w-ps[i][0]];
		}
	}	

	double ** g = new double*[n];
	for(j=0; j<n; j++)
	{
		g[j] = new double[p1+1];
		for(w=0; w<=p1; w++)
		{
			g[j][w] = INF;
		}
	}
	for(w=ps[0][0]; w<=ps[0][1]; w++)
	{
		g[0][w] = c[w][0];
	}
	for(j=1; j<n; j++)
	{
		for(w=1; w<=k+j-1; w++)
		{
			g[j][w] = INF;
			for(i=1; i<=w-j+1; i++)
			{
				if (g[j][w] > g[j-1][w-i] + c[i][j])
				{
					//printf("", j, w, i);
					//printf("%d  %d  %d\n", j, w, i);
					g[j][w] = g[j-1][w-i] + c[i][j];
				}
			}
		}
	}
	j=n-1;
	w=p1;
	g[j][w] = INF;
	for(i=1; i<=k; i++)
	{
		if (g[j][w] > g[j-1][w-i] + c[i][j])
		{
			g[j][w] = g[j-1][w-i] + c[i][j];
		}
	}

	if(0)
	{
		FILE *ff = fopen("lll.txt","w");
		for(j=0; j<n; j++)
		{
			fprintf(ff, "%4d-%3d-%3d:", j, ps[j][0], ps[j][1]);
			for(w=0; w<=p1; w++)
			{
				if(	g[j][w] >= INF)
					fprintf(ff, "%4d=%8.1f(-)", w, 0.);
				else
					fprintf(ff, "%4d=%8.1f(-)", w, g[j][w]);
			}
			fprintf(ff, "\n");
		}
		fclose(ff);
	}

	int p;
	for(p=p0; p<=p1; p++)
	{
		values[p-p0] = lb[p-p0] = g[n-1][p];
	}

	/*char ss[500];
	for(p=p0; p<=p1; p++)
	{
		sprintf(ss, "f(%d) = %.4f\n", p, values[p-p0]);
		print(ss);
	}*/

	
	for(w=0; w<=k; w++)
		delete[] c[w];
	delete[] c;

	for(j=0; j<n; j++)
		delete[] g[j];
	delete[] g;

	sprintf(ss, "%.2f seconds\n", CPU_time() - tt);
	//print("OK in ", ss);
	return 0;
}

int Commom::DP()
{
	double tt = CPU_time();
	char ss[500];
	//print("\nKnapsack problem optimize...");
	int w=0, i=0, j=0;
	//int k = p1 - n + 1;
	double INF = 1.e30;

	double ** c = new double*[p1+1];
	for(w=0; w<=p1; w++)
	{
		c[w] = new double[n];
		for(i=0; i<n ;i++)
			c[w][i] = INF;
	}
	for(i=0; i<n; i++)
	{
		for(w=ps[i][0]; w<=ps[i][1]; w++)
		{
			c[w][i]= costs[i][w-ps[i][0]];
		}
	}	

	double ** g = new double*[n];
	int ** s = new int*[n];
	for(j=0; j<n; j++)
	{
		s[j] = new int[p1+1];
		g[j] = new double[p1+1];
		for(w=0; w<=p1; w++)
		{
			g[j][w] = INF;
			s[j][w] = -1;
		}
	}
	for(w=ps[0][0]; w<=ps[0][1]; w++)
	{
		g[0][w] = c[w][0];
		s[0][w] = w;
	}
	for(j=1; j<n; j++)
	{
		for(w=1; w<=p1; w++)
		{
			g[j][w] = INF;
			for(i=1; i<=w-j+1; i++)
			{
				if (g[j][w] > g[j-1][w-i] + c[i][j])
				{
					//printf("", j, w, i);
					//printf("%d  %d  %d\n", j, w, i);
					s[j][w] = i;
					g[j][w] = g[j-1][w-i] + c[i][j];
				}
			}
		}
	}
	j=n-1;
	w=p1;
	g[j][w] = INF;
	for(i=1; i<=p1; i++)
	{
		if (g[j][w] > g[j-1][w-i] + c[i][j])
		{
			g[j][w] = g[j-1][w-i] + c[i][j];
			s[j][w] = i;
		}
	}

	if(0)
	{
		FILE *ff = fopen("ll.txt","w");
		for(j=0; j<n; j++)
		{
			fprintf(ff, "%4d-%3d-%3d:", j, ps[j][0], ps[j][1]);
			for(w=0; w<=p1; w++)
			{
				if(	g[j][w] >= INF)
					fprintf(ff, "%4d=%8.1f(%4d)", w, 0., s[j][w]);
				else
					fprintf(ff, "%4d=%8.1f(%4d)", w, g[j][w], s[j][w]);
			}
			fprintf(ff, "\n");
		}
		fclose(ff);
	}

	int p;
	for(p=p0; p<=p1; p++)
	{
		values[p-p0] = lb[p-p0] = g[n-1][p];
		w = p;
		solutions[p-p0][n-1] = s[n-1][w];

		for(j=n-2; j>=0; j--)
		{
			w -= s[j+1][w];
			solutions[p-p0][j] = s[j][w];
		}
	}

	/*char ss[500];
	for(p=p0; p<=p1; p++)
	{
		sprintf(ss, "f(%d) = %.4f\n", p, values[p-p0]);
		print(ss);
	}*/

	
	for(w=0; w<=p1; w++)
		delete[] c[w];
	delete[] c;

	for(j=0; j<n; j++)
		delete[] g[j];
	delete[] g;

	for(j=0; j<n; j++)
		delete[] s[j];
	delete[] s;

	sprintf(ss, "%.2f seconds\n", CPU_time() - tt);
	//print("OK in ", ss);
	return 0;
}

/////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Commom::Commom(int n_, double** costs_, int** ps_, int p0_, int p1_): 
		n(n_), costs(costs_), ps(ps_), p0(p0_), p1(p1_)
{
	values = new double[p1-p0+1];
	lb = new double[p1-p0+1];
	tmp = new int[n];
	solutions = new int*[p1-p0+1];
	for(int p=p0; p<=p1; p++)
	{
		solutions[p-p0] = new int[n];
		lb[p-p0] = .0;
	}
}

Commom::~Commom()
{
	delete[] tmp;
	delete[] lb;
	delete[] values;
	for(int p=p0; p<=p1; p++)
	{
		delete[] solutions[p-p0];
	}
	delete[] solutions;
}

int Commom::printdata()
{
	print("\nKnapsack problem\n*******************************************************\n");
	String s;

	sprintf(s, "n=%d\n", n);
	print(s);
	sprintf(s, "p0=%d\n", p0);
	print(s);
	sprintf(s, "p1=%d\n", p1);
	print(s);
	for(int i=0; i<n; i++)
	{
		sprintf(s, "%3d:%3d-%3d\n", i, ps[i][0], ps[i][1]);
		print(s);
	}
	print("*******************************************************\n");

	return 0;

}

int Commom::optimize()
{
    print("\nKnapsack problem optimize\n-------------------------------------------------------\n");
	String s;


	values[0] = 0.;
	int cnt = 0;
	for(int i=0; i<n; i++)
	{
		values[0] += costs[i][0];
		cnt += ps[i][0];
		solutions[0][i] = ps[i][0];
	}

	sprintf(s, "f(%d)=%.2f\n", p0, values[0]);
	print(s);

	lb[0] = values[0];

	for(int p=p0+1; p<=p1; p++)
	{
		double min = 0;
		int best = -1;
		for(int i=0; i<n; i++)
		{
			if(solutions[p-1-p0][i]!=ps[i][1] && 
				min>costs[i][ solutions[p-1-p0][i] - ps[i][0] +1]-costs[i][ solutions[p-1-p0][i] - ps[i][0] ])
			{
				best = i;
				min = costs[i][ solutions[p-1-p0][i] - ps[i][0] +1]-costs[i][ solutions[p-1-p0][i] - ps[i][0] ];
			}
			solutions[p-p0][i] = solutions[p-p0-1][i];
		}

		values[p-p0] = values[p-p0-1] + min;
		solutions[p-p0][best]++;
		sprintf(s, "f(%d)=%.2f", p, values[p-p0]);
		print(s);

		if(p == p0+1)
		{
			lb[1] = values[1];
		}
		if(p > p0+1)
		{
			LH(p);
			sprintf(s, " lb(%d)=%.2f GAP = %.4f", p, lb[p-p0], (values[p-p0]-lb[p-p0])/values[p-p0]*100.);
			print(s);

		}
		print("\n");
	}


	print("-------------------------------------------------------\n");
	return 0;
}

int Commom::LH(int p)
{
	// The lagrangean multipliers
	double lambda = 0.;
	// The subgradient
	double grad = 0.;
	
	// The solution
	int* sol = new int[n];
	int i;
	for( i=0; i<n; i++)
	{
		sol[i] = solutions[p-p0][i];
	}

	// the lower bound
	double lb = 0.;

	// 
	//CPUTIMER* timer=new CPUTIMER;

	//timer->start();
	String sstop;
	double phi=2.;

	double phid = 1.01;
	
	int istop;
	int maxitn=100000;

	double epsg=0.005;
	double epsgap=0.01;
	double epsh=0.001;

	values[p-p0] = 100000000;

	funLH(p, lambda, grad, lb, values[p-p0], sol);

	String s;
	sprintf(s, "\n%5d%15.2f%15.6f\n", 0, lb, phi);
	//print(s);

	double h;

	double dg;
	//main circle
	int itn;
	{for(itn=1;itn<=maxitn;itn++)
	{
		

		phi/=phid;
		dg = grad*grad;
		istop=2;
		if(sqrt(dg)<epsg)
			break;
		
		h= phi*(1.05*values[p-p0]-lb)/dg;
		lambda += h*grad;
		

		funLH(p, lambda, grad, lb, values[p-p0], sol);

		if(itn%1==0)
		{
			sprintf(s, "%5d%15.2f%15.6f%15.6f %.0f %.2f\n", itn, lb, values[p-p0], h, grad, lambda);
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

	sprintf(s, "%5d %15.2f %15.6f\n", 0, lb, h);
	//print(s);

	//print("Stopped by ", sstop, "\n");

	//timer->stop();

	//print("time ");
	//timer->printTimer();
	//print("\n");

	Commom::lb[p-p0] = lb;
	for(i=0; i<n; i++)
	{
		solutions[p-p0][i] = sol[i];
	}

	
	delete[] sol;
	//error("Break point in GraphComponent::LH()");
	return 0;

}

int Commom::funLH(int p, double lambda, double &grad, double &lb, double &ub, int *sol)
{

	lb = p*lambda;
	grad = p;
	double ub1 = 0.;

	int pp = 0;

	for(int i=0; i<n; i++)
	{
		double min = 1e30;
		int best;
		for(int j = ps[i][0]; j<=ps[i][1]; j++)
		{
			if(min > costs[i][j - ps[i][0]] - j*lambda)
			{
				best = j;
				min = costs[i][j - ps[i][0]] - j*lambda;
			}
		}
		lb += min;
		grad -= best;
		tmp[i] = best;
		pp += best;
		ub1 += costs[i][best-ps[i][0]];
	}

	//error("Break point in funLH");

	if( pp < p )
	{
		for( pp = pp+1; pp<=p; pp++)
		{
			double min = 0;
			int best = -1;
			for(int i=0; i<n; i++)
			{
				if(tmp[i]!=ps[i][1] && min>costs[i][ tmp[i] - ps[i][0] +1]-costs[i][ tmp[i] - ps[i][0] ])
				{
					best = i;
					min = costs[i][ tmp[i] - ps[i][0] +1]-costs[i][ tmp[i] - ps[i][0] ];
				}
			}

			ub1 += min;
			tmp[best]++;
		}
	}

	if(pp > p)
	{
		for( pp = pp-1; pp>=p; pp--)
		{
			double max = 1e30;
			int best = -1;
			for(int i=0; i<n; i++)
			{
				if(tmp[i]!=ps[i][0] && max > costs[i][ tmp[i] - ps[i][0] -1]-costs[i][ tmp[i] - ps[i][0] ])
				{
					best = i;
					max = costs[i][ tmp[i] - ps[i][0] -1]-costs[i][ tmp[i] - ps[i][0] ];
				}
			}

			ub1 += max;
			tmp[best]--;
		}
	}

	if( ub1+0.001 < ub)
	{
		ub = ub1;
		for(int i=0; i<n; i++)
		{
			sol[i] = tmp[i];
		}
	}

	return 0;
}
