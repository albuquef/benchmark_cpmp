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
#include <stdlib.h>




#include "Pmed.h"

int main(int argc, char* argv[])
{
	if(argc <2)
		error("Give at least 1 argument.");

	switch(argv[1][0])
	{
	case 'V':// cheching solution
		{
			Pmed pmed;
			try
			{
				pmed.readprob(argv[2], false);
				pmed.p = atoi(argv[3]);
			}
			catch(...)
			{
				error("Wrong command arguments.");
			}
			pmed.ver();
			break;
		}
	case 'L':// only subgradient
		{
			Pmed pmed;
			double tt = CPU_time();
			try
			{
				pmed.readprob(argv[2],false);
				pmed.p = atoi(argv[3]);
				pmed.bub = atof(argv[4]);
				pmed.MEMSIZE = atoll_(argv[5]);
			}
			catch(...)
			{
				error("Wrong command arguments.");
			}
			pmed.lang_set1();
			pmed.PRIMALFREQ = 0;
			pmed.writedist();
			//pmed.readdist();
			double tt1 = CPU_time();
			pmed.lang_init();
			pmed.lang();
			pmed.lang_final();

			printf("Total time %.2f\n", CPU_time() - tt);
			printf("Lag time %.2f\n", CPU_time() - tt1);
			FILE *ff = fopen("res.txt", "a");
			fprintf(ff, "%s\t%d\t%d\t%d\t%.2lf\t%.2lf\t%.4f\t%.2f%.2f\n", pmed.probname, 
				pmed.m, pmed.n, pmed.p, pmed.blb, pmed.bub, 
				(pmed.bub-pmed.blb)/pmed.bub * 100., CPU_time() - tt,
				CPU_time() - tt1);
			fclose(ff);
			break;
		}
	case 'C':// Core heuristic
		{
			Pmed pmed;
			double tt = CPU_time();
			try
			{
				pmed.readprob(argv[2], false);
				/* pmed.readprobmatr(argv[2]); */
				pmed.p = atoi(argv[3]);
				pmed.bub = atof(argv[4]);
				pmed.MEMSIZE = atoll_(argv[5]);
				pmed.CORETIME = atoi(argv[8]);
				pmed.MAXNODE = atof(argv[6]);
				pmed.MAXCORE = atoi(argv[7]) * (double) pmed.m;
				pmed.MINARCS = 10;
				pmed.r = atoi(argv[9]);
			}
			catch(...)
			{
				error("Wrong command arguments.");
			}
			double tt1 = CPU_time();
			pmed.lang_set1c();
			pmed.PRIMALFREQ = 1;
			pmed.writedist();
			//pmed.readdist();
			pmed.lang_init();
			pmed.lang();
			pmed.writecore4();
			pmed.lang_set2();
			pmed.lang();
			pmed.writecore4();
			pmed.lang_set3();
			pmed.lang();
			pmed.lang_final();

			printf("Total time %.2f\n", CPU_time() - tt);
			printf("Heur time %.2f\n", CPU_time() - tt1);
			FILE *ff = fopen("res.txt", "a");
			fprintf(ff, "%s\t%d\t%d\t%d\t%d\t%.2lf\t%.2lf\t%.4f\t%.2f\t%.2f\n", pmed.probname, 
				pmed.m, pmed.n, pmed.p, pmed.r, pmed.blb, pmed.bub, 
				(pmed.bub-pmed.blb)/pmed.bub * 100., CPU_time() - tt,
				CPU_time() - tt1);
			fclose(ff);
			pmed.savesol1();
			break;
		}
	case 'A':// Aggregation heuristic
		{
			Pmed pmed;
			int pp;
			double tt = CPU_time();
			try
			{
				pmed.readprob(argv[2], false);
				pp = atoi(argv[3]);
				pmed.bub = atof(argv[4]);
				pmed.MEMSIZE = atoll_(argv[5]);
				pmed.CORETIME = atoi(argv[8]);
				pmed.MAXNODE = atoi(argv[6]);
				pmed.MAXCORE = atoi(argv[7]) * pmed.m;
				pmed.p = atoi(argv[9]);
				pmed.MINARCS = 3;
			}
			catch(...)
			{
				error("Wrong command arguments.");
			}
			
			pmed.lang_set1c();
			pmed.PRIMALFREQ = 10000000;
			pmed.writedist();
			//pmed.readdist();
			pmed.lang_init();
			pmed.lang();
			pmed.writecore4();
			pmed.lang_set2();
			pmed.lang();
			pmed.writecore4();
			pmed.lang_set3();
			pmed.lang();
			
			pmed.savesol(pp);
			pmed.lang_final();

			pmed.lang_set1();
			pmed.PRIMALFREQ = 0;
			
			pmed.lang_init();
			pmed.lang();
			pmed.lang_final();

			printf("Total time %.2f\n", CPU_time() - tt);
			FILE *ff = fopen("res.txt", "a");
			fprintf(ff, "%s\t%d\t%d\t%d\t%.2lf\t%.2lf\t%.4f\t%.2f\t%d\t%d\t%.2f\n", pmed.probname, 
				pmed.m, pmed.n, pmed.p, pmed.blb, pmed.bub, 
				(pmed.bub-pmed.blb)/pmed.bub * 100., CPU_time() - tt,
				pmed.MAXNODE, pmed.MAXCORE, pmed.CORETIME);
			fclose(ff);
			pmed.savesol1();
			break;
		}
	default:
		{
			error("Unknown method lable.");
		}
	}
	return 0;
}


