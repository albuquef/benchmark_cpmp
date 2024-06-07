
#include "Convert.h"
#include "ODMProblem.h"
#include "PARAMFILE.h"

double CPU_time();
int nprob1 = 0;
int nprob2 = 0;
int nprob3 = 0;

int main(int argc, char* argv[])
{

	int _nc_ = 1;
#ifndef _OPENMP
	printf("OpenMP is not supported!\n"); 
#endif 

#ifdef _OPENMP
	_nc_ = omp_get_max_threads();
	printf("Max number of threads = %d\n", _nc_);
	if(argc >= 4)
		_nc_ = atoi(argv[3]);
	printf("Number of threads to use = %d\n", _nc_);
	omp_set_num_threads(_nc_);
#endif

	int dp, i;
	int p1 = atoi(argv[2]);
	double beg_t = CPU_time();
	ODMProblem* prob= new ODMProblem(argv[1], p1);

	double t0 = CPU_time();
	prob->optimize();

	double t1 = CPU_time();
	prob->reopt();

	double t2 = CPU_time();
	prob->optBB();

	double t3 = CPU_time();

	double end_t = CPU_time();

	char ss[500];
	sprintf(ss, "Init time = %8.2f\n", t0 - beg_t);
	print(ss);
	sprintf(ss, "Phase 1: #subproblems = %6d time = %8.2f\n", nprob1, t1 - t0);
	print(ss);
	sprintf(ss, "Phase 2: #subproblems = %6d time = %8.2f\n", nprob2, t2 - t1);
	print(ss);
	sprintf(ss, "Phase 3: #subproblems = %6d time = %8.2f\n", nprob3, t3 - t2);
	print(ss);
	sprintf(ss, "Total time = %8.2f\n", end_t - beg_t);
	print(ss);
	
	printf("*********************************************************\n");
	printf("*********************************************************\n");
	printf("wall time = %.2f\n", end_t-beg_t);
	printf("*********************************************************\n");
	printf("*********************************************************\n");

	if(1)
	{
		FILE * ff = fopen("res.txt", "a");
		fprintf(ff, "%s\t%d\t%d\t%d\t%d\t%d\t", prob->probname, prob->ngraphcomponents, prob->nmodels, prob->narcs, prob->p0, prob->p1);
		fprintf(ff, "%d\t%d\t%d\t", nprob1, nprob2, nprob3);
		fprintf(ff, "%.2f\t%.2f\t%.2f\t%.2f\t%d\n", t1 - t0, t2 - t1, t3 - t2, end_t-beg_t, _nc_);
		fclose(ff);
	}

	// saving solution
	//prob->writesol(convert->probname, 2);
	//prob->writesol1("ww.txt");

	// saving graph
	//prob->writegraph1();
	

	//deleting prob
	delete prob;

	//deleting convert
	//delete convert;


	return 0;
}


int main2(int argc, char* argv[])
{


	int dp, i;
	int p1 = atoi(argv[2]);
	ODMProblem* prob= new ODMProblem(argv[1], p1);

	//optimizing
	int method = atoi(argv[3]);//par.intParam("method");


	time_t beg_t = time(NULL);
	
	switch(method)
	{
	case 0:
		break;
	case 1:
		prob->optimize(0,0);
		break;
	case 2:
		prob->optimize(1,0);
		break;
	case 3:
		prob->optimize(1, 1);
		break;
	case 4:
		prob->optimize(1, 2);
		break;
	case 5:
		prob->optimize(1, 3);
		break;
	default:
		print("Unknown method!");
		break;
	}

	time_t end_t = time(NULL);
	printf("*********************************************************\n");
	printf("*********************************************************\n");
	printf("wall time = %d\n", end_t-beg_t);
	printf("*********************************************************\n");
	printf("*********************************************************\n");

	

	// saving solution
	//prob->writesol(convert->probname, 2);
	prob->writesol1("ww.txt");

	// saving graph
	//prob->writegraph1();
	

	//deleting prob
	delete prob;

	//deleting convert
	//delete convert;


	return 0;
}

int main1(int argc, char* argv[])
{

	// Cheking program parameters
	/*if(argc!=2)
	{
		print("\nConverting files...\n");
		error("Enter exactly four program paramenters: parameter file name,\n", 
			"Log is writen in .log.");
	}*/

	PARAMFILE par(argv[1]);

	String fn1, fn2, fn3, fn4, fn5, fn6;

	par.strParam("probname", fn1);
	par.strParam("namefile", fn2);
	par.strParam("costfile", fn3);
	par.strParam("modelfile", fn4);
	par.strParam("modelnamefile", fn5);
	par.strParam("constrfile", fn6);

	int dp = par.intParam("dp");

	//Initialize object convert, where argv[1] is the problem name
	Convert* convert=new Convert(fn1, fn2, fn3, fn4, fn5, fn6);

 
	// processing
	convert->process();

	int** solutions=new int*[dp+1];
	int i;
	for( i=0; i<dp+1; i++)
	{
		solutions[i]=new int[convert->nmodels];
	}
	double* values=new double[dp+1];

	//Get data to solve
	ODMProblem* prob= new ODMProblem(convert->probname, convert->nmodels, convert->modelmasks,
									 convert->modelcosts, convert->modeldemands, 
									 convert->realmodelnames, convert->noptions,
									 convert->hardoptions, dp, solutions, values,
									 convert->constructioncost);





	//optimizing
	int method = atoi(argv[2]);//par.intParam("method");


	time_t beg_t = time(NULL);
	
	switch(method)
	{
	case 0:
		break;
	case 1:
		prob->optimize();
		break;
	case 2:
		prob->optimize(1, 0);
		break;
	case 3:
		prob->optimize(1, 1);
		break;
	case 4:
		prob->optimize(1, 2);
		break;
	case 5:
		prob->optimize(1, 3);
		break;
	default:
		print("Unknown method!");
		break;
	}

	time_t end_t = time(NULL);
	printf("*********************************************************\n");
	printf("*********************************************************\n");
	printf("wall time = %d\n", end_t-beg_t);
	printf("*********************************************************\n");
	printf("*********************************************************\n");

	

	// saving solution
	prob->writesol(convert->probname, 2);

	// saving graph
	 //prob->writegraph();
	

	delete[] values;
	for( i=0; i<dp+1; i++)
	{
		delete[] solutions[i];
	}
	delete[] solutions;

	//deleting prob
	delete prob;

	//deleting convert
	delete convert;


	return 0;
}
