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
 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <string>

//#include <omp.h>

//#include <process.h>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_01.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#define NC_ 32

#define INFTY 1.e30

#include <omp.h>



inline double CPU_time()
{

	double tt = (double)clock()/(double)CLOCKS_PER_SEC;
#ifdef _OPENMP
	tt = omp_get_wtime();
#endif
#ifdef MPI_INCLUDED
	tt = MPI_Wtime();
#endif
	return tt;

}

//inline void error(const char* s1, const char *s2 = "")
//{
//	printf("\nERROR:\n%s %s\n", s1, s2);
//	exit(1);
//}

inline long long atoll_(const char* str)
{
    long long i = 0;
    while(strcmp(str, "") != 0)
    {
       i*=10;
       i += *str++ - '0';
    }
    return i;
}




