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
#include <memory.h>
#include <string>

#ifndef _MSC_VER
#include <chrono>
#else
#include <time.h>
#endif

#define INFTY 1.e30




inline double CPU_time()
{
#ifndef _MSC_VER
  using namespace std::chrono;
  auto now = time_point_cast<milliseconds>(system_clock::now());
  using sys_milliseconds = decltype(now);
  return now.time_since_epoch().count()/1000.0;
#else
	return (double)clock()/(double)CLOCKS_PER_SEC;
#endif
}

inline void error(const char* s1, const char *s2 = "")
{
	printf("\nERROR:\n%s %s\n", s1, s2);
	exit(1);
}

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

extern int KK;



