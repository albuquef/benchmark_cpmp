//
// Copyright (c) 2003 Igor Vasil'ev (igor@diima.unisa.it)
// All rights reserved.
//
// logerr.cpp: implementation of fuctions in logerr.h
//


#include "logerr.h"
#include <string.h>
#include <stdlib.h>




String LOGFILE=".log";

void setlogfile(String filename)
{
#ifndef _PARALLEL_
	strcpy(LOGFILE,filename);

	std::ofstream to(LOGFILE);
	to.close();
#endif
}

void print(const char* s1, const char* s2, const char* s3)
{
	printf("%s%s%s", s1, s2, s3);
	/*fflush(stdout);
	std::cout << s1 << s2 << s3 << std::flush;
	fflush(stdout);
	std::ofstream to(LOGFILE,std::ios::app);
	to << s1 << s2 << s3;
	to.close();*/
}

void error(const char* s1, const char* s2, const char* s3)
{
	print("\nError:\n");
	print(s1,s2,s3);
	print("\nProgram is terminated\n");
	exit(1);
}

int compareIntArrayAZ( const void *arg1, const void *arg2 )
{
   if( (* ( int* ) arg1) > * ( int* ) arg2 )
	   return 1;
   else
	   return -1;
}
