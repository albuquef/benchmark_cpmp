//
// Copyright (c) 2003 Igor Vasil'ev (igor@diima.unisa.it)
// All rights reserved.
//
// logerr.h: declaration for logging and error and some usefull macros and functions

#ifndef LOGERROR_H
#define LOGERROR_H

//#include <algorithm>

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#ifdef _MSC_VER
#else
	using namespace std;
#endif


// MIN and MAX macroses
#define MIN(a,b)               ((a) < (b) ? (a) : (b))
#define MAX(a,b)               ((a) > (b) ? (a) : (b))



//type of string
typedef char String[300];


// procedure is set log file with filename
void setlogfile(String filename);

// print s1, s2 and s3 in the standard output and log file
void print(const char* s1, const char* s2="", const char* s3="");

// print s1, s2 and s3 in the standard output and log file
// with error message and exit the programm with exit code 1
void error( const char* s1, const char* s2="", const char* s3="");

// *arg1 and arg2 are integer, if *arg1 > *arg2 it returns 1, else -1.
// this function is used to order int array in qsort.
int compareIntArrayAZ( const void *arg1, const void *arg2 );

#endif  LOGERROR_H
