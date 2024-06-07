//
// Copyright (c) 2003 Igor Vasil'ev (igor@diima.unisa.it)
// All rights reserved.
//
// CPUTIMER.cpp: implementation of the CPUTIMER class.
//
//////////////////////////////////////////////////////////////////////


#include "CPUTIMER.h"
#include <stdio.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

unsigned long CPUTIMER::centiseconds()
{
	return centiseconds_+getCurTime();
}
CPUTIMER::CPUTIMER(unsigned long beg)
{
	stopped_=true;
	centiseconds_=beg;

}

CPUTIMER::~CPUTIMER()
{

}

unsigned long CPUTIMER::getCurTime()
{
	if(stopped_)
		return 0;
	else
	{
		clock_t end_=clock()*100/ CLOCKS_PER_SEC;
		return end_-beg_;
	}
	
}

void CPUTIMER::start(bool reset)
{
	if(!stopped_)
	{
		error("Starting started timer.");
	}
	stopped_=false;
	beg_=clock()*100/ CLOCKS_PER_SEC;
	if(reset)
		centiseconds_=0;
}

void CPUTIMER::stop()
{
	if(stopped_)
	{
		error("Stopping stopped timer.");
	}
	centiseconds_+=getCurTime();
	stopped_=true;
}

int CPUTIMER::getString(String s, bool ext)
{
	unsigned long d_=days();

	unsigned long h_=hours()-d_*24;
	if(!ext)
	{
		d_=0;
		h_=hours();
	}
	unsigned long m_=minuts()-h_*60-d_*60*24;
	unsigned long s_=seconds()-m_*60-h_*60*60-d_*60*60*24;
	unsigned long c_=centiseconds()-s_*100-m_*60*100-
					h_*60*60*100-d_*60*60*24*100;
	if(ext)
		sprintf(s,"%3.2dd%2.2d:%2.2d:%2.2d.%2.2d",d_,h_,m_,s_,c_);
	else
		sprintf(s,"%6.4d:%2.2d:%2.2d.%2.2d",h_,m_,s_,c_);

	return 0;
}

int CPUTIMER::printTimer(bool ext)
{
	String s;
	getString(s,ext);
	print(s);
	return 0;
}
