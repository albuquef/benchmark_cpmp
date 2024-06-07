//
// Copyright (c) 2003 Igor Vasil'ev (igor@diima.unisa.it)
// All rights reserved.
//
// CPUTIMER.h: interface for the CPUTIMER class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_CPUTIMER_H__645CB0CE_5688_45CD_98C1_AB17DC97C7E7__INCLUDED_)
#define AFX_CPUTIMER_H__645CB0CE_5688_45CD_98C1_AB17DC97C7E7__INCLUDED_


#include "time.h"
#include "logerr.h"

//class to time the CPU time
class CPUTIMER  
{
public:
	// print the time, if ext=true, with number of days
	int printTimer(bool ext=false);
	// return the formated time, if ext=true, with number of days
	int getString(String s, bool ext=false);
	//stop timer
	void stop();
	// start timer, if reset then start from 0
	void start(bool reset=false);
	//number of centiseconds
	unsigned long centiseconds();
	//number of seconds
	unsigned long seconds(){return centiseconds()/100;};
	//number of minuts
	unsigned long minuts(){return seconds()/60;};
	//number of hours
	unsigned long hours(){return minuts()/60;};
	//number of days
	unsigned long days(){return hours()/24;};

	//constructor, beg is a initial number of centiseconds
	CPUTIMER(unsigned long beg=0);
	//destructor
	virtual ~CPUTIMER();
private:
	bool stopped_;
	unsigned long centiseconds_;
	clock_t beg_;
	


protected:
	virtual unsigned long getCurTime();
};

#endif // !defined(AFX_CPUTIMER_H__645CB0CE_5688_45CD_98C1_AB17DC97C7E7__INCLUDED_)
