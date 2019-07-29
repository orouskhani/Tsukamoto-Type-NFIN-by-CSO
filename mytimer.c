#include <sys/time.h>
#include <sys/types.h>
#include <sys/resource.h>
#include "mytimer.h"
#include "stdlib.h"

mytimer *mytimer_create(void)
{
	mytimer *pTimer;
	
	pTimer=(mytimer *)malloc(sizeof(mytimer));
	pTimer->abs_stTime=0;
	pTimer->state=MYTIMER_STATE_STOP;
	
	return pTimer;
}

double mytimer_get_abs(void)
{
	struct timeval tp;
	struct timezone tz;
    int ret;

	ret=gettimeofday(&tp,&tz);
	if(ret!=0)
		return -1.0;

	return (double)tp.tv_sec+(double)tp.tv_usec/1e6;
}

int mytimer_start(mytimer *pTimer)
{
	double abs_time;

	abs_time=mytimer_get_abs();
	if(abs_time<0.0)
		return -1;

	pTimer->abs_stTime=abs_time;

	return 0;
}

double mytimer_hold(mytimer *pTimer,int index)
{
	double abs_time;

	abs_time=mytimer_get_abs();
	if(abs_time<0.0)
		return -1;

	pTimer->abs_hldTime[index]=abs_time;

	return 0;
}

double mytimer_stop(mytimer *pTimer)
{
	double abs_time;

	abs_time=mytimer_get_abs();
	if(abs_time<0.0)
		return -1;

	pTimer->abs_curTime=abs_time;

	return 0;
}

double mytimer_get_hold(mytimer *pTimer,int index)
{
	return pTimer->abs_hldTime[index]-pTimer->abs_stTime;
}

double mytimer_get_curr(mytimer *pTimer)
{
	return pTimer->abs_curTime-pTimer->abs_stTime;
}

void mytimer_free(mytimer *pTimer)
{
  	free(pTimer);
  	return;
}


