/*===================================================================
 algorthm.c

 Version 1.1

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  May 18, 1999
===================================================================*/

#include "algorthm.h"

double average(double x, double y)
{
    return x + (y - x) / 2.0;
}

void pairwdiff(double* first1, double* last1, double* first2, double* last2, double* dest)
{
    double* i;
    double* j;

    for(i = first1; i < last1; ++i){
        for(j = first2; j < last2; ++j){
            *dest = *i - *j;
            ++dest;
        }
    }
}

void walshavg(double* first, double* last, double* dest)
{
    double *i, *j;

    for(i = first; i < last; ++i){
        for(j = i; j < last; ++j){
            *dest = average(*i, *j);
            ++dest;
        }
    }
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
 Version 1.1 - 05/18/1999 - Added pairwdiff.
===================================================================*/
