/*===================================================================
 correlat.c

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare
 
 Created:  August 28, 1998
 Revised:  
===================================================================*/

#include <math.h>

#include "correlat.h"
#include "descript.h"
#include "ranking.h"
#include "vecutils.h"

NUMERICS_EXPORT double covariance(double* first1, double* last1, double* first2)
{
    return sumproduct(first1, last1, first2) / length(first1, last1);
}

NUMERICS_EXPORT double correlation(double* first1, double* last1, double* first2)
{
    return sumproduct(first1, last1, first2) / sqrt(sumsquare(first1, last1) * sumsquare(first2, first2 + (last1 - first1)));
}

NUMERICS_EXPORT double kendall(double* first1, double* last1, double* first2)
{
    int nc, nd, ey, ex;
    int i, j;
    int n = length(first1, last1);
    
    nc = nd = ey = ex = 0;
    
    for(i = 0; i < n; ++i){
        for(j = 0; j < i; ++j){
            if(*(first1 + i) < *(first1 + j)){
                if(*(first2 + i) < *(first2 + j)){
                    ++nc;
                } else if(*(first2 + i) > *(first2 + j)){
                    ++nd;
                } else {
                    ++ex;
                }
            } else if(*(first1 + i) > *(first1 + j)){
                if(*(first2 + i) < *(first2 + j)){
                    ++nd;
                } else if(*(first2 + i) > *(first2 + j)){
                    ++nc;
                } else {
                    ++ex;
                }
            } else {
                ++ey;
                if(*(first2 + i) == *(first2 + j)){
                    ++ex;
                }
            }
        }
    }
    
    return (double)(nc - nd) / sqrt((double)((nc + nd + ey) * (nc + nd + ex)));
}

NUMERICS_EXPORT double spearman(double* first1, double* last1, double* first2)
{
    int n = length(first1, last1);
    double *r1 = dvector(0, n - 1);
    double *r2 = dvector(0, n - 1);
    double ret;
    
    rank(first1, last1, r1, RANK_MEAN);
    rank(first2, first2 + n, r2, RANK_MEAN);
    
    ret = correlation(r1, r1 + n, r2);
    
    free_dvector(r1, 0);
    free_dvector(r2, 0);
    
    return ret;
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/

