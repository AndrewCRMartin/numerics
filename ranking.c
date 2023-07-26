/*===================================================================
 ranking.c

 Version 1.1

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  May 18, 1999
===================================================================*/

#include <math.h>

#include "algorthm.h"
#include "datamanp.h"
#include "descript.h"
#include "numerror.h"
#include "ranking.h"
#include "sort.h"
#include "vecutils.h"

NUMERICS_EXPORT void index(double *first, double* last, int *indx)
{
    int i, n = length(first, last);
    double* cp = dvector(0, n - 1);
    double* ci = dvector(0, n - 1);
    
    copy(first, last, cp);
    for(i = 0; i < n; ++i){
        *(ci + i) = i;
    }
    
    qsort2(cp, cp + n, ci);
    
    for(i = 0; i < n; ++i){
        *(indx + i) = (int)(*(ci + i));
    }
    
    free_dvector(cp, 0);
    free_dvector(ci, 0);
}

NUMERICS_EXPORT void rank(double* first, double* last, double* irank, int ties)
{
    int n = length(first, last);
    int* indx = ivector(0, n - 1);
    int i, j, start, end;
    double rnk;
    
    index(first, last, indx);
     
    j = 0;
    while(j < n){
        start = end = j;
        while(end < n - 1 && *(first + *(indx + j)) == *(first + *(indx + end + 1))){
            ++end;
        }
        if(ties == RANK_MEAN){
            rnk = average(start, end);
        } else if(ties == RANK_LOW){
            rnk = start;
        } else if(ties == RANK_HIGH){
            rnk = end;
        } else {
            NUMERICS_ERROR("rank", "Invalid ties parameter.");
            rnk = start;
        }
        for(i = start; i <= end; ++i){
            *(irank + *(indx + i)) = rnk + 1;
        }
        j = end + 1;
    }
    
    free_ivector(indx, 0);
}

NUMERICS_EXPORT void arank(double* first, double* last, double* irank, int ties)
{
    int i = 0, n = length(first, last);
    double* vec = dvector(0, n - 1);
    
    while(first < last){
        *(vec + i) = fabs(*first);
        ++first;
        ++i;
    }

    rank(vec, vec + n, irank, ties);
    
    free_dvector(vec, 0);
}

NUMERICS_EXPORT void abrank(double* first, double* last, double* irank, int ties)
{
    int n = last - first;
    double m = (n + 1.0) / 2.0;
    int i;
    
    rank(first, last, irank, ties);
    for(i = 0; i < n; ++i){
        if(*irank > m){
            *irank = n - *irank + 1.0;
        }
        ++irank;
    }
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
 Version 1.1 - 05/18/1999 - Added abrank.
===================================================================*/
