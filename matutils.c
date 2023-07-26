/*===================================================================
 matutils.c

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#include <malloc.h>
#include <stddef.h>

#include "matutils.h"

NUMERICS_EXPORT double **dmatrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
    double **m;
    
    m = (double **)malloc((size_t)((nrh - nrl + 1)*sizeof(double*)));
    if(!m) return NULL;
    m -= nrl;
    for(i = nrl; i <= nrh; i++){
        m[i] = (double *)malloc((size_t)((nch - ncl + 1)*sizeof(double)));
        if(!(m[i])){
            free_dmatrix(m, nrl, i-1, ncl);
            return NULL;
        }
        m[i] -= ncl;
    }
    return m;
}

NUMERICS_EXPORT void free_dmatrix(double **m, int nrl, int nrh, int ncl)
{
    int i;
    
    for(i = nrh; i >= nrl; i--) free((char*)(m[i]+ncl));
    free((char*)(m+nrl));
}

NUMERICS_EXPORT int **imatrix(int nrl, int nrh, int ncl, int nch)
{
    int i;
    int **m;
    
    m = (int **)malloc((size_t)((nrh - nrl + 1)*sizeof(int*)));
    if(!m) return NULL;
    m -= nrl;
    for(i = nrl; i <= nrh; i++){
        m[i] = (int *)malloc((size_t)((nch - ncl + 1)*sizeof(int)));
        if(!(m[i])){
            free_imatrix(m, nrl, i-1, ncl);
            return NULL;
        }
        m[i] -= ncl;
    }
    return m;
}

NUMERICS_EXPORT void free_imatrix(int **m, int nrl, int nrh, int ncl)
{
    int i;
    
    for(i = nrh; i >= nrl; i--) free((char*)(m[i]+ncl));
    free((char*)(m+nrl));
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/
