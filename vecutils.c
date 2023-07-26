/*===================================================================
 vecutils.c

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

#include "vecutils.h"

NUMERICS_EXPORT double *dvector(int nl, int nh)
{
    double *v;
    
    v = (double *)malloc((size_t)((nh - nl + 1)*sizeof(double)));
    if(!v) return NULL;
    return v - nl;
}

NUMERICS_EXPORT void free_dvector(double *v, int nl)
{
    free((char *)(v + nl));
}

NUMERICS_EXPORT int *ivector(int nl, int nh)
{
    int *v;
    
    v = (int *)malloc((size_t)((nh - nl + 1)*sizeof(int)));
    if(!v) return NULL;
    return v - nl;
}

NUMERICS_EXPORT void free_ivector(int *v, int nl)
{
    free((char *)(v + nl));
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/98 - New.
===================================================================*/

