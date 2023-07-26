/*===================================================================
 residual.c

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

#include "residual.h"

NUMERICS_EXPORT void resid(double* first, double* last, double val)
{
    while(first != last){
        *first -= val;
        ++first;
    }
}

NUMERICS_EXPORT void aresid(double* first, double* last, double val)
{
    while(first != last){
        *first = fabs(*first - val);
        ++first;
    }
}

NUMERICS_EXPORT void nresid(double* first, double* last, double val, double n)
{
    if(n != 1.0){
        while(first != last){
            *first = pow(*first - val, n);
            ++first;
        }
    } else {
        resid(first, last, val);
    }
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/
