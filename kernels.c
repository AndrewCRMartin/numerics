/*===================================================================
 kernels.cpp

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

#include "kernels.h"

double biweight(double t)
{
    double tmp;
    
    if(-1.0 <= t && t <= 1.0){
        tmp = (1.0 - t * t);
        return .3125 * tmp * tmp;
    }
    
    return 0.0;
}

double cosine(double x)
{
    if(-1.0 <= x && x <= 1.0) return .5 + .5 * cos(NUMERICS_PI * x);
    return 0.0;
}

double epanech(double x)
{
    if(-1.0 <= x && x <= 1.0) return .75 - .75 * x * x;
    return 0.0;
}

double gaussian(double x)
{
    return 0.3989422804014327 * exp(-x*x*.5);
}

double rectangle(double x)
{
    if(-1.0 <= x && x <= 1.0) return 0.5;
    return 0.0;
}

double triangle(double x)
{
    if(-1.0 <= x && x <= 1.0) return fabs(1.0 - x);
    return 0.0;
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/


