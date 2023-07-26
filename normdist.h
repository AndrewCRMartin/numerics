/*===================================================================
 normdist.h

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#ifndef _NORMDIST_H_
#define _NORMDIST_H_

#include "numerics.h"

NUMERICS_EXPORT double chisqp(double x2, long v);
/*-------------------------------------------------------------------
 Returns the value of the Chi Square random variable distribution
 function with v degrees of freedom at the value x2.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double chisqv(double p, long v);
/*-------------------------------------------------------------------
 Returns the critical point of the Chi Square random variable
 distribution function with v degrees of freedom associated with
 probability p.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double fratiop(double f, long v1, long v2);
/*-------------------------------------------------------------------
 Returns the value of the F-ratio random variable distribution
 function with v1 degrees of freedom numerator and v2 degrees of
 freedom denomenator at the value f.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double fratiov(double p, long v1, long v2);
/*-------------------------------------------------------------------
 Returns the critical point of the F-ratio random variable
 distribution function with v1 degrees of freedom numerator and v2
 degrees of freedom denomenator associatied with probability p.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double normalp(double x, double m, double v);
/*-------------------------------------------------------------------
 Returns the value of the Normal random variable distribution
 function with mean m and variance v at the value x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double normalv(double p, double m, double v);
/*-------------------------------------------------------------------
 Returns the critical point of the Normal random variable
 distribution function with mean m and variance v associated
 with probability p.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double studtp(double t, double v);
/*-------------------------------------------------------------------
 Returns the value of Student's T random variable distribution
 function with v degrees of freedom at the value x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double studtv(double p, double v);
/*-------------------------------------------------------------------
 Returns the critical point of Student's T random variable
 distribution function with v degrees of freedom associated with
 probability p.
-------------------------------------------------------------------*/

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/
