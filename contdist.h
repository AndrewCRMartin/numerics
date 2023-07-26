/*===================================================================
 contdist.h

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#ifndef _CONTDIST_H_
#define _CONTDIST_H_

#include "numerics.h"

NUMERICS_EXPORT double betap(double x, double a, double b);
/*-------------------------------------------------------------------
 Returns the value of the Beta random variable distribution
 function with alpha a and beta b at the value x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double betav(double p, double a, double b);
/*-------------------------------------------------------------------
 Returns the critical point of the Beta random variable
 distribution function with alpha a and beta b associated with
 probability p.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double cauchyp(double x, double m, double s);
/*-------------------------------------------------------------------
 Returns the value of the Cauchy random variable distribution
 function with median m and scale parameter s at the value x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double cauchyv(double p, double t, double s);
/*-------------------------------------------------------------------
 Returns the critical point of the Cauchy random variable
 distribution function with median m and scale parameter s
 associated with probability p.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double dblexpp(double x, double m, double s);
/*-------------------------------------------------------------------
 Returns the value of the Double Exponential random variable
 distribution function with mean m and scale parameter s at the
 value x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double dblexpv(double p, double m, double s);
/*-------------------------------------------------------------------
 Returns the critical point of the Double Exponential random
 variable distribution function with mean m and scale parameter s
 associated with probability p.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double expp(double x, double b);
/*-------------------------------------------------------------------
 Returns the value of the Exponential random variable distribution
 function with mean b at the value x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double expv(double p, double b);
/*-------------------------------------------------------------------
 Returns the critical point of the Exponential random variable
 distribution function with mean b associated with probability p.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double gammap(double x, double a, double b);
/*-------------------------------------------------------------------
 Returns the value of the Gamma random variable distribution
 function with alpha a and beta b at the value x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double gammav(double p, double a, double b);
/*-------------------------------------------------------------------
 Returns the critical point of the Gamma random variable
 distribution function with alpha a and beta b associated with
 probability p.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double logisticp(double x, double m, double b);
/*-------------------------------------------------------------------
 Returns the value of the Logistic random variable distribution
 function with mean m and beta b at the value x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double logisticv(double p, double m, double b);
/*-------------------------------------------------------------------
 Returns the critical point of the Logistic random variable
 distribution function with mean m and beta b associated with
 probability p.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double lognormp(double x, double m, double v);
/*-------------------------------------------------------------------
 Returns the value of the Lognormal random variable distribution
 function with mean m and scale parameter v at the value x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double lognormv(double p, double m, double v);
/*-------------------------------------------------------------------
 Returns the critical point of the Lognormal random variable
 distribution function with mean m and scale parameter v associated
 with probability p.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double rayleighp(double x, double s);
/*-------------------------------------------------------------------
 Returns the value of the Rayleigh random variable distribution
 function with scale parameter s at the value x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double rayleighv(double p, double s);
/*-------------------------------------------------------------------
 Returns the critical point of the Rayleigh random variable
 distribution function with scale parameter s associated with
 probability p.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double uniformp(double u, double a, double b);
/*-------------------------------------------------------------------
 Returns the value of the Uniform random variable distribution
 function on the interval [a,b] at the value u.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double uniformv(double p, double a, double b);
/*-------------------------------------------------------------------
 Returns the critical point of the Uniform random variable
 distribution function on the interval [a,b] associated with
 probability p.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double weibullp(double x, double g, double b);
/*-------------------------------------------------------------------
 Returns the value of the Weibull random variable distribution
 function with gamma g and beta b at the value x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double weibullv(double p, double g, double b);
/*-------------------------------------------------------------------
 Returns the critical point of the Weibull random variable
 distribution function with gamma g and beta b associated with
 probability p.
-------------------------------------------------------------------*/

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/

