/*===================================================================
 discdist.h

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#ifndef _DISCDIST_H_
#define _DISCDIST_H_

#include "numerics.h"

NUMERICS_EXPORT double binomialp(int x, int n, double p);
/*-------------------------------------------------------------------
 Returns the value of the Binomial random variable distribution
 function with n trials and p chance of success at the value x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void binomialv(double pp, int n, double p, int *b0, int *b1);
/*-------------------------------------------------------------------
 Returns two critical points of the Binomial random variable
 distribution function with n trials and p chance of success such
 that P(B <= b0) <= pp <= P(B <= b1).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double geomp(int x, double p);
/*-------------------------------------------------------------------
 Returns the value of the Geometric random variable distribution
 function with probability of success p at the value x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void geomv(double pp, double p, int *x0, int *x1);
/*-------------------------------------------------------------------
 Returns two critical points of the Geometric random variable
 distribution function with probability of success p such that
 P(X <= x0) <= pp <= P(X <= x1).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double hyperp(int x, int n, int N, int M);
/*-------------------------------------------------------------------
 Returns the probability from the Hypergeometric random variable
 distribution function with a population size M, containing N
 successes, with a sample of size n, and observing x or less
 successes.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void hyperv(double p, int n, int N, int M, int *x0, int *x1);
/*-------------------------------------------------------------------
 Returns two critical points of the Hypergeometric random variable
 distribution function with a population size M, containing N
 successes, with a sample of size n, such that P(X <= x0) <= p
 <= P(X <= x1).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double negbnlp(int x, int r, double p);
/*-------------------------------------------------------------------
 Returns the value of the negative binomial random variable
 distribution function with r successes each with probability of p
 at the value of x failures.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void negbnlv(double p, int r, double pp, int *x0, int *x1);
/*-------------------------------------------------------------------
 Returns two critical points of the negative binomial random
 variable distribution function with r successes each with
 probability of pp, such that P(X <= x0) <= p <= P(X <= x1).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double poissonp(int x, double l);
/*-------------------------------------------------------------------
 Returns the value of the Poisson random variable distribution
 function with mean l at the value x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void poissonv(double p, double l, int *x0, int *x1);
/*-------------------------------------------------------------------
 Returns two critical points of the Poisson random variable
 distribution function with mean l such that P(X <= x0) <= p <=
 P(X <= x1).
-------------------------------------------------------------------*/

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/
