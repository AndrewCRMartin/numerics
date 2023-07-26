/*===================================================================
 deviate.h

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#ifndef _DEVIATE_H_
#define _DEVIATE_H_

#include "numerics.h"

NUMERICS_EXPORT double ran0(long *idum);
/*--------------------------------------------------------------------
 "Minimal" random number generator of Park and Miller.  Returns a
 uniform random deviate between 0.0 and 1.0.  Set or reset idum to
 any integer value to initialize the sequence; idum must not be
 altered between calls for successive deviates in a sequence.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double ran1(long *idum);
/*--------------------------------------------------------------------
 "Minimal" random number generator of Park and Miller with
 Bays-Durham shuffle and added safeguards.  Returns a uniform random
 deviate between 0.0 and 1.0 (exclusive of the endpoint values).
 Call with idum a negative integer to initialize; thereafter, do not
 alter idum between successive deviates in a sequence.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double ran2(long *idum);
/*--------------------------------------------------------------------
 Long period (> 2 X 10^18) random number generator of L'Ecuyer with
 Bays-Durham shuffle and added safeguards.  Returns a uniform random
 deviate between 0.0 and 1.0 (exclusive of the endpoint values).
 Call with idum a negative integer to initialize; thereafter, do not
 alter idum between successive deviates in a sequence.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double ran3(long *idum);
/*--------------------------------------------------------------------
 Returns a uniform random deviate between 0.0 and 1.0.  Set idum to
 any negative value to initialize or reinitialize the sequence.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double berdev(double p, long *idum);
/*--------------------------------------------------------------------
 Returns as a doubleing-point number an integer value that is a
 random deviate drawn from a bernoulli distribution with probability
 p, using RAND(idum) as a source of uniform deviates.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double betadev(double a, double b, long *idum);
/*--------------------------------------------------------------------
 Returns a beta distributed deviate with parameters a and b, using
 RAND(idum) as the source of uniform deviates.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double bnldev(int n, double p, long *idum);
/*--------------------------------------------------------------------
 Returns as a doubleing-point number an integer value that is a
 random deviate drawn from a binomial distribution of n trials each
 of probability p, using RAND(idum) as a source of uniform deviates.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double caudev(double m, double s, long *idum);
/*--------------------------------------------------------------------
 Returns a cauchy distributed deviate with median m and scale s,
 using RAND(idum) as the source of uniform deviates.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double chidev(double v, long *idum);
/*--------------------------------------------------------------------
 Returns a chi-squared distributed, random deviatate with v degrees
 of freedom, using gamdev(v/2, 2, idum) as the source of gamma
 deviates.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double dexpdev(double m, double s, long *idum);
/*--------------------------------------------------------------------
 Returns a double exponential distributed, random deviate of mean m
 and scale parameter s, using RAND(idum) as the source of uniform
 deviates.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double expdev(double b, long *idum);
/*--------------------------------------------------------------------
 Returns an exponentially distributed, random deviate of mean b,
 using RAND(idum) as the source of uniform deviates.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double gamdev(double a, double b, long *idum);
/*--------------------------------------------------------------------
 Returns a gamma distributed, random deviatate with parameters a and
 b, using RAND(idum) as the source of uniform deviates.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double fdev(double v1, double v2, long *idum);
/*--------------------------------------------------------------------
 Returns a F distributed, random deviatate with v1 degrees of
 freedom numerator and v2 degrees of freedom denominator using
 RAND(idum) as the source of uniform deviates.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double gasdev(double m, double v, long *idum);
/*--------------------------------------------------------------------
 Returns a normal distributed deviate with mean m and variance v,
 using RAND(idum) as the source of uniform deviates.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double logdev(double m, double b, long *idum);
/*-------------------------------------------------------------------
 Returns a logistic dostrobited deviate with mean m and beta b,
 using RAND(idum) as the source of uniform deviates.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double nchidev(double df, double xnonc, long *idum);
/*-------------------------------------------------------------------
 Returns a non-central chi-squared distributed deviate with df
 degrees of freedom and non-centrallity parameter xnonc using
 RAND(idum) as the source of uniform deviates.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double nfdev(double dfn, double dfd, double xnonc, long *idum);
/*-------------------------------------------------------------------
 Returns a non-central F distributed deviate with dfn degrees of
 freedom numerator, dfd degrees of freedom denominator, and
 non-centrallity parameter xnonc using RAND(idum) as the source of
 uniform deviates.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void orddev(double *x, int n, long *idum);
/*--------------------------------------------------------------------
 Returns as the array x[0..n-1], n order statistic random deviates
 from a uniform distibution (0,1).
--------------------------------------------------------------------*/

NUMERICS_EXPORT double poisdev(double m, long *idum);
/*--------------------------------------------------------------------
 Returns an poisson distributed deviate with mean m, using
 RAND(idum) as the source of uniform deviates.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double tdev(int v, long *idum);
/*--------------------------------------------------------------------
 Returns a Student's t distributed deviate with v degrees of
 freedom, using RAND(idum) as the source of uniform deviates.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double unifdev(double lo, double hi, long *idum);
/*--------------------------------------------------------------------
 Returns an Uniform distributed deviate in the range of (lo, hi),
 using RAND(idum) as the source of uniform(0,1) deviates.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double weibdev(double a, double b, long *idum);
/*--------------------------------------------------------------------
 Returns an weibull distributed deviate with parameters a and b,
 using RAND(idum) as the source of uniform deviates.
--------------------------------------------------------------------*/

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/
