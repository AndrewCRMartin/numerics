/*===================================================================
 mathx.h

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#ifndef _MATHX_H_
#define _MATHX_H_

#include "numerics.h"

NUMERICS_EXPORT double acosh(double x);
/*--------------------------------------------------------------------
 Returns the inverse hyperbolic cosine of x.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double asinh(double x);
/*--------------------------------------------------------------------
 Returns the inverse hyperbolic sine of x.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double atanh(double x);
/*--------------------------------------------------------------------
 Returns the inverse hyperbolic tangent of x.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double beta(double z, double w);
/*--------------------------------------------------------------------
 Returns the value of the beta function B(z,w).
--------------------------------------------------------------------*/

NUMERICS_EXPORT double betai(double a, double b, double x);
/*--------------------------------------------------------------------
 Returns the incomplete beta function Ix(a,b).
--------------------------------------------------------------------*/

NUMERICS_EXPORT double bico(int n, int x);
/*--------------------------------------------------------------------
 Returns the binomial coefficient nCx as a floating-point number.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double bicoln(int n, int k);
/*--------------------------------------------------------------------
 Returns the Ln(nCx).
--------------------------------------------------------------------*/

NUMERICS_EXPORT double erff(double x);
/*--------------------------------------------------------------------
 Returns the error function erf(x).
--------------------------------------------------------------------*/

NUMERICS_EXPORT double erffc(double x);
/*--------------------------------------------------------------------
 Returns the complementary error function erfc(x).
--------------------------------------------------------------------*/

NUMERICS_EXPORT double erfcc(double x);
/*--------------------------------------------------------------------
 Returns the complementary error function erfc(x) with fractional
 error everywhere less than 1.2 X 10^-7.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double expint(int n, double x);
/*--------------------------------------------------------------------
 Evaluates the exponential integral En(x).
--------------------------------------------------------------------*/

NUMERICS_EXPORT double factrl(int n);
/*--------------------------------------------------------------------
 Returns the value n! as a floating-point number.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double gammln(double xx);
/*--------------------------------------------------------------------
 Returns the value Ln(Gamma(xx)).
--------------------------------------------------------------------*/

NUMERICS_EXPORT double gammp(double a, double x);
/*--------------------------------------------------------------------
 Returns the incomplete gamma function P(a,x).
--------------------------------------------------------------------*/

NUMERICS_EXPORT double gammq(double a, double x);
/*--------------------------------------------------------------------
 Returns the incomplete gamma function Q(a,x)=1-P(a,x).
--------------------------------------------------------------------*/

NUMERICS_EXPORT double factln(int n);
/*--------------------------------------------------------------------
 Returns Ln(n!).
--------------------------------------------------------------------*/

NUMERICS_EXPORT double log1x(double x);
/*--------------------------------------------------------------------
 Returns Ln(1 + x).
--------------------------------------------------------------------*/

NUMERICS_EXPORT double pythag(double a, double b);
/*--------------------------------------------------------------------
 Returns sqrt(a^2 + b^2) without destructive overflow or underflow.
--------------------------------------------------------------------*/

NUMERICS_EXPORT double sign(double x);
/*--------------------------------------------------------------------
 Returns:
 -1 if x < 0
  0 if x = 0
  1 if x > 0
--------------------------------------------------------------------*/

#endif

/*====================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
====================================================================*/

