/*===================================================================
 descript.h

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:
===================================================================*/

#ifndef _DESCRIPT_H_
#define _DESCRIPT_H_

#include "numerics.h"

NUMERICS_EXPORT double centmom(double* first, double* last, double n);
/*-------------------------------------------------------------------
 Returns the n-th central moment of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double coeffvar(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the coefficient of variation of the elements in
 [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double geomean(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the geometric mean of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double harmean(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the harmonic mean of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double iqr(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the interquartile range of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double kurtosis(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the kurtosis of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT int length(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the kurtosis of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double maximum(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the maximum of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double mean(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the mean of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double meandev(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the mean deviation of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double meddev(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the median deviation of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double median(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the median of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double minimum(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the minimum of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double quantile(double* first, double* last, double q);
/*-------------------------------------------------------------------
 Returns the q*100 percentile of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double quartile1(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the first quartile of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double quartile3(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the third quartile of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double range(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the range of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double rms(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the root mean square of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double skewness(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the skewness of the elements in
 [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double stddev(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the sample standard deviation of the elements in
 [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double stderrmean(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the standard error of the mean of the elements in
 [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double sum(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the sum of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double sum2(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the sum of the squares of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double sump(double* first1, double* last1, double* first2);
/*-------------------------------------------------------------------
 Returns the sum of the products of the corresponding elements in
 [first1, last1) and [first2, first2 + (last1 - first1) ).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double sumproduct(double* first1, double* last1, double* first2);
/*-------------------------------------------------------------------
 Returns the sums of products of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double sumsquare(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the sums of squares of the elements in [first, last).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double trimmean(double* first, double* last, double f);
/*-------------------------------------------------------------------
 Returns the mean of the elements in [first, last) with f*100
 percent of the lowest and highest terms removed.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double trimmean2(double* first, double* last, double f1, double f2);
/*-------------------------------------------------------------------
 Returns the mean of the elements in [first, last) with f1*100
 percent of the lowest terms and f2*100 percent of the highest
 terms removed.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double variance(double* first, double* last);
/*-------------------------------------------------------------------
 Returns the sample variance of the elements in [first, last).
-------------------------------------------------------------------*/

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/
