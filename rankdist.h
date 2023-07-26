/*===================================================================
 rankdist.h

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#ifndef _RANKDIST_H_
#define _RANKDIST_H_

#include "numerics.h"

NUMERICS_EXPORT double wilcoxonp(int w, int n);
/*-------------------------------------------------------------------
 Returns the value of the null distribution function of the
 Wilcoxon Signed Rank Test at the value w with sample size n.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void wilcoxonv(double p, int n, int *w0, int *w1);
/*-------------------------------------------------------------------
 Returns two critical points of the Wilcoxon Signed Rank Test null
 distribution function with sample size of n trials such that
 P(W <= w0) <= p <= P(W <= w1).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double wilmannwhitp(int w, int m, int n);
/*-------------------------------------------------------------------
 Returns the value of the null distribution function of the
 Wilcoxon-Mann-Whitney Rank Sum Test at the value w with sample
 sizes m and n.  w is a value from the Wilcoxon form of the test.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void wilmannwhitv(double p, int m, int n, int *w0, int *w1);
/*-------------------------------------------------------------------
 Returns two critical points of the Wilcoxon-Mann-Whitney Rank Sum
 Test null distribution function with sample size of n trials such
 that P(W <= w0) <= p <= P(W <= w1).  w0 and w1 are values from the
 Wilcoxon from of the test.
-------------------------------------------------------------------*/

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/
