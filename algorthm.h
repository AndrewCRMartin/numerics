/*===================================================================
 algorthm.h

 Version 1.1

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  May 18, 1999
===================================================================*/

#ifndef _ALGORTHM_H_
#define _ALGORTHM_H_

#include "numerics.h"

NUMERICS_EXPORT double average(double x, double y);
/*-------------------------------------------------------------------
 Returns the average of elements in x and y.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void pairwdiff(double* first1, double* last1, double* first2, double* last2, double* dest);
/*-------------------------------------------------------------------
 Computes the pairwise differences between the elements in
 [first1, last1) and the elements in [first2, last2) and places them
 in dest.  dest must be large enough to hold all of the m * n
 differences, where m = last1 - first1 and n = last2 - first2.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void walshavg(double* first, double* last, double* dest);
/*-------------------------------------------------------------------
 Computes the Walsh averages of the elements in [first, last) and
 places them in dest.  dest must be large enough to hold all of the
 n averages ( n = k * (k + 1) / 2, where k = last - first).
-------------------------------------------------------------------*/

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
 Version 1.1 - 05/18/1999 - Added pairwdiff.
===================================================================*/


