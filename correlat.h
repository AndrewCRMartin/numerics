/*===================================================================
 correlat.h

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare
 
 Created:  August 28, 1998
 Revised:  
===================================================================*/

#ifndef _CORRELAT_H_
#define _CORRELAT_H_

#include "numerics.h"

NUMERICS_EXPORT double covariance(double* first1, double* last1, double* first2);
/*-------------------------------------------------------------------
 Returns the covariance of the corresponding elements in
 [first1, last1) and [first2, first2 + (last1 - first1) ).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double correlation(double* first1, double* last1, double* first2);
/*-------------------------------------------------------------------
 Returns the Pearson's product moment correlation coefficient of
 the corresponding elements in [first1, last1) and [first2,
 first2 + (last1 - first1) ).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double kendall(double* first1, double* last1, double* first2);
/*-------------------------------------------------------------------
 Returns Kenall's Tau correlation coefficient of the corresponding
 values in [first1, last1) and [first2, first2 + (last1 - first1)).
-------------------------------------------------------------------*/

NUMERICS_EXPORT double spearman(double* first1, double* last1, double* first2);
/*-------------------------------------------------------------------
 Returns Spearman's correlation coefficient of the corresponding
 values in [first1, last1) and [first2, first2 + (last1 - first1)).
 Note:  Ties in the values are not accounted for.
-------------------------------------------------------------------*/

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/

