/*===================================================================
 sort.h

 Version 1.1

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  April 10, 1999
===================================================================*/

#ifndef _SORT_H_
#define _SORT_H_

#include "numerics.h"

NUMERICS_EXPORT void isort1(double* first, double* last);
/*-------------------------------------------------------------------
 Sorts the elements of [first, last) into ascending numerical order
 using the insertion sort.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void isort2(double* first1, double* last1, double* first2);
/*-------------------------------------------------------------------
 Sorts the elements of [first, last) into ascending numerical order
 using the insertion sort while making the same rearrangement of the
 elements of [first2, first2 + (last1 - first1) ).
-------------------------------------------------------------------*/

NUMERICS_EXPORT void qsort1(double* first, double* last);
/*-------------------------------------------------------------------
 Sorts the elements of [first, last) into ascending numerical order
 using the Quicksort.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void qsort2(double* first1, double* last1, double* first2);
/*-------------------------------------------------------------------
 Sorts the elements of [first, last) into ascending numerical order
 using the Quicksort while making the corresponding rearrangement
 of the elements of [first2, first2 + (last1 - first1) ).
-------------------------------------------------------------------*/

NUMERICS_EXPORT void qsort2key(double* first1, double* last1, double* first2);
/*-------------------------------------------------------------------
 Sorts the element of [first1, last1) and [first2, first2 + (last1
 - first1) ) into ascending order using [first1, last1) as the
 primary key and [first2, first2 + (last1 - first1) ) as the
 secondary key.
-------------------------------------------------------------------*/

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
 Version 1.1 - 04/10/1999 - Added isort1 and isort2.
===================================================================*/
