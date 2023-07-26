/*===================================================================
 ranking.h

 Version 1.1

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  May 18, 1999
===================================================================*/

#ifndef _RANKING_H_
#define _RANKING_H_

#include "numerics.h"

NUMERICS_EXPORT void index(double *first, double* last, int *indx);
/*-------------------------------------------------------------------
 Indexes the elements in [first, last) and places the index values
 in [indx, indx + (last - first)).  The input elements of [first,
 last) are not changed.
-------------------------------------------------------------------*/

enum {RANK_MEAN, RANK_LOW, RANK_HIGH};

NUMERICS_EXPORT void rank(double* first, double* last, double* irank, int ties);
/*-------------------------------------------------------------------
 Ranks the elements in [first, last) by returning the table of
 ranks in [irank, irank + (last - first)).  The input elements of
 [first, last) are not changed.  The parameter ties determines the
 method used to assign ranks to tied values.  RANK_MEAN uses the
 mean of the corresponding ranks.  RANK_HIGH uses the largest of
 the corresponding ranks.  RANK_LOW uses the smallest of the
 corresponding ranks.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void arank(double* first, double* last, double* irank, int ties);
/*-------------------------------------------------------------------
 Ranks the absolute values of the elements in [first, last) by
 returning the table of ranks in [irank, irank + (last - first)).
 The input elements of [first, last) are not changed.  The
 parameter ties determines the method used to assign ranks to tied
 values as described above.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void abrank(double* first, double* last, double* irank, int ties);
/*-------------------------------------------------------------------
 Ranks the values of the elements in [first, last) using the scoring
 method used in the Ansari-Bradley test returning the table of ranks
 in [irank, irank + (last - first)).  The input elements of
 [first, last) are not changed.  The parameter ties determines the
 method used to assign ranks to tied values as described above.
-------------------------------------------------------------------*/

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
 Version 1.1 - 05/18/1999 - Added abrank.
===================================================================*/
