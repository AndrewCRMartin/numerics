/*===================================================================
 matutils.h

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#ifndef _MATUTILS_H_
#define _MATUTILS_H_

#include "numerics.h"

NUMERICS_EXPORT double **dmatrix(int nrl, int nrh, int ncl, int nch);
/*-------------------------------------------------------------------
 Allocate an double matrix with subscript range
 m[nrl..nrh][ncl..nch].
-------------------------------------------------------------------*/

NUMERICS_EXPORT void free_dmatrix(double **m, int nrl, int nrh, int ncl);
/*-------------------------------------------------------------------
 Deallocate an double matrix associated with dmatrix().
-------------------------------------------------------------------*/

NUMERICS_EXPORT int **imatrix(int nrl, int nrh, int ncl, int nch);
/*-------------------------------------------------------------------
 Allocate an int matrix with subscript range m[nrl..nrh][ncl..nch].
-------------------------------------------------------------------*/

NUMERICS_EXPORT void free_imatrix(int **m, int nrl, int nrh, int ncl);
/*-------------------------------------------------------------------
 Deallocate an int matrix associated with imatrix().
-------------------------------------------------------------------*/

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/98 - New.
===================================================================*/

