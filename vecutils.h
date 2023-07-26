/*===================================================================
 vecutils.h

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#ifndef _VECUTILS_H_
#define _VECUTILS_H_

#include "numerics.h"

NUMERICS_EXPORT double *dvector(int nl, int nh);
/*-------------------------------------------------------------------
 Allocate an double vector with subscript range v[nl..nh].
-------------------------------------------------------------------*/

NUMERICS_EXPORT void free_dvector(double *v, int nl);
/*-------------------------------------------------------------------
 Deallocate an double vector associated with dvector().
-------------------------------------------------------------------*/

NUMERICS_EXPORT int *ivector(int nl, int nh);
/*-------------------------------------------------------------------
 Allocate an int vector with subscript range v[nl..nh].
-------------------------------------------------------------------*/

NUMERICS_EXPORT void free_ivector(int *v, int nl);
/*-------------------------------------------------------------------
 Deallocate an int vector associated with ivector().
-------------------------------------------------------------------*/

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/98 - New.
===================================================================*/

