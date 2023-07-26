/*===================================================================
 residual.h

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#ifndef _RESIDUAL_H_
#define _RESIDUAL_H_

#include "numerics.h"

NUMERICS_EXPORT void resid(double* first, double* last, double val);
/*-------------------------------------------------------------------
 Replaces the elements in [first, last) with thier deviations from
 val.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void aresid(double* first, double* last, double val);
/*-------------------------------------------------------------------
 Replaces the elements in [first, last) with thier absolute
 deviations from val.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void nresid(double* first, double* last, double val, double n);
/*-------------------------------------------------------------------
 Replaces the elements in [first, last) with thier deviations from
 val raised to the n-th power.
-------------------------------------------------------------------*/

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/
