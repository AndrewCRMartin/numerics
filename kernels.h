/*===================================================================
 kernels.h

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#ifndef _KERNELS_H_
#define _KERNELS_H_

#include "numerics.h"

NUMERICS_EXPORT double biweight(double x);
/*-------------------------------------------------------------------
 Returns the value of the Biweight kernel at the point x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double cosine(double x);
/*-------------------------------------------------------------------
 Returns the value of the Cosine kernel at the point x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double epanech(double x);
/*-------------------------------------------------------------------
 Returns the value of the Epanechnikov kernel at the point x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double gaussian(double x);
/*-------------------------------------------------------------------
 Returns the value of the Gaussian kernel at the point x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double rectangle(double x);
/*-------------------------------------------------------------------
 Returns the value of the Rectangular kernel at the point x.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double triangle(double x);
/*-------------------------------------------------------------------
 Returns the value of the Triangular kernel at the point x.
-------------------------------------------------------------------*/

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/

