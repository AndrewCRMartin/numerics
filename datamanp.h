/*===================================================================
 datamanp.h

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#ifndef _DATAMANP_H_
#define _DATAMANP_H_

#include "numerics.h"

NUMERICS_EXPORT void copy(double* first, double* last, double* dest);
/*-------------------------------------------------------------------
 Copies the elements in [first1, last1) to
 [dest, dest + (last - first)).
-------------------------------------------------------------------*/

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/
