/*===================================================================
 noncdist.h

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#ifndef _NONCDIST_H_
#define _NONCDIST_H_

#include "numerics.h"

NUMERICS_EXPORT double ncstudtp(double x, double df, double delta);
/*-------------------------------------------------------------------
 Returns the value of Non-central Student's T random variable
 distribution function with df degrees of freedom and
 non-centrality parameter delta at the value x.
-------------------------------------------------------------------*/

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/
