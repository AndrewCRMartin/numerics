/*===================================================================
 numerics.h

 Version 1.1

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  December 22, 1998
===================================================================*/

#ifndef _NUMERICS_H_
#define _NUMERICS_H_

#ifdef NUMERICS_DLL
	#define NUMERICS_EXPORT __declspec (dllexport)
#else
	#define NUMERICS_EXPORT
#endif

#include <float.h>

#define NUMERICS_PI        3.14159265358979323846
#define NUMERICS_E         2.71828182845904523536
#define NUMERICS_EULER     0.5772156649
#define NUMERICS_ITMAX     100
#define NUMERICS_MAX_ERROR 5.0e-9
#define NUMERICS_FLOAT_MIN DBL_MIN
#define NUMERICS_FLOAT_MAX DBL_MAX

#ifndef BOOL
#define BOOL  int
#define FALSE 0
#define TRUE  1
#endif

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
 Version 1.1 - 12/22/1998 - Fixed bug in NUMERICS_EXPORT definition.
===================================================================*/
