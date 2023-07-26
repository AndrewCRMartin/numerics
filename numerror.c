/*===================================================================
 algorthm.hpp

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#include <stdio.h>

#include "numerror.h"

NUMERICS_EXPORT void NUMERICS_ERROR(const char *func, const char *msg)
{
    fprintf(stderr, "Numerics Error in routine %s\n", func);
    fprintf(stderr, "    %s\n", msg);
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/
