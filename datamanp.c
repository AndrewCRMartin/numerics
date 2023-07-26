/*===================================================================
 datamanp.c

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#include "datamanp.h"

NUMERICS_EXPORT void copy(double* first, double* last, double* dest)
{
    while(first < last){
        *dest = *first;
        ++first;
        ++dest;
    }
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/
