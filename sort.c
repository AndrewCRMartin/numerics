/*===================================================================
 sort.c

 Version 1.1

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  April 10, 1999
===================================================================*/

#include "sort.h"

#define MAX_QSORT 10

#define iter_swap(a, b) \
{                       \
    double tmp;         \
    tmp = *(a);         \
    *(a) = *(b);        \
    *(b) = tmp;         \
}

NUMERICS_EXPORT void isort1(double* first, double* last)
{
    double* iter;
    double* curr;
    double value;
    
    for(iter = first + 1; iter != last; ++iter){
		value = *iter;
		curr = iter - 1;
        while(curr >= first && *curr > value){
			*(curr + 1) = *curr;
            curr--;
        }
        *(curr + 1) = value;
    }
}

NUMERICS_EXPORT void isort2(double* first1, double* last1, double* first2)
{
    double* iter1;
    double* iter2;
    double* j1;
    double* j2;
    double m1;
    double m2;
    
    for(iter1 = first1 + 1, iter2 = first2 + 1; iter1 != last1; ++iter1, ++iter2){
        m1 = *iter1;
        m2 = *iter2;
        j1 = iter1 - 1;
        j2 = iter2 - 1;
        while(j1 >= first1 && *j1 > m1){
            *(j1 + 1) = *j1;
            *(j2 + 1) = *j2;
            --j1;
            --j2;
        }
        *(j1 + 1) = m1;
        *(j2 + 1) = m2;
    }
}

NUMERICS_EXPORT void qsort1(double* first, double* last)
{
    double* loSwap;
    double* hiSwap;
    
    if(last - first > MAX_QSORT){
        loSwap = first + 1;
        hiSwap = last - 1;
        
        do {
            while(loSwap <= hiSwap && *loSwap <= *first){
                ++loSwap;
            }
            while(*hiSwap > *first){
                --hiSwap;
            }
            if(loSwap < hiSwap){
                iter_swap(loSwap, hiSwap);
            }
        } while(loSwap < hiSwap);
        
        iter_swap(first, hiSwap);
        
        if(first < hiSwap - 1){
            qsort1(first, hiSwap);
        }
        if(hiSwap + 2 < last){
			qsort1(hiSwap + 1, last);
        }
    } else {
		isort1(first, last);
	}
}

NUMERICS_EXPORT void qsort2(double* first1, double* last1, double* first2)
{
    double* loSwap, *hiSwap;
    double* lo2, *hi2;
    
    if(last1 - first1 > MAX_QSORT){
        loSwap = first1 + 1;
        hiSwap = last1 - 1;
        lo2 = first2 + 1;
        hi2 = first2 + (last1 - first1) - 1;
        
        do {
            while(loSwap <= hiSwap && *loSwap <= *first1){
                ++loSwap;
                ++lo2;
            }
            while(*hiSwap > *first1){
                --hiSwap;
                --hi2;
            }
            if(loSwap < hiSwap){
                iter_swap(loSwap, hiSwap);
                iter_swap(lo2, hi2);
            }
        } while(loSwap < hiSwap);
        
        iter_swap(first1, hiSwap);
        iter_swap(first2, hi2);
        
        if(first1 < hiSwap - 1){
            qsort2(first1, hiSwap, first2);
        }
        if(hiSwap + 2 < last1){
			qsort2(hiSwap + 1, last1, hi2 + 1);
        }
    } else {
		isort2(first1, last1, first2);
	}
}

NUMERICS_EXPORT void qsort2key(double* first1, double* last1, double* first2)
{
	double* key;
	double* lo, *hi;
     
	qsort2(first1, last1, first2);

	hi = first2;
	while(first1 < last1){
         key = first1;
         lo = hi;
         while(first1 < last1 && *first1 == *key){
             ++first1;
             ++hi;
         }
         if(hi - lo > 0){
             qsort1(lo, hi);
         }
	}
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
 Version 1.1 - 04/10/1999 - Added isort1 and isort2.
                            Applied isort1 and isort2 to small groups
                            in qsort1 and qsort2.
===================================================================*/
