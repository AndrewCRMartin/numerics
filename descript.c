/*===================================================================
 descript.c

 Version 1.1

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  May 09, 1999
===================================================================*/

#include <math.h>

#include "datamanp.h"
#include "descript.h"
#include "residual.h"
#include "sort.h"
#include "vecutils.h"

NUMERICS_EXPORT double centmom(double* first, double* last, double n)
{
    int len = length(first, last);
    double* r = dvector(0, len - 1);
    double ret;
    
    copy(first, last, r);
    nresid(r, r + len, mean(first, last), n);
    ret = mean(r, r + len);
    
    free_dvector(r, 0);
    
    return ret;
}

NUMERICS_EXPORT double coeffvar(double* first, double* last)
{
    return stddev(first, last) / mean(first, last) * 100.0;
}

NUMERICS_EXPORT double geomean(double* first, double* last)
{
    double p = 1.0 / length(first, last);
    double m = 1.0;
     
    while(first != last){
        m *= pow(*first, p);
        ++first;
    };

    return m;
}

NUMERICS_EXPORT double harmean(double* first, double* last)
{
    double n = length(first, last);
    double m = 0.0;
    
    while(first != last){
        if(*first){
            m += 1.0 / *first;
        }
        ++first;
    }
    
    return m;
}

NUMERICS_EXPORT double iqr(double* first, double* last)
{
    return quartile3(first, last) - quartile1(first, last);
}

NUMERICS_EXPORT double kurtosis(double* first, double* last)
{
    return centmom(first, last, 4) / pow(variance(first, last), 2);
}

NUMERICS_EXPORT int length(double* first, double* last)
{
    return last - first;
}

NUMERICS_EXPORT double maximum(double* first, double *last)
{
    double tmp = *first;
    
    while(first < last){
	   if(tmp < *first){
            tmp = *first;
	   }
	   ++first;
    }
    
    return tmp;
}

NUMERICS_EXPORT double mean(double* first, double* last)
{
	double m = 0.0;
	double n = 0.0;
	while(first != last){
		n += 1.0;
		m += (*first - m) / n;
		++first;
	}
    return m;
}

NUMERICS_EXPORT double meandev(double* first, double* last)
{
    int len = length(first, last);
    double* r = dvector(0, len - 1);
    double ret;
    
    copy(first, last, r);
    aresid(r, r + len, mean(first, last));
    ret = mean(r, r + len);
    
    free_dvector(r, 0);
    
    return ret;
}

NUMERICS_EXPORT double meddev(double* first, double* last)
{
    int len = length(first, last);
    double* r = dvector(0, len - 1);
    double ret;
    
    copy(first, last, r);
    aresid(r, r + len, median(first, last));
    ret = median(r, r + len);
    
    free_dvector(r, 0);
    
	return ret;
}

NUMERICS_EXPORT double median(double* first, double* last)
{
    double *tmp, a, b;
    int n = length(first, last);
    
    tmp = dvector(0, n-1);
    
    copy(first, last, tmp);
    
    qsort1(tmp, tmp + n);
    a = tmp[n/2-1];
    b = tmp[n/2];
    
    free_dvector(tmp, 0);
    
    if(n % 2 == 0){
	   return a + 0.5 * (b - a);
    } else {
        return b;
    }
}

NUMERICS_EXPORT double minimum(double* first, double* last)
{
    double tmp = *first;
    
    while(first < last){
	   if(tmp > *first){
            tmp = *first;
	   }
	   ++first;
    }

    return tmp;
}

NUMERICS_EXPORT double quantile(double* first, double* last, double q)
{
    int i, n = length(first, last);
    double* vec = dvector(0, n - 1);
    double u, f, v, ans;
    
    if(q < 0.0){
        q = 0.0;
    }
    if(q > 1.0){
        q = 1.0;
    }
    
    copy(first, last, vec);
    qsort1(vec, vec + n);
    if(q < (v = 1.0 / (2.0 * n))){
        ans = *vec;
    } else if(q > 1.0 - v){
        ans = *(vec + n - 1);
    } else {
        u = n * q + .5;
        i = (int)u;
        f = u - i;
        ans = (1.0 - f) * *(vec + i - 1) + f * *(vec + i);
    }
    
    free_dvector(vec, 0);
    
    return ans;
}

NUMERICS_EXPORT double quartile1(double* first, double* last)
{
    return quantile(first, last, .25);
}

NUMERICS_EXPORT double quartile3(double* first, double* last)
{
    return quantile(first, last, .75);
}

NUMERICS_EXPORT double range(double* first, double* last)
{
    return maximum(first, last) - minimum(first, last);
}

NUMERICS_EXPORT double rms(double* first, double* last)
{
    return sqrt(sum2(first, last) / length(first, last));
}

NUMERICS_EXPORT double skewness(double* first, double* last)
{
    return centmom(first, last, 3) / pow(stddev(first, last), 3);
}

NUMERICS_EXPORT double stddev(double* first, double* last)
{
    return sqrt(variance(first, last));
}

NUMERICS_EXPORT double stderrmean(double* first, double* last)
{
    return sqrt(variance(first, last) / (double)length(first, last));
}

NUMERICS_EXPORT double sum(double* first, double* last)
{
    double s = 0.0;

    while(first < last){
        s += *first;
        ++first;
    }

	return s;
}

NUMERICS_EXPORT double sum2(double* first, double* last)
{
    double s = 0.0;

    while(first < last){
        s += *first * *first;
        ++first;
    }

    return s;
}

NUMERICS_EXPORT double sump(double* first1, double* last1, double* first2)
{
    double s = 0.0;

    while(first1 < last1){
        s += *first1 * *first2;
        ++first1;
        ++first2;
    }
    
    return s;
}

NUMERICS_EXPORT double sumproduct(double* first1, double* last1, double* first2)
{
    int n = length(first1, last1);

    return sump(first1, last1, first2) - sum(first1, last1) *
        sum(first2, first2 + n) / (double)n;
}

NUMERICS_EXPORT double sumsquare(double* first, double* last)
{
    int n = length(first, last);
    double* r = dvector(0, n - 1);
    double ep, s2;
    
    copy(first, last, r);
    resid(r, r + n, mean(first, last));
    
    ep = sum(r, r + n);
    s2 = sum2(r, r + n);
    
    free_dvector(r, 0);
    
    return s2 - ep * ep / (double)n;
};

NUMERICS_EXPORT double trimmean(double* first, double* last, double f)
{
    return trimmean2(first, last, f, f);
}

NUMERICS_EXPORT double trimmean2(double* first, double* last, double f1, double f2)
{
    int n = length(first, last);
    double* v = dvector(0, n - 1);
    int n1 = (int)floor(f1 * n), n2 = (int)floor(f2 * n);
    double ret;
    
    copy(first, last, v);
    qsort1(v, v + n);
    
    ret = mean(v + n1, v + n - n2);
    
    free_dvector(v, 0);
    
    return ret;
}

NUMERICS_EXPORT double variance(double* first, double* last)
{
    return sumsquare(first, last) / (double)(length(first, last) - 1);
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
             - 05/09/1999 - Change implementation of mean
===================================================================*/
