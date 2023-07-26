/*===================================================================
 discdist.c

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

#include "discdist.h"
#include "domain.h"
#include "mathx.h"
#include "normdist.h"
#include "numerror.h"

NUMERICS_EXPORT double binomialp(int x, int n, double p)
{
	if(!isZeroOne(p)){
		NUMERICS_ERROR("binomialp", "Invalid p parameter");
	} else if(n < 1){
		NUMERICS_ERROR("binomialp", "Invalid n parameter");
	} else if(isNegative(x)){
		return 0.0;
	} else if(x >= n){
		return 1.0;
	}
	return 1.0 - betai(x + 1, n - x, p);
}

NUMERICS_EXPORT void binomialv(double pp, int n, double p, int *b0, int *b1)
{
    int b;
    
	if(!isZeroOne(p)){
		NUMERICS_ERROR("binomialv", "Invalid p parameter");
	} else if(!isZeroOne(pp)){
		NUMERICS_ERROR("binomialv", "Invalid pp parameter");
	} else if(n < 1){
		NUMERICS_ERROR("binomialv", "Invalid n parameter");
	} else {
		b = (int)(n * p);
		while(binomialp(b, n, p) > pp) b--;
		while(binomialp(b, n, p) <= pp) b++;
		if(b == 0){
	        *b0 = b;
		} else *b0 = b-1;
		*b1 = b;
	}
}

NUMERICS_EXPORT double geomp(int x, double p)
{
	if(!isZeroOne(p)){
		NUMERICS_ERROR("geomp", "Invalid p parameter");
	} else if(isNonPositive(x)){
		return 0.0;
	}
	return (1.0 - pow(1.0 - p, x));
}

NUMERICS_EXPORT void geomv(double pp, double p, int *x0, int *x1)
{
    double x;
    
	if(!isZeroOne(p)){
		NUMERICS_ERROR("geomv", "Invalid p parameter");
	} else if(!isZeroOne(pp)){
		NUMERICS_ERROR("geomv", "Invalid pp parameter");
	} else {
		x = log(1.0 - pp) / log(1.0 - p);
		*x0 = (int)floor(x);
		if(*x0 == 0) *x0 = 1;
		*x1 = (int)ceil(x);
	}
}

NUMERICS_EXPORT double hyperp(int x, int n, int N, int M)
{
    int i1, i2, i3, i4, i5, mnk1, i, j, k, l, m, nn, k1, n1;
    double ret, r1, r2, mn, p, pt, sig;
    BOOL dir;
    
    k = N + 1;
    l = x + 1;
    m = M + 1;
    nn = n + 1;
    dir = TRUE;
    
    ret = 0.0;
    if(nn < 1 || m < nn || k < 1 || k > m) {
        return ret;
    }
    
    if(l < 1 || k - l > m - nn) {
        return ret;
    }

    ret = 1.0;
    if (l > nn || l > k) {
        return ret;
    }
    if (k == 1 || k == m || nn == 1 || nn == m) {
        return ret;
    }
    if(x == ((N < n) ? N : n)) {
        return ret;
    }
    
    p = (double)n / (double)(M - n);
    i1 = N; i2 = M - N;
    r1 = p; r2 = 1.0 / (double)p;
    if(((i1 < i2) ? i1 : i2) > ((r1 > r2) ? r1 : r2) * 16.0 && M > 1000) {
        mn = (double)(N * n) / (double)M;
        sig = sqrt(mn * ((double)(M - n) / (double)M) * ((double)(M - N) / (double)(M - 1.0)));
        r1 = (double)(x + .5 - mn) / sig;
        ret = normalp(r1, 0.0, 1.0);
    } else {
        i1 = k - 1, i2 = m - k;
        i3 = nn - 1, i4 = m - nn;
        if(((i1 < i2) ? i1 : i2) > ((i3 < i4) ? i4 : i3)){
            i = k;
            k = nn;
            nn = i;
        }
        if(m - k < k - 1) {
            dir = ! dir;
            l = nn - l + 1;
            k = m - k + 1;
        }
        if(M > 600){
            i1 = M - N;
            i2 = M - n;
            i3 = n - x;
            i4 = N - x;
            i5 = M - n - N + x;
            p = factln(n) - factln(M) + factln(i1) + factln(N) +
                factln(i2) - factln(x) - factln(i3) - factln(i4) - factln(i5);
            ret = 0.0;
            if (p >= -88.0) {
                ret = exp(p);
            }
        } else {
            i1 = l - 1;
            for(i = 1; i <= i1; ++i) {
                ret *= (double)((k - i) * (nn - i)) / (double)((l - i) * (m - i));
            }
            if (l != k) {
                j = m - nn + l;
                i1 = k - 1;
                for (i = l; i <= i1; ++i) {
                    ret *= (double) (j - i) / (double) (m - i);
                }
            }
        }
        if (ret == 0.0) {
            if (M <= 600) {
                i1 = M - n;
                i2 = n - x;
                i3 = N - x;
                i4 = M - n - N + x;
                i5 = M - N;
                p = factln(n) - factln(M) + factln(N) + factln(i1)
                    - factln(x) - factln(i2) - factln(i3) -
                    factln(i4) + factln(i5);
            }
            p += log(1.0e35);
            if (p < -88.0) {
                if (x > (n * N + n + N + 1) / (M + 2)){
                    ret = 1.0;
                }
                return ret;
            } else {
                p = exp(p);
            }
        } else {
            p = ret * 1.0e35;
        }
        pt = 0.0;
        n1 = nn - l;
        k1 = k - l;
        mnk1 = m - nn - k1 + 1;
        if (l <= k1) {
            i1 = l - 1;
            for (i = 1; i <= i1; ++i) {
                p *= (double)((l - i) * (mnk1 - i)) / (double)((n1 + i) * (k1 + i));
                pt += p;
            }
        } else {
            dir = ! dir;
            i1 = k1 - 1;
            for (j = 0; j <= i1; ++j) {
                p *= (double)((n1 - j) * (k1 - j)) / (double)((l + j) * (mnk1 + j));
                pt += p;
            }
        }
        if(dir){
            ret += pt / 1.0e35;
        } else {
            ret = 1.0 - pt / 1.0e35;
        }
    }
    
    return ret;
}

NUMERICS_EXPORT void hyperv(double p, int n, int N, int M, int *x0, int *x1)
{
    int x;
    
	if(!isZeroOne(p)){
		NUMERICS_ERROR("hyperv", "Invalid p parameter");
	} else {
		x = (int)((double)(n * N) / (double)M);
		while(hyperp(x, n, N, M) > p) x--;
		while(hyperp(x, n, N, M) <= p) x++;
		if(x <= 0){
			*x0 = x;
		} else {
			*x0 = x-1;
		}
		*x1 = x;
	}
}

NUMERICS_EXPORT double negbnlp(int x, int r, double p)
{
	if(!isZeroOne(p)){
		NUMERICS_ERROR("negbnlp", "Invalid p parameter");
	} else if(isNonPositive(r)){
		NUMERICS_ERROR("negbnlp", "Invalid r parameter");
	} else if(isNegative(x)){
		return 0.0;
	}
	return 1.0 - betai(x + 1.0, r - 1.0, p);
}

NUMERICS_EXPORT void negbnlv(double p, int r, double pp, int *x0, int *x1)
{
    int x;
    
	if(!isZeroOne(p)){
		NUMERICS_ERROR("negbnlv", "Invalid p parameter");
	} else if(!isZeroOne(pp)){
		NUMERICS_ERROR("negbnlv", "Invalid pp parameter");
	} else if(isNonPositive(r)){
		NUMERICS_ERROR("negbnlv", "Invalid r parameter");
	} else {
		x = (int)(r*(1.0-p)/p);
		while(negbnlp(x, r, pp) > p) x--;
		while(negbnlp(x, r, pp) <= p) x++;
		if(x == 0){
			*x0 = x;
		} else *x0 = x-1;
		*x1 = x;
	}
}

NUMERICS_EXPORT double poissonp(int x, double l)
{
	if(!isPositive(l)){
		NUMERICS_ERROR("poissonp", "Invalid l parameter");
	}
    return gammq(x + 1.0, l);
}

NUMERICS_EXPORT void poissonv(double p, double l, int *x0, int *x1)
{
    int x;
    
	if(!isZeroOne(p)){
		NUMERICS_ERROR("poissonv", "Invalid p parameter");
	} else if(!isPositive(l)){
		NUMERICS_ERROR("poissonv", "Invalid l parameter");
	} else {
		x = (int)l;
		while(poissonp(x, l) > p) x--;
		while(poissonp(x, l) <= p) x++;
		if(x == 0){
			*x0 = x;
		} else *x0 = x-1;
		*x1 = x;
	}
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
 Version 1.1 - 05/09/1999 - Added use of domain calls.
===================================================================*/
