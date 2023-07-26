/*===================================================================
 rankdist.cpp

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#include <math.h>

#include "normdist.h"
#include "rankdist.h"
#include "vecutils.h"

int wilcox1(int x, int n, int k)
{
    int c, h, m, w, v, j, y;
    
    if(k == 1){
        if(x < n) return x;
        return n;
    }
    c = 1;
    h = k-1;
    m = n-1;
    w = 0;
    for(j = 1; j <= h; j++) c = c*(m-j+1)/j;
    v = h*(m-h)+h*(h+1)/2;
    y=x-h-1;
    while(y >= v){
        w += c;
        c = (c*(m-h))/m;
        m--;
        v-=h;
        y=y-h-1;
    }
    do {
        w += wilcox1(y, m, h);
        m--;
        y = y-h-1;
    } while(2*y >= h*(h+1));
    
    return w;
}

NUMERICS_EXPORT double wilcoxonp(int x, int n)
{
    int c, k, u, v, w, xx;
    BOOL a;
    double m=n*(n+1.0)/4.0, p, cv;
    
    if(n <= 25){
        if(x < 0) return 0.0;
        m = n*(n+1)/2;
        if(x >= m) return 1.0;
        if(2*x > m){
            xx = (int)(m)-x-1;
            a = FALSE;
        } else {
            xx = x;
            a = TRUE;
        }
        c = 1; k = 1; u = 1; v = n; w = 1;
        while(xx >= v){
            c = (c*(n-k+1))/k;
            w += c;
            k++;
            u += k;
            v += n+1-k;
        }
        do {
            w += wilcox1(xx, n, k);
            k++;
            u += k;
        } while(xx >= u);
        if(a) p = ldexp((double)w, -(int)(n));
        else p = 1.0 - ldexp((double)w, -(int)(n));
    } else {
        cv = n*(n+1.0)*(2.0*n+1.0)/24.0;
        p = normalp(x+.5, m, cv);
    }

    return p;
}

NUMERICS_EXPORT void wilcoxonv(double p, int n, int *w0, int *w1)
{
    double m, v, wd;
    int w;
    
    m = (double)(n*(n+1))/4.0;
    v = (double)(n*(n+1)*(2.0*n+1))/24.0;
    wd = normalv(p, m, v);
    if(p < 0.0 || p > 1.0 || n < 1){
        *w0 = *w1 = -1;
    } else if(n > 25){
        if(wd < 1.0) wd = 1.0;
        else if(wd > 2.0*m) wd = 2.0*m-.5;
        *w0 = (int)floor(wd);
        *w1 = (int)ceil(wd);
    } else {
        w = (int)wd;
        while(wilcoxonp(w, n) > p) w--;
        while(wilcoxonp(w, n) <= p) w++;
        if(w <= 1){
            *w0 = w;
        } else *w0 = w-1;
        *w1 = w;
    }
}

NUMERICS_EXPORT double wilmannwhitp(int w, int m, int n)
{
    int u = w - m * (m + 1) / 2;
    int minmn = ((m < n) ? m : n);
    int mn1 = m*n+1, i;
    int maxmn = ((m > n) ? m : n), n1 = maxmn + 1;
    double ret = 0.0;
    double *freq;
    double *work;
    
    if(u < 0){
        return ret;
    }
    
    if(minmn < 1){
        return ret;
    }

    freq = dvector(0, mn1 - 1);
    
    for(i = 1; i <= n1; ++i){
        freq[i-1] = 1.0;
    }
    
    if(minmn != 1){
        double sum;
        int in = maxmn, l, k, j;
        
        work = dvector(0, (mn1 + 1) / 2 + minmn - 1);
        
        n1 = n1 + 1;
        for(i = n1; i <= mn1; ++i){
            freq[i-1] = 0.0;
        }
        
        work[1-1] = 0.0;
        for(i = 2; i <= minmn; ++i){
            work[i-1] = 0.0;
            in += maxmn;
            n1 = in + 2;
            l = 1 + in/2;
            k = i;
            for(j = 1; j <= l; ++j){
                ++k;
                --n1;
                sum = freq[j-1] + work[j-1];
                freq[j-1] = sum;
                work[k-1] = sum - freq[n1-1];
                freq[n1-1] = sum;
            }
        }
        sum = 0.0;
        for(i = 1; i <= mn1; ++i){
            sum += freq[i-1];
            freq[i-1] = sum;
        }
        for(i = 1; i <= u+1; ++i){
            freq[i-1] /= sum;
        }
        ret = freq[u];
        
        free_dvector(work, 0);
    }
    
    free_dvector(freq, 0);

    return ret;
}

NUMERICS_EXPORT void wilmannwhitv(double p, int m, int n, int *w0, int *w1)
{
    double mn, v, wd;
    int w;
    
    if(p < 0.0 || p > 1.0 || n < 1){
        *w0 = *w1 = -1;
    } else {
        mn = (double)(m*(m+n+1))/2.0;
        v = (double)(m*n*(m+n+1))/12.0;
        wd = normalv(p, mn, v);
        w = (int)wd;
        while(wilmannwhitp(w, m, n) > p) w--;
        while(wilmannwhitp(w, m, n) <= p) w++;
        if(w <= m*(m+1)/2){
            *w0 = m*(m+1)/2;
        } else {
            *w0 = w-1;
        }
        *w1 = w;
    }
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/
