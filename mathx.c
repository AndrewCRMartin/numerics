/*===================================================================
 mathx.h

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

#include "mathx.h"
#include "numerror.h"

NUMERICS_EXPORT double acosh(double x)
{
    return ((x <= 1.0) ? 0.0 :
            ((x > 1.0e10) ? 0.69314718055995 + log(x) :
             log(x + sqrt((x-1.0)*(x+1.0)))));
}

NUMERICS_EXPORT double asinh(double x)
{
    double ax = fabs(x);
    if(ax > 1.0e10){
        return ((x > 0.0) ? 0.69314718055995 + log(ax) :
                -.69314718055995 + log(ax));
    } else {
        double y = x*x;
        return ((x == 0.0) ? 0.0 :
                ((x > 0.0) ? log1x(ax+y/(1.0+sqrt(1.0+y))) :
                 -log1x(ax+y/(1.0+sqrt(1.0+y)))));
    }
}

NUMERICS_EXPORT double atanh(double x)
{
    double ax = fabs(x);
    
    if(ax >= 1.0){
        return ((x > 0.0) ? FLT_MAX : -FLT_MAX);
    } else {
        return ((x == 0.0) ? 0.0 :
                ((x > 0.0) ? .0*log1x(2.0*ax/(1.0-ax)) :
                 -.5*log1x(2.0*ax/(1.0-ax))));
    }
}

NUMERICS_EXPORT double log1x(double x)
{
    if(x == 0.0){
        return 0.0;
    } else if(x < -.2928 || x > .4142){
        return log(1.0+x);
    } else {
        double z = x/(x+2.0);
        double y = z*z;
        return (z*(2.0+y*(0.666666666663366+y*(0.400000001206045+y*
               (0.285714091590488+y*(0.22223823332791+y*
               (0.1811136267967+y*0.16948212488)))))));
    }
}

NUMERICS_EXPORT double betacf(double a, double b, double x)
{
    int m,m2;
    double aa,c=1.0,del,h,qab=a+b,qam=a-1.0,qap=a+1.0, d=1.0-qab*x/qap;
    
    if(fabs(d) < NUMERICS_FLOAT_MIN) d=NUMERICS_FLOAT_MIN;
    d=1.0/d;
    h=d;
    for(m=1; m <= NUMERICS_ITMAX; m++){
        m2=2*m;
        aa=m*(b-m)*x/((qam+m2)*(a+m2));
        d=1.0+aa*d;
        if(fabs(d) < NUMERICS_FLOAT_MIN) d = NUMERICS_FLOAT_MIN;
        c=1.0+aa/c;
        if(fabs(c) < NUMERICS_FLOAT_MIN) c = NUMERICS_FLOAT_MIN;
        d=1.0/d;
        h*=d*c;
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
        d=1.0+aa*d;
        if(fabs(d) < NUMERICS_FLOAT_MIN) d = NUMERICS_FLOAT_MIN;
        c=1.0+aa/c;
        if(fabs(c) < NUMERICS_FLOAT_MIN) c = NUMERICS_FLOAT_MIN;
        d=1.0/d;
        del=d*c;
        h*=del;
        if(fabs(del-1.0) < NUMERICS_MAX_ERROR) break;
    }
    if(m > NUMERICS_ITMAX){
        NUMERICS_ERROR("betacf", "a or b too big, or NUMERICS_ITMAX too small");
    }
    return h;
}

void gser(double *gamser, double a, double x, double *gln)
{
    int n;
    double sum, del, ap;
    
    *gln = gammln(a);
    if(x <= 0.0){
        if(x < 0.0){
            NUMERICS_ERROR("gser", "x less than 0");
        }
        *gamser=0.0;
        return;
    } else {
        ap = a;
        del = sum = 1.0/a;
        for(n = 1; n <=NUMERICS_ITMAX; n++){
            ++ap;
            del *= x / ap;
            sum += del;
            if(fabs(del) < fabs(sum)*NUMERICS_MAX_ERROR){
                *gamser = sum * exp(-x+a*log(x)-(*gln));
                return;
            }
        }
    }
    NUMERICS_ERROR("gser", "a too large, NUMERICS_ITMAX too small");
    return;
}

void gcf(double *gammcf, double a, double x, double *gln)
{
    int i;
    double an, b, c, d, del, h;
    
    *gln = gammln(a);
    b = x + 1.0 - a;
    c = 1.0 / NUMERICS_FLOAT_MIN;
    d = 1.0 / b;
    h = d;
    for(i = 1; i <= NUMERICS_ITMAX; i++){
        an = -i*(i-a);
        b += 2.0;
        d = an*d+b;
        if(fabs(d) < NUMERICS_FLOAT_MIN) d = NUMERICS_FLOAT_MIN;
        c = b + an / c;
        if(fabs(c) < NUMERICS_FLOAT_MIN) c = NUMERICS_FLOAT_MIN;
        d = 1.0 / d;
        del = d*c;
        h *= del;
        if(fabs(del-1.0) < NUMERICS_MAX_ERROR) break;
    }
    if(i > NUMERICS_ITMAX){
        NUMERICS_ERROR("gcf", "a too large, NUMERICS_ITMAX too small");
    }
    *gammcf = exp(-x+a*log(x)-(*gln))*h;
}

NUMERICS_EXPORT double factrl(int n)
{
    static int ntop = 4;
    static double a[33] = {1.0, 1.0, 2.0, 6.0, 24.0};
    int j;
    
    if(n < 0){
        NUMERICS_ERROR("factrl", "Negative factorial");
        return 0.0;
    }
    if(n > 32) return exp(gammln(n + 1.0));
    while(ntop < n){
        j = ntop++;
        a[ntop] = a[j] * ntop;
    }
    return a[n];
}

NUMERICS_EXPORT double betai(double a, double b, double x)
{
    double bt;
    
    if(x < 0.0 || x > 1.0){
        NUMERICS_ERROR("betai", "Bad x");
        return 0.0;
    }
    if(x == 0.0 || x == 1.0) bt = 0.0;
    else bt = exp(gammln(a + b) - gammln(a) - gammln(b) + a * log(x) + b *
                  log(1.0 - x));
    if(x < (a + 1.0) / (a + b + 2.0)) return bt * betacf(a, b, x) / a;
    else return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
}

NUMERICS_EXPORT double gammln(double xx)
{
    double x, y, tmp, ser;
    static double cof[6]={76.18009172947146, -86.50532032941677,
        24.01409824083091, -1.231739572460166, 0.1208650973866179e-2,
        -0.5395239384953e-5};
    int j;

    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for(j = 0; j <= 5; j++) ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005 * ser / x);
}

NUMERICS_EXPORT double gammp(double a, double x)
{
    double gamser, gammcf, gln;
    
    if(x < 0.0 || a <= 0.0){
        NUMERICS_ERROR("gammp", "Invalid arguments");
        return 0.0;
    }
    if(x < (a+1.0)){
        gser(&gamser, a, x, &gln);
        return gamser;
    } else {
        gcf(&gammcf, a, x, &gln);
        return 1.0 - gammcf;
    }
}

NUMERICS_EXPORT double gammq(double a, double x)
{
    double gamser, gammcf, gln;
    
    if(x < 0.0 || a <= 0.0){
        NUMERICS_ERROR("gammq", "Invalid arguments");
        return 0.0;
    }
    if(x < (a+1.0)){
        gser(&gamser, a, x, &gln);
        return 1.0 - gamser;
    } else {
        gcf(&gammcf, a, x, &gln);
        return gammcf;
    }
}

NUMERICS_EXPORT double beta(double z, double w)
{
    return exp(gammln(z) + gammln(w) - gammln(z+w));
}

NUMERICS_EXPORT double erff(double x)
{
    return x < 0.0 ? -gammp(.5, x * x) : gammp(.5, x * x);
}

NUMERICS_EXPORT double bico(int n, int k)
{
    return floor(0.5 + exp(factln(n) - factln(k) - factln(n - k)));
}

NUMERICS_EXPORT double bicoln(int n, int k)
{
    return factln(n) - factln(k) - factln(n - k);
}

NUMERICS_EXPORT double factln(int n)
{
    static double a[101];
    
    if(n < 0){
        NUMERICS_ERROR("factrl", "Negative factorial");
        return 0.0;
    }
    if(n <= 1) return 0.0;
    if(n <= 100) return a[n] ? a[n] : (a[n] = gammln(n + 1.0));
    else return gammln(n + 1.0);
}

NUMERICS_EXPORT double erffc(double x)
{
    return x < 0.0 ? 1.0 + gammp(0.5, x*x) : gammq(0.5, x*x);
}

NUMERICS_EXPORT double erfcc(double x)
{
    double t,z, ans;

    z = fabs(x);
    t=1.0/(1.0+.5*z);
    ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*
        (.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t*
        (1.48851587+t*(-.82215223+t*.17087277)))))))));
    
    return x >= 0.0 ? ans : 2.0-ans;
}

NUMERICS_EXPORT double expint(int n, double x)
{
    int i, ii, nm1;
    double a,b,c,d,del,fact,h,psi,ans;
    
    nm1 = n-1;
    if(n<0||x<0.0||(x==0.0 &&(n==0||n==1))) return 0.0;  /*Error condition*/
    else {
        if(n==0)ans=exp(-x)/x;
        else {
            if(x==0.0) ans=1.0/nm1;
            else {
                if(x > 1.0){
                    b=x+n;
                    c=1.0/NUMERICS_FLOAT_MIN;
                    d=1.0/b;
                    h=d;
                    for(i=1;i<=NUMERICS_ITMAX;i++){
                        a = -i*(nm1+1);
                        b+=2.0;
                        d=1.0/(a*d+b);
                        c=b+a/c;
                        del =c*d;
                        h*=del;
                        if(fabs(del-1.0) < NUMERICS_MAX_ERROR){
                            ans=h*exp(-x);
                            return ans;
                        }
                    }
                    return ans;  /* Error condition */
                } else {
                    ans =(nm1!=0 ? 1.0/nm1 : -log(x)-NUMERICS_EULER);
                    fact=1.0;
                    for(i=1;i<=NUMERICS_ITMAX;i++){
                        fact*=-x/i;
                        if(i!=nm1) del=-fact/(i-nm1);
                        else {
                            psi = -NUMERICS_EULER;
                            for(ii=1;ii<=nm1;ii++) psi += 1.0/ii;
                            del=fact*(-log(x)+psi);
                        }
                        ans+=del;
                        if(fabs(del)<fabs(ans)*NUMERICS_MAX_ERROR) return ans;
                    }
                    return ans;  /*Error condition*/
                }
            }
        }
    }
    return ans;
}

NUMERICS_EXPORT double pythag(double a, double b)
{
    double aa = fabs(a), ab = fabs(b);
    
    if(aa > ab) return (aa * sqrt(1.0 + (ab/aa * ab/aa)));
    else return (ab == 0.0 ? 0.0 : ab * sqrt(1.0 + (aa/ab * aa/ab)));
}

NUMERICS_EXPORT double sign(double x)
{
    if(x < 0.0) return -1.0;
    if(x > 0.0) return 1.0;
    
    return 0.0;
}

/*====================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
====================================================================*/
