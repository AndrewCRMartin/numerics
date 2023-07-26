/*===================================================================
 deviate.cpp

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#include <math.h>

#include "contdist.h"
#include "deviate.h"
#include "mathx.h"

/*--------------------------------------------------------------------
 Function to let user specifiy the Uniform(0,1) generator of thier
 choice.  RAND(seed) is then used by the other random variable
 generators.
--------------------------------------------------------------------*/
#define RAND(seed) ( ran1(seed) )

#define IA   16087L
#define IM   2147483647L
#define AM   ( 1.0 / (double)IM ) 
#define IQ   127773L
#define IR   2836L
#define MASK 123459876L

NUMERICS_EXPORT double ran0(long *idum)
{
    long k;
    double ans;
    
    *idum ^= MASK;
    k = (*idum)/IQ;
    *idum = IA*(*idum-k*IQ)-IR*k;
    if(*idum < 0) *idum += IM;
    ans = AM*(*idum);
    *idum ^= MASK;
    
    return ans;
}

#define NTAB 32L
#define NDIV ( 1.0 + (double)(IM-1.0) / (double)NTAB )
#define RNMX ( 1.0 - NUMERICS_MAX_ERROR )

NUMERICS_EXPORT double ran1(long *idum)
{
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    double temp;
    
    if(*idum <= 0 || !iy){
        if(-(*idum) < 1){
            *idum = 1;
        } else {
            *idum = -(*idum);
        }
        for(j = NTAB+7; j >= 0; j--){
            k = (*idum) / IQ;
            *idum = IA * (*idum - k * IQ) - IR * k;
            if(*idum < 0){
                *idum += IM;
            }
            if(j < NTAB){
                iv[j] = *idum;
            }
        }
        iy = iv[0];
    }
    k = (*idum) / IQ;
    *idum = IA * (*idum - k * IQ) - IR * k;
    if(*idum < 0){
        *idum += IM;
    }
    j = (int)(iy / NDIV);
    iy = iv[j];
    iv[j] = *idum;
    if((temp = AM * iy) > RNMX){
        return RNMX;
    } else {
        return temp;
    }
}

#define IA1   40014L
#define IA2   40692L
#define IM1   2147483563L
#define AM1   ( 1.0 / (double)IM1 )
#define IM2   2147483399L
#define IMM1  (IM1 - 1)
#define IQ1   53668L
#define IQ2   52774L
#define IR1   12211L
#define IR2   3791
#define NDIV1 (1.0+ (double)IMM1 / (double)NTAB )

NUMERICS_EXPORT double ran2(long *idum)
{
    int j;
    long k;
    static long idum2 = 123456789L;
    static long iy=0;
    static long iv[NTAB];
    double temp;
    
    if(*idum <= 0){
        if(-(*idum) < 1) *idum = 1;
        else *idum = -(*idum);
        idum2=*idum;
        for(j=NTAB+7;j>=0;j--){
            k=(*idum)/IQ1;
            *idum =IA1*(*idum-k*IQ1)-k*IR1;
            if(*idum < 0) *idum += IM1;
            if(j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if(*idum < 0) *idum += IM1;
    k = idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if(idum2 < 0) idum2 += IM2;
    j = (int)(iy/NDIV1);
    iy=iv[j]-idum2;
    iv[j]=*idum;
    if(iy < 1)iy +=IMM1;
    if((temp=AM1*iy) > RNMX) return RNMX;
    else return temp;
}

#define MBIG  1000000000L
#define FAC   ( 1.0 / (double)MBIG )
#define MSEED 161803398L
#define MZ    0L

NUMERICS_EXPORT double ran3(long *idum)
{
    static int inext, inextp;
    static long ma[56];
    static int iff = 0;
    long mj, mk;
    int i, ii, k;
    
    if(*idum < 0 || iff == 0){
        iff = 1;
        mj = MSEED - (*idum < 0 ? -*idum : *idum);
        mj %= MBIG;
        ma[55] = mj;
        mk = 1;
        for(i = 1; i <= 54; i++){
            ii = (21 * i) % 55;
            ma[ii] = mk;
            mk = mj - mk;
            if(mk < MZ) mk += MBIG;
            mj = ma[ii];
        }
        for(k = 1; k <= 4; k++){
            for(i = 1; i <= 55; i++){
                ma[i] -= ma[1+(i+30) % 55];
                if(ma[i] < MZ) ma[i] += MBIG;
            }
        }
        inext = 0;
        inextp = 31;
        *idum = 1;
    }
    if(++inext == 56) inext = 1;
    if(++inextp == 56) inextp = 1;
    mj = ma[inext] - ma[inextp];
    if(mj < MZ) mj += MBIG;
    ma[inext] = mj;
    
    return mj * FAC;
}

NUMERICS_EXPORT double berdev(double p, long *idum)
{
    if(RAND(idum) <= p) return 1.0;

    return 0.0;
}

NUMERICS_EXPORT double betadev(double a, double b, long *idum)
{
    double aa, bb, c, l, m, s, y, x;
    
    if(a <= 1.0 || b <= 1.0) return betav(RAND(idum), a, b);
    aa = a-1.0; bb = b-1.0; c = aa + bb; l = c * log(c); m = aa / c; s = .5 / sqrt(c);
    do {
        do {
            y = gasdev(0.0, 1.0, idum);
            x = s * y + m;
        } while(x < 0.0 || x > 1.0);
    } while(log(RAND(idum)) > aa * log(x/aa) + bb*log((1.0-x)/bb) + l + .5*y*y);
    
    return x;
}

NUMERICS_EXPORT double bnldev(int n, double pp, long *idum)
{
    int i;
    static int nold = -1;
    double am, em, g, angle, p, bnl, sq, t, y;
    static double pold = -1.0, pc, plog, pclog, en, oldg;
    
    p = (pp <= .5 ? pp : 1.0 - pp);
    am = n * p;
    if(n < 25){
        bnl = 0.0;
        for(i = 0; i < n; i++) bnl += berdev(p, idum);
    } else if(am < 1.0){
        g = exp(-am);
        t = 1.0;
        for(i = 0; i <= n; i++){
            t *= RAND(idum);
            if(t < g) break;
        }
        bnl = (i <= n ? i : n);
    } else {
        if(n != nold){
            en = n;
            oldg = gammln(en + 1.0);
            nold = n;
        }
        if(p != pold){
            pc = 1.0 - p;
            plog = log(p);
            pclog = log(pc);
            pold = p;
        }
        sq = sqrt(2.0 * am * pc);
        do {
            do {
                angle = NUMERICS_PI * RAND(idum);
                y = tan(angle);
                em = sq * y + am;
            } while(em < 0.0 || em >= (en + 1.0));
            em = floor(em);
            t = 1.2 * sq * (1.0 + y * y) * exp(oldg - gammln(em + 1.0) - gammln(en - em + 1.0) + em * plog + (en - em) * pclog);
        } while(RAND(idum) > t);
        bnl = em;
    }
    if(p != pp) bnl = n - bnl;
    
    return bnl;
}

NUMERICS_EXPORT double caudev(double m, double s, long *idum)
{
    return m + s*tan(NUMERICS_PI*(RAND(idum)-.5));
}

NUMERICS_EXPORT double chidev(double v, long *idum)
{
    return gamdev(v * .5, 2.0, idum);
}

NUMERICS_EXPORT double dexpdev(double m, double s, long *idum)
{
    double u = RAND(idum), sn = sign(2.0 * u - 1.0);
    
    return m - s * sn * log(1.0-sn*(2.0*u-1.0));
}

NUMERICS_EXPORT double expdev(double b, long *idum)
{
    double dum;
    
    do
        dum=RAND(idum);
    while(dum == 0.0);
    return -log(dum)*b;
}

NUMERICS_EXPORT double fdev(double dfn, double dfd, long *idum)
{
    double xnum = chidev(dfn, idum) / dfn, xden;
    
    do {
        xden = chidev(dfd, idum) / dfd;
    } while(xden <= 9.999999999998E-39*xnum);
    
    return (xnum / xden);
}

NUMERICS_EXPORT double gamdev(double a, double b, long *idum)
{
    double u, bb, p, x, y;
    
    if(a < 1.0){
        do {
            u = RAND(idum);
            bb = (NUMERICS_E + a) / NUMERICS_E;
            if((p = bb * u) > .1){
                x = -log((bb-p) / a);
                if(RAND(idum) <= pow(x, a - 1.0)) return b * x;
            } else {
                x = pow(p, 1.0/a);
                if(RAND(idum) <= exp(-x)) return b * x;
            }
        } while(1);
    } else {
        do {
            y = expdev(1.0, idum);
        } while(RAND(idum) > pow(y / exp(y-1.0), a - 1.0));
        
        return a * b * y;
    }
}

NUMERICS_EXPORT double gasdev(double m, double v, long *idum)
{
    static int iset=0;
    static double gset;
    double fac, rsq, v1, v2;
    
    if(iset==0){
        do{
            v1=2.0*RAND(idum)-1.0;
            v2=2.0*RAND(idum)-1.0;
            rsq=v1*v1+v2*v2;
        } while(rsq>=1.0 || rsq==0.0);
        fac=sqrt(-2.0*log(rsq)/rsq);
        gset=v1*fac;
        iset=1;
        return v2*fac*sqrt(v)+m;
    } else {
        iset = 0;
        return gset*sqrt(v)+m;
    }
}

NUMERICS_EXPORT double logdev(double m, double b, long *idum)
{
    return m-b*log(1.0/RAND(idum)-1.0);
}

NUMERICS_EXPORT double nchidev(double df, double xnonc, long *idum)
{
    double gas = gasdev(sqrt(xnonc), 1.0, idum);
    
    return (chidev(df-1.0, idum) + gas * gas);
}

NUMERICS_EXPORT double nfdev(double dfn, double dfd, double xnonc, long *idum)
{
    double xden, xnum;

    do {
        xnum = nchidev(dfn, xnonc, idum) / dfn;
        xden = chidev(dfd, idum) / dfd;
    } while(xden < 1.0 && xnum > NUMERICS_FLOAT_MAX * xden);
    
    return (xnum/xden);
}

NUMERICS_EXPORT void orddev(double *x, int n, long *idum)
{
    double s = 0.0;
    int i;
    
    for(i = 0; i < n; i++) x[i] = (s += expdev(1.0, idum));
    s += expdev(1.0, idum);
    for(i = 0; i < n; i++) x[i] /= s;
}

NUMERICS_EXPORT double poisdev(double xm, long *idum)
{
    static double sq, alxm, g, oldm = -1.0;
    double em, t, y;
    
    if(xm < 12.0){
        if(xm != oldm){
            oldm = xm;
            g = exp(-xm);
        }
        em = -1;
        t = 1.0;
        do {
            ++em;
            t *= RAND(idum);
        } while(t > g);
    } else {
        if(xm != oldm){
            oldm = xm;
            sq = sqrt(2.0 * xm);
            alxm = log(xm);
            g = xm * alxm - gammln(xm + 1.0);
        }
        do {
            do {
                y = tan(NUMERICS_PI * RAND(idum));
                em = sq * y + xm;
            } while (em < 0.0);
            em = floor(em);
            t = .9 * (1.0 + y*y) * exp(em *alxm-gammln(em+1.0)-g);
        } while(RAND(idum) > t);
    }
    
    return em;
}

double funt(double x, int v)
/* internal function needed for tdev(v, idum) */
{
    return pow(1.0 + x * x / v, -(v + 1.0) * .5);
}

NUMERICS_EXPORT double tdev(int v, long *idum)
{
    double c, u1, u2, u3, x, s;
    
    c = (v + 1.0) * .5 / (sqrt(NUMERICS_PI * v) * v * .5);
    u1 = RAND(idum);
    if(u1 >= sqrt(2.0 / NUMERICS_PI) || u1 >= 2.0 * c){
        do {
            u2 = RAND(idum);
            if(u2 <= .3622520694){
                if(u2 <= .05300969080){
                    s = sign(7.840088159 * u2 - .2078);
                    x = s * (fabs(7.840088159 * u2 - .2078) + 1.7922);
                    u3 = RAND(idum);
                    if(.2 * u3 <= funt(x, v) - 1.0 + fabs(x) * .5) return x;
                } else {
                    x = 11.5909050257 * u2 - 2.406629332;
                    u3 = RAND(idum);
                    if(.13528 * u3 <= funt(x, v) - 1.0 + fabs(x) * .5) return x;
                }
            } else {
                x = 1.0 / (1.0680176321 - 1.5680176321 * u2);
                u3 = RAND(idum);
                if(u3 <= x * x * funt(x, v)) return x;
            }
        } while(1);
    }
    
    return 2.0 * (RAND(idum) + RAND(idum) - 1.0);
}

NUMERICS_EXPORT double unifdev(double lo, double hi, long *idum)
{
    return (lo + RAND(idum) * (hi - lo));
}

NUMERICS_EXPORT double weibdev(double a, double b, long *idum)
{
    return b*pow(log(1.0/RAND(idum)),1.0/a);
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/

