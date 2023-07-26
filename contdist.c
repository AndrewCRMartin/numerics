/*===================================================================
 contdist.c

 Version 1.1

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  April 10, 1999
===================================================================*/

#include <math.h>

#include "contdist.h"
#include "domain.h"
#include "mathx.h"
#include "numerror.h"
#include "utility.h"

NUMERICS_EXPORT double betap(double x, double a, double b)
{
    double ret = -1.0;

    if(!isPositive(a)){
        NUMERICS_ERROR("betap", "Invalid a parameter");
    } else if(!isPositive(b)){
        NUMERICS_ERROR("betap", "Invalid b parameter");
    } else if(!isPositive(x)){
        ret = 0.0;
    } else if(x >= 1.0){
        ret = 1.0;
    } else {
        ret = betai(a, b, x);
    }
    
    return ret;
}

NUMERICS_EXPORT double betav(double p, double a, double b)
{
    double b0, b1, b2;
    int i;
    double ret = -1.0;
    
    if(!isZeroOne(p)){
        NUMERICS_ERROR("betav", "Invalid p parameter");
	} else if(!isPositive(a)){
		NUMERICS_ERROR("betav", "Invalid a parameter");
    } else if(!isPositive(b)){
		NUMERICS_ERROR("betav", "Invalid b parameter");
    } else {
        b0 = b1 = .5;
        while(betap(b0, a, b) > p){
            b0 -= .1;
        }
        while(betap(b1, a, b) < p){
            b1 += .1;
        }
        if(b0 < 0.0){
            b0 = 0.0;
        }
        if(b1 > 1.0){
            b1 = 1.0;
        }
        
        i = 0;
        do {
            b2 = b0 + (b1 - b0) / 2.0;
            if(betap(b2, a, b) > p){
                b1 = b2;
            } else {
                b0 = b2;
            }
        } while(!areEqual(b0, b1, NUMERICS_MAX_ERROR) && ++i < NUMERICS_ITMAX);
        
        if(i >= NUMERICS_ITMAX){
            NUMERICS_ERROR("betav", "Iteration failed to converge");
        } else {
            ret = b2;
        }
    }
    
    return ret;
}

NUMERICS_EXPORT double cauchyp(double x, double m, double s)
{
    double ret = -1.0;
    
    if(!isPositive(s)){
        NUMERICS_ERROR("cauchyp", "Invalid s parameter");
    } else {
        ret = atan((x - m) / s) / NUMERICS_PI + .5;
    }
    
    return ret;
}

NUMERICS_EXPORT double cauchyv(double p, double m, double s)
{
    double ret = 0.0;
    
    if(!isPositive(s)){
        NUMERICS_ERROR("cauchyv", "Invalid s parameter");
    } else if(!isZeroOne(p)){
        NUMERICS_ERROR("cauchyv", "Invalid p parameter");
    } else {
        ret = m + s * tan(NUMERICS_PI * (p - .5));
    }
    
    return ret;
}

NUMERICS_EXPORT double dblexpp(double x, double m, double s)
{
    double ret = 0.0;
    double sn = sign(x - m);
    
    if(!isPositive(s)){
        NUMERICS_ERROR("dblexpp", "Invalid s parameter");
    } else {
        ret = (1.0 + sn * (1.0 - exp(-sn * (x - m) / s))) * .5;
    }
    
    return ret;
}

NUMERICS_EXPORT double dblexpv(double p, double m, double s)
{
    double sn = sign(2.0 * p - 1.0);
    double ret = 0.0;
    
    if(!isPositive(s)){
        NUMERICS_ERROR("dblexpv", "Invalid s parameter");
    } else if(!isZeroOne(p)){
        NUMERICS_ERROR("dblexpv", "Invalid p parameter");
    } else {
        ret = m - s * sn * log(1.0 - sn * (2.0 * p - 1.0));
    }
    
    return ret;
}

NUMERICS_EXPORT double expp(double x, double b)
{
    double ret = -1.0;
    
    if(!isPositive(b)){
        NUMERICS_ERROR("expp", "Invalid b parameter");
    } else if(!isPositive(x)){
        ret = 0.0;
    } else {
        ret = 1.0 - exp(-x / b);
    }
    
    return ret;
}

NUMERICS_EXPORT double expv(double p, double b)
{
    double ret = -1.0;
    
    if(!isZeroOne(p)){
        NUMERICS_ERROR("expv", "Invalid p parameter");
    } else if(!isPositive(b)){
        NUMERICS_ERROR("expv", "Invalid b parameter");
    } else {
        ret = -log(1 - p) * b;
    }
    
    return ret;
}

NUMERICS_EXPORT double gammap(double x, double a, double b)
{
    double ret = -1.0;

    if(!isPositive(a)){
        NUMERICS_ERROR("gammp", "Invalid a parameter");
    } else if(!isPositive(b)){
        NUMERICS_ERROR("gammp", "Invalid b parameter");
    } else if(!isPositive(x)){
        ret = 0.0;
    } else {
        ret = 1.0 - gammq(a, x / b);
    }
    
    return ret;
}

NUMERICS_EXPORT double gammav(double p, double a, double b)
{
    int i;
    double g0, g1, g2;
    double ret = -1.0;
    
    if(!isPositive(a)){
        NUMERICS_ERROR("gammav", "Invalid a parameter");
    } else if(!isPositive(b)){
        NUMERICS_ERROR("gammav", "Invalid b parameter");
    } else if(!isZeroOne(p)){
        NUMERICS_ERROR("gammav", "Invalid p parameter");
    } else {
        g0 = g1 = a*b;
        while(gammap(g0, a, b) > p){
            g0 -= .5;
        }
        while(gammap(g1, a, b) < p){
            g1 += .5;
        }
        if(g0 < 0.0){
            g0 = 0.0;
        }
        i = 0;
        do {
            g2 = g0 + (g1-g0)/2.0;
            if(gammap(g2,a,b) > p){
                g1 = g2;
            } else {
                g0 = g2;
            }
        } while(!areEqual(g0, g1, NUMERICS_MAX_ERROR) && ++i < 100);
        if(i >= NUMERICS_ITMAX){
            NUMERICS_ERROR("gammav", "Iteration failed to converge");
        } else {
            ret = g2;
        }
    }
    
    return ret;
}

NUMERICS_EXPORT double logisticp(double x, double m, double b)
{
    double ret = -1.0;
    
    if(!isPositive(b)){
        NUMERICS_ERROR("logisticp", "Invalid b parameter");
    } else {
        ret = 1.0 / (1.0 + exp((m - x) / b));
    }
    
    return ret;
}

NUMERICS_EXPORT double logisticv(double p, double m, double b)
{
    double ret = -1.0;
    
    if(!isPositive(b)){
        NUMERICS_ERROR("logisticp", "Invalid b parameter");
    } else if(!isZeroOne(p)){
        NUMERICS_ERROR("logisticv", "Invalid p parameter");
    } else {
        ret = m - b * log(1.0 / p - 1.0);
    }
    
    return ret;
}

NUMERICS_EXPORT double lognormp(double x, double m, double v)
{
    double ret = -1.0;
    
    if(!isPositive(v)){
        NUMERICS_ERROR("lognormp", "Invalid v parameter");
    } else if(!isPositive(x)){
        ret = 0.0;
    } else {
        ret = .5*(erff((log(x)-m)/(sqrt(2.0*v)))+1.0);
    }
    
    return ret;
}

NUMERICS_EXPORT double lognormv(double p, double m, double v)
{
    double l0, l1, l2;
    int i;
    double ret = -1.0;
    
    if(!isPositive(v)){
        NUMERICS_ERROR("lognormv", "Invalid v parameter");
    } else if(!isZeroOne(p)){
        NUMERICS_ERROR("lognormv", "Invalid p parameter");
    } else {
        l0 = l1 = exp(m+v*v*.5);
        while(lognormp(l0, m, v) > p){
            l0 -= .5;
        }
        while(lognormp(l1, m, v) < p){
            l1 += .5;
        }
        if(l0 < 0.0){
            l0 = 0.0;
        }
        i = 0;
        do {
            l2 = l0 + (l1-l0)/2.0;
            if(lognormp(l2,m,v) > p){
                l1 = l2;
            } else {
                l0 = l2;
            }
        } while(!areEqual(l0, l1, NUMERICS_MAX_ERROR) && ++i < 100);
        if(i >= NUMERICS_ITMAX){
            NUMERICS_ERROR("lognormv", "Iteration failed to converge");
        } else {
            ret = l2;
        }
    }
    
    return l2;
}

NUMERICS_EXPORT double rayleighp(double x, double s)
{
    return 1.0 - exp(-x * x/ (2.0 * s * s));
}

NUMERICS_EXPORT double rayleighv(double p, double s)
{
    double ret = -1.0;
    
    if(!isZeroOne(p)){
        NUMERICS_ERROR("rayleighv", "Invalid p parameter");
    } else {
        ret = s * sqrt(log( 1.0 / ((1.0 - p) * (1.0 - p))));
    }
    
    return ret;
}

NUMERICS_EXPORT double uniformp(double u, double a, double b)
{
    double ret = -1.0;
    
    if(a > b){
        NUMERICS_ERROR("uniformp", "Invalid a and b parameters");
    } else if(u <= a){
        ret = 0.0;
    } else if(u > b){
        ret = 1.0;
    } else {
        ret = (u - a) / (b - a);
    }
    
    return ret;
}

NUMERICS_EXPORT double uniformv(double p, double a, double b)
{
    double ret = -1.0;
    
    if(!isZeroOne(p)){
        NUMERICS_ERROR("uniformv", "Invalid p parameter");
    } else if(a > b){
        NUMERICS_ERROR("uniformv", "Invalid a and b parameters");
    } else {
        ret = p * (b - a) + a;
    }
    
    return ret;
}

NUMERICS_EXPORT double weibullp(double x, double g, double b)
{
    double ret = -1.0;
    
    if(!isPositive(g)){
        NUMERICS_ERROR("weibullp", "Invalid g parameter");
    } else if(!isPositive(b)){
        NUMERICS_ERROR("weibullp", "Invalid b parameter");
    } else if(!isPositive(x)){
        ret = 0.0;
    } else {
        ret = 1.0 - exp(-pow(x / b, g));
    }
    
    return ret;
}

NUMERICS_EXPORT double weibullv(double p, double g, double b)
{
    double ret = -1.0;
    
    if(!isPositive(g)){
        NUMERICS_ERROR("weibullv", "Invalid g parameter");
    } else if(!isPositive(b)){
        NUMERICS_ERROR("weibullv", "Invalid b parameter");
    } else if(!isZeroOne(p)){
        NUMERICS_ERROR("weibullv", "Invalid p parameter");
    } else {
        ret = b * pow(-log(1.0 - p), (1.0 / g));
    }
    
    return ret;
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
 Version 1.1 - 04/10/1999 - Use areEqual to test for floating-point
                            equality.
                            Added use of domain calls.
===================================================================*/
