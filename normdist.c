/*===================================================================
 normdist.cpp

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
#include "normdist.h"
#include "numerror.h"
#include "utility.h"

NUMERICS_EXPORT double chisqp(double x2, long v)
{
	if(!isPositive(v)){
		NUMERICS_ERROR("chisqp", "Invalid v parameter");
	} else if(x2 <= 0.0){
		return 0.0;
	}
    return gammp(.5 * v, .5 * x2);
}

NUMERICS_EXPORT double chisqv(double p, long v)
{
	double ch = 0.0;

	if(!isZeroOne(p)){
		NUMERICS_ERROR("chisqv", "Invalid p parameter");
	} else if(!isPositive(v)){
		NUMERICS_ERROR("chisqv", "Invalid v parameter");
	} else {
	    double g = gammln(.5*v), aa = .6931471806, e = 5e-7, r1, a, b, c, q, t, x;
		double p1, p2, s1, s2, s3, s4, s5, s6, xx;
		long i;
		xx = .5 * v;
		c = xx - 1.0;
		if(v < -1.24 * log(p)){
			ch = pow(p * xx * exp(g + xx * aa), 1.0 / xx);
			if (ch < e) return ch;
		} else {
			if(v < .32){
				ch = .4;
				a = log(1.0 - p);
				do {
					q = ch;
					p1 = 1.0 + ch * (4.67 + ch);
					p2 = ch * (6.73 + ch * (6.66 + ch));
					t = -.5 + (4.67 + 2.0 * ch) / p1 - (6.73 + ch * (13.32 + 3.0 * ch)) / p2;
					ch -= (1.0 - exp(a + g + .5 * ch + c * aa) * p2 / p1) / t;
				} while(!areEqual(q / ch, 1.0, NUMERICS_MAX_ERROR));
			} else {
				x = normalv(p, 0.0, 1.0);
				p1 = .222222 / v;
				r1 = x * sqrt(p1) + 1.0 - p1;
				ch = v * r1 * r1 * r1;
				if(ch > 2.2 * v + 6.0) ch = -2.0 * (log(1.0 - p) - c * log(.5 * ch) + g);
			}
		}
		i = 0;
		do {
			q = ch;
			p1 = .5 * ch;
			p2 = p - gammp(xx, p1);
			t = p2 * exp(xx * aa + g + p1 - c * log(ch));
			b = t / ch;
			a = .5 * t - b * c;
			s1 = (210.0 + a * (140.0 + a * (105.0 + a * (84.0 + a * (70.0 + 60.0 * a))))) / 420.0;
			s2 = (420.0 + a * (735.0 + a * (966.0 + a * (1141.0 + 1278.0 * a)))) / 2520.0;
			s3 = (210.0 + a * (462.0 + a * (707.0 + 932.0 * a))) / 2520.0;
			s4 = (252.0 + a * (672.0 + 1182.0 * a) + c * (294.0 + a * (889.0 + 1740.0 * a))) / 5040.0;
			s5 = (84.0 + 2264.0 * a + c * (1175.0 + 606.0 * a)) / 2520.0;
			s6 = (120. + c * (346.0 + 127.0 * c)) / 5040.0;
			ch += t * (1.0 + .5 * t * s1 - b * c * (s1 - b * (s2 - b * (s3 - b * (s4 - b * (s5 - b * s6))))));
	    } while(!areEqual(q / ch, 1.0, NUMERICS_MAX_ERROR) && i++ < 20);
	}
    return ch;
}

NUMERICS_EXPORT double fratiop(double f, long v1, long v2)
{
	if(!isPositive(v1)){
		NUMERICS_ERROR("fratiop", "Invalid v1 parameter");
	} else if(!isPositive(v2)){
		NUMERICS_ERROR("fratiop", "Invalid v2 parameter");
	} else if(isNonPositive(f)){
		return 0.0;
	}
    return 1.0 - betai(v2 * .5, v1 * .5, v2 / (v2 + v1 * f));
}

NUMERICS_EXPORT double fratiov(double p, long v1, long v2)
{
    double f0, f1, f2;
    long i = 0;
    
	if(!isPositive(v1)){
		NUMERICS_ERROR("fratiov", "Invalid v1 parameter");
	} else if(!isPositive(v2)){
		NUMERICS_ERROR("fratiov", "Invalid v2 parameter");
	} else if(!isZeroOne(p)){
		NUMERICS_ERROR("fratiov", "Invalid p parameter");
	} else {
		f0 = f1 = 5.0;
		while(fratiop(f0,v1,v2) > p) f0 -= .5;
		while(fratiop(f1,v1,v2) < p) f1 += .5;
		if(f0 < 0.0) f0 = 0.0;
		do {
			f2 = f0 + (f1-f0)/2.0;
			if(fratiop(f2,v1,v2) > p) f1 = f2;
			else f0 = f2;
	    } while(!areEqual(f0, f1, NUMERICS_MAX_ERROR) && ++i < 100);
	}
    return f2;
}

NUMERICS_EXPORT double normalp(double x, double m, double v)
{
	if(!isPositive(v)){
		NUMERICS_ERROR("normalp", "Invalid v parameter");
	} else if(areEqual(x, m, NUMERICS_MAX_ERROR)){
		return .5;
	}
    return .5*(erff((x-m)/(sqrt(2.0 * v))) + 1.0);
}

NUMERICS_EXPORT double normalv(double p, double m, double v)
{
    double q, r, z;
    
    q = p - .5;
    if(fabs(q) <= .42){
        r = q * q;
        z = q * (((-25.4410604963 * r + 41.39119773534) * r - 18.61500062529) * r + 2.50662823884)
            / ((((3.13082909833 * r -21.06224101826) * r + 23.08336743743) * r - 8.47351093090) * r + 1.0);
        return z * sqrt(v) + m;
    }
    r = p;
    if(q > 0.0) r = 1.0 - p;
    if(r > 0.0){
        r = sqrt(-log(r));
        z = (((2.32121276858 * r + 4.85014127135) * r - 2.29796479134) * r - 2.78718931138)
            / ((1.63706781897 * r + 3.54388924762) * r + 1.0);
        if(q < 0.0) return -z * sqrt(v) + m;
        return z * sqrt(v) + m;
    }
    return 0.0;
}

NUMERICS_EXPORT double studtp(double t, double v)
{
	if(!isPositive(v)){
		NUMERICS_ERROR("studtp", "Invalid v parameter");
	} else if(areEqual(t, 0.0, NUMERICS_MAX_ERROR)){
		return .5;
	} else if(!isPositive(t)){
		return 0.5 * betai(.5 * v, .5, v / (v + t * t));
	}
    return 1.0 - 0.5 * betai(.5 * v, .5, v / (v + t * t));
}

NUMERICS_EXPORT double studtv(double p, double v)
{
    double t0, t1, t2, xx, sum, x;
    int i;
    
    x = fabs(normalv(p, 0.0, 1.0));
    xx = x * x;
    sum = x + ((xx + 1.0) * x) / (4.0 * v);
    sum += (((5.0 * xx + 16.0) * xx + 3.0) * x) / (96.0 * v * v);
    sum += ((((3.0 * xx + 19.0) * xx + 17.0) * xx - 15.0) * x) / (384.0 * v * v * v);
    sum += (((((79.0 * xx + 776.0) * xx + 1482.0) * xx - 1920.0) * xx - 945.0) * x) / (92160.0 * v * v * v * v);
    if (p >= .5) t1 = sum;
    else t1 = -sum;
    t0 = t1;
    while(studtp(t0,v) > p) t0 -= .5;
    while(studtp(t1,v) < p) t1 += .5;
    i = 0;
    do {
        t2 = t0 + (t1-t0)/2.0;
        if(studtp(t2,v) > p) t1 = t2;
        else t0 = t2;
    } while(!areEqual(t0, t1, NUMERICS_MAX_ERROR) && ++i < 100);
    return t2;
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
 Version 1.1 - 04/10/1999 - Use areEqual to test for floating-point
                            equality.
                            Added use of domain calls.
===================================================================*/
