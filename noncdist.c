/*===================================================================
 noncdist.cpp

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

#include "domain.h"
#include "mathx.h"
#include "noncdist.h"
#include "normdist.h"
#include "numerics.h"
#include "numerror.h"

#define R2DIVPI ( sqrt(2.0 / NUMERICS_PI) )
#define LOGRPI  ( log(sqrt(NUMERICS_PI)) )

NUMERICS_EXPORT double ncstudtp(double x, double df, double delta)
{
    double tnc = 0.0;
    double del = delta;
    BOOL negdel = FALSE;
    int en = 1;
    double xx;
    
	if(!isPositive(df)){
        NUMERICS_ERROR("ncstudtp", "Invalid df");
    } else {
	    if(x < 0.0){
	        negdel = TRUE;
		    del *= -1;
		}
    
		en = 1;
		xx = x * x / (x * x + df);
    
		if(xx > 0.0){
			double lambda = del*del;
			double p = .5*exp(-.5*lambda);
			double q = R2DIVPI * p * del;
			double s = .5 - p, a = .5, b = .5*df;
			double rxb = pow(1.0-xx, b);
			double albeta = LOGRPI + gammln(b) - gammln(a+b);
			double xodd = betai(a, b, xx);
			double godd = 2.0*rxb*exp(a*log(x)-albeta);
			double xeven = 1.0-rxb, geven = b*xx*rxb;
			double errbd;
			tnc = p * xodd + q * xeven;
			do {
				a += 1.0;
				xodd -= godd;
				xeven -= geven;
				godd *= (xx*(a+b-1.0)/a);
				geven *= (xx*(a+b-.5)/(a+.5));
				p *= (lambda/(2.0*en));
				q *= (lambda/(2.0*en+1.0));
				s -= p;
				++en;
				tnc += p*xodd + q*xeven;
				errbd = 2.0 * s * (xodd-godd);
			} while(errbd > NUMERICS_MAX_ERROR && en <= NUMERICS_ITMAX);
		}
    
		if(en > NUMERICS_ITMAX){
			NUMERICS_ERROR("ncstudtp", "Iteration failed to converge");
		} else {
			tnc += (1.0 - normalp(del, 0.0, 1.0));
			if(negdel){
				tnc = 1.0 - tnc;
			}
		}
	}    
    return tnc;
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
 Version 1.1 - 04/10/1999 - Added use of domain calls.
===================================================================*/
