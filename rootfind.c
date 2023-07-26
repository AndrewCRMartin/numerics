//===================================================================
// rootfind.c
//
// Version 1.0
//
// Written by:
//   Brent Worden
//   WordenWare
//   email:  Brent@Worden.org
//
// Copyright (c) 1999 WordenWare
//
// Created:  May 09, 1999
// Revised:  
//===================================================================

#include "algorthm.h"
#include "numerror.h"
#include "rootfind.h"
#include "utility.h"

double bisection(Function func, double x0, double x1, double tol)
{
	double f0 = func(x0);
	double f1 = func(x1);

	if(f0 * f1 > 0.0){
		NUMERICS_ERROR("bisection", "func(x0) and func(x1) must have oppisite signs");
	} else {
		double fmid;
		double xmid;

		while(!areEqual(x0, x1, tol)){
			xmid = average(x0, x1);
			fmid = func(xmid);
			if(fmid * f0 > 0.0){
				f0 = fmid;
				x0 = xmid;
			} else {
				f1 = fmid;
				x1 = xmid;
			}
		}
	}
	return average(x0, x1);
}

double falsePosition(Function func, double x0, double x1, double tol)
{
	double f0 = func(x0);
	double f1 = func(x1);

	if(f0 * f1 > 0.0){
		NUMERICS_ERROR("falsePosition", "func(x0) and func(x1) must have oppisite signs");
	} else {
		double f;
		double x;
		int n = 0;

		do {
			x = x1 - f1 * (x1 - x0) / (f1 - f0);
			f = func(x);
			if(f * f1 < 0.0){
				x0 = x1;
				f0 = f1;
			}
			x1 = x;
			f1 = f;
			++n;
		} while(n < NUMERICS_ITMAX && !areEqual(x0, x1, tol));
		if(n >= NUMERICS_ITMAX){
			NUMERICS_ERROR("falsePosition", "Failed to converge to root.");
		}
	}
	return x1;
}

double newton(Function func, Function der, double x0, double tol)
{
	int n = 0;
	double x1 = x0;

	do {
		x0 = x1;
		x1 = x0 - func(x0) / der(x0);
		++n;
	} while(n < NUMERICS_ITMAX && !areEqual(x0, x1, tol));
	if(n >= NUMERICS_ITMAX){
		NUMERICS_ERROR("newton", "Failed to converge to root.");
	}

	return x1;
}

double secant(Function func, double x0, double x1, double tol)
{
	double f0 = func(x0);
	double f1 = func(x1);
	int n = 0;
	double x;

	do {
		x = x1 - f1 * (x1 - x0) / (f1 - f0);
		x0 = x1;
		f0 = f1;
		x1 = x;
		f1 = func(x);
		++n;
	} while(n < NUMERICS_ITMAX && !areEqual(x0, x1, tol));
	if(n >= NUMERICS_ITMAX){
		NUMERICS_ERROR("secant", "Failed to converge to root.");
	}

	return x1;
}

//===================================================================
// Revision History
//
// Version 1.0 - 05/09/1999 - New.
//===================================================================
