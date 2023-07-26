//===================================================================
// rootfind.h
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

#ifndef _ROOTFIND_H_
#define _ROOTFIND_H_

#include "numerics.h"

typedef double (*Function)(double);

double bisection(Function func, double x0, double x1, double tol);
//-------------------------------------------------------------------
// Returns the root of func bracketed between x0 and x1 using the
// bisection method.  The returned root is within tol of the true
// value of the root.
//-------------------------------------------------------------------

double falsePosition(Function func, double x0, double x1, double tol);
//-------------------------------------------------------------------
// Returns the root of func using the false position method with
// initial approximations x0 and x1.  The returned root is within tol
// of the true value of the root.
//-------------------------------------------------------------------

double newton(Function func, Function der, double x0, double tol);
//-------------------------------------------------------------------
// Returns the root of func with first derivative der using Newton's
// method with initial approximation x0.  The returned root is within
// tol of the true value of the root.
//-------------------------------------------------------------------

double secant(Function func, double x0, double x1, double tol);
//-------------------------------------------------------------------
// Returns the root of func using the secant method with initial
// approximations x0 and x1.  The returned root is within tol of the
// true value of the root.
//-------------------------------------------------------------------

#endif

//===================================================================
// Revision History
//
// Version 1.0 - 05/09/1999 - New.
//===================================================================
