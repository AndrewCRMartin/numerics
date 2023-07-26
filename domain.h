//===================================================================
// domain.h
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

#ifndef _DOMAIN_H_
#define _DOMAIN_H_

#include "numerics.h"

#define isPositive(x) ((x) > 0)
//-------------------------------------------------------------------
// Evaluates to true if x is positive, false otherwise.
//-------------------------------------------------------------------

#define isNegative(x) ((x) < 0)
//-------------------------------------------------------------------
// Evaluates to true if x is positive, false otherwise.
//-------------------------------------------------------------------

#define isNonNegative(x) ((x) >= 0)
//-------------------------------------------------------------------
// Evaluates to true if x is positive or zero, false otherwise.
//-------------------------------------------------------------------

#define isNonPositive(x) ((x) <= 0)
//-------------------------------------------------------------------
// Evaluates to true if x is negative or zero, false otherwise.
//-------------------------------------------------------------------

#define isZeroOne(x) (0 <= (x) && (x) <= 1)
//-------------------------------------------------------------------
// Evaluates to true if x is an element of [0.0, 1.0], false
// otherwise.
//-------------------------------------------------------------------

#endif

//===================================================================
// Revision History
//
// Version 1.0 - 05/09/1999 - New.
//===================================================================
