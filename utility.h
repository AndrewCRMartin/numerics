//===================================================================
// utility.h
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

#ifndef _UTILITY_H_
#define _UTILITY_H_

#include "numerics.h"

#define areEqual(a, b, tol) ((b) - (tol) <= (a)) && ((a) <= (b) + (tol))
//--------------------------------------------------------------------
// Evaluates to true if a and b are within tol of each other.
// Otherwise, return false.
//--------------------------------------------------------------------

#endif

//===================================================================
// Revision History
//
// Version 1.0 - 05/09/1999 - New.
//===================================================================
