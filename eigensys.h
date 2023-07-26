/*====================================================================
 eigensys.hpp

 Version 1.1

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  May 11, 1999
====================================================================*/

#ifndef _EIGENSYS_H_
#define _EIGENSYS_H_

#include "numerics.h"

NUMERICS_EXPORT void tridred(double** a, int n, double *d, double *e);
/*--------------------------------------------------------------------
 Householder reduction of a real, symmetric matrix a[0..n-1][0..n-1].
 On output, a is replaced by the orthogonal matrix effecting the
 transformation.  d[0..n-1] returns the diagonal elements of the
 tridiagonal matrix, and e[0..n-1] the off-diagonal elements, with
 e[0] = 0.
--------------------------------------------------------------------*/

NUMERICS_EXPORT void trideig(double *dfirst, double *dlast, double *e);
/*--------------------------------------------------------------------
 Algorithm to determine the eigenvalues of a real, symetric,
 tridiagonal matrix, or of a real, symmetric matrix previously
 reduced by tridred().  On input, d[0..n-1] contains the diagonal
 elements of the tridaigonal matrix.  On output, it returns the
 eigenvalues.  e[0..n-1] inputs the subdiagonal elements of the
 tridiagonal matrix, with e[0] arbitrary.  On output, e is destroyed.
--------------------------------------------------------------------*/

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
 Version 1.1 - 05/11/1999 - Removed calcution of eigenvectors from
                            tridred and trideig
===================================================================*/

