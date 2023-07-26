/*===================================================================
 linalg.h

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#ifndef _LINALG_H_
#define _LINALG_H_

#include "numerror.h"

NUMERICS_EXPORT BOOL choldcmp(double **a, int n, double *p);
/*-------------------------------------------------------------------
 Given a positive-definit symmetric matrix a[0..n-1][0..n-1], this
 routine constructs its Colesky decomposition, A = LL'.  On input,
 only the upper triangle of a need be given; it is not modified.
 the Cholesky factor L is returned in the lower traingle of a,
 except for its diagonal elements which are returned in p[0..n-1].
 The routine returns true if the inversion was successful, otherwise
 it returns FALSE.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void cholsl(double **a, int n, double *p, double *b, double *x);
/*-------------------------------------------------------------------
 Solves the set of n linear equations Ax = b, where a is a positive-
 definite symmetric matrix.  a[0..n-1][0..n-1] and p[0..n-1] are
 input as the output of the routine choldcmp.  Only the lower
 triangle of a is accessed.  b[0..n-1] is input as the right-hand
 side vector.  The solution vector is returned in x[0..n-1].  a, n,
 and p are not modified and can be left in place for successive
 calls with different right-hand sides b.  b is not modified unless
 you identify b as x in the calling sequence, which is allowed.
-------------------------------------------------------------------*/

NUMERICS_EXPORT double determinant(double **a, int n);
/*-------------------------------------------------------------------
 Returns the determinant of the matrix a[0..n-1][0..n-1].  a is
 destroyed by this routine.
-------------------------------------------------------------------*/

NUMERICS_EXPORT BOOL inverse(double **a, int n);
/*-------------------------------------------------------------------
 Replaces the matrix a[0..n-1][0..n-1] with its inverse. The
 routine returns TRUE if the inversion was successful, otherwise it
 returns FALSE.
-------------------------------------------------------------------*/

enum {LINSOLVE_LU = 0x1};

NUMERICS_EXPORT BOOL linsolve(double **m, double *b, int n, int method);
/*-------------------------------------------------------------------
 Given a matrix mat[0..n-1][0..n-1] and vector b[0..n-1], the
 solution vector x is found for the linear system m.x = b.  The
 solution vector is returned in b.  The routine returns TRUE if the
 solution vector is successfully found, otherwise it returns FALSE.
 Both m and b are destroyed by this routine.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void lubksb(double **mat, int n, int *indx, double *b);
/*-------------------------------------------------------------------
 Given a matrix mat[0..n-1][0..n-1] and permutation vector
 indx[0..n-1] returned from ludcmp, this routines solves the set of
 n linear equations mat.x = b.  b[0..n-1] is input as the
 right-hand side vector and returns the solution vector x.  a, n,
 and indx are not modified by this routine and can be left in place
 for successive colls with different right-hand sides b.
-------------------------------------------------------------------*/

NUMERICS_EXPORT BOOL ludcmp(double **mat, int n, int *indx, double *d);
/*-------------------------------------------------------------------
 Given a matrix a[0..n-1][0..n-1], this routines replaces it by the
 LU decomposition of a rowwise permutation of itself.  a and n are
 input.  a is output, with the diagonal elements of the lower
 triangular matrix are equal to 1.  indx[0..n-1] is an output
 vector that records the row permutation effected by the partial
 pivoting.  d is output as 1 or -1 depending on whether the number
 of row interchanges was even or odd, respectively.  This routine
 is used in combination with lubksb to solve linear equations or
 invert a matrix.  The routine returns TRUE if the decomposition
 was successful, otherwise it returns FALSE.
-------------------------------------------------------------------*/

NUMERICS_EXPORT void matmat(double **a, int nra, int nca, double **b, int ncb, double **prod);
/*-------------------------------------------------------------------
 Postmultiplies the matrix a[0..nra-1][0..nca-1] by the matrix
 b[0..nca-1][0..ncb-1] and returns the product in the matrix
 prod[0..nra-1][0..ncb-1].
-------------------------------------------------------------------*/

NUMERICS_EXPORT void matvec(double **a, int nra, int nca, double *x, double *b);
/*-------------------------------------------------------------------
 Postmultiplies the matrix a[0..nra-1][0..nca-1] by the vector
 x[0..nca-1] and returns the product in the vector b[0..nra-1].
-------------------------------------------------------------------*/

NUMERICS_EXPORT void transpose(double **a, int nr, int nc, double **at);
/*-------------------------------------------------------------------
 Returns the transpose of a[0..nr-1][0..nc-1] as
 at[0..nc-1][0..nr-1].
-------------------------------------------------------------------*/

NUMERICS_EXPORT void vecmat(double *x, double **a, int nra, int nca, double *b);
/*-------------------------------------------------------------------
 Premultiplies the matrix a[0..nra-1][0..nca-1] by the vector
 x[0..nra-1] and returns the product in the vector b[0..nca-1].
-------------------------------------------------------------------*/

NUMERICS_EXPORT double vecvec(double *first1, double* last1, double* first2);
/*-------------------------------------------------------------------
 Returns the inner product between the vectors u[0..n-1] and
 v[0..n-1].
-------------------------------------------------------------------*/

#endif

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/
