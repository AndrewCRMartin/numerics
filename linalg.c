/*===================================================================
 linalg.c

 Version 1.0

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  
===================================================================*/

#include <math.h>

#include "linalg.h"
#include "matutils.h"
#include "numerror.h"
#include "vecutils.h"

NUMERICS_EXPORT BOOL choldcmp(double **a, int n, double *p)
{
    int i, j, k;
    double sum;
    
    *p = sqrt(a[0][0]);
    for(i = 1; i < n; i++){
        a[i][0] /= *p;
    }
    for(i = 1; i < n; i++){
        sum = a[i][i];
        for(k = 0; k < i; ++k){
            sum -= a[i][k] * a[i][k];
        }
        if(sum <= 0.0){
            return FALSE;
        }
        *(p + i) = sqrt(sum);
        for(j = n - 1; j > i; --j){
            sum = a[j][i];
            for(k = 0; k < i; ++k){
                sum -= a[j][k] * a[i][k];
            }
            a[j][i] = sum / *(p + i);
        }
    }
    
    return TRUE;
}

NUMERICS_EXPORT void cholsl(double **a, int n, double *p, double *b, double *x)
{
    int i, k;
    double sum;
    
    for(i = 0; i < n; i++){
        for(sum = *(b + i), k = i - 1; k >= 0; --k){
            sum -= a[i][k] * *(x + k);
        }
        *(x + i) = sum / *(p + i);
    }
    for(i = n - 1; i >= 0; --i){
        for(sum = *(x + i), k = i + 1; k < n; ++k){
            sum -= a[k][i] * *(x + k);
        }
        *(x + i) = sum / *(p + i);
    }
}

NUMERICS_EXPORT double determinant(double **a, int n)
{
    double ret;
    int i;
    int *indx = ivector(0, n - 1);
    
    if(ludcmp(a, n, indx, &ret)){
        for(i = 0; i < n; i++){
            ret *= a[i][i];
        }
    } else {
        ret = 0.0;
    }
    
    free_ivector(indx, 0);
    
    return ret;
}

NUMERICS_EXPORT BOOL inverse(double **a, int n)
{
    double d;
    int i, j;
    BOOL ret = FALSE;
    double** ai = dmatrix(0, n - 1, 0, n - 1);
    double* col = dvector(0, n - 1);
    int* indx = ivector(0, n - 1);
    
    if(ludcmp(a, n, indx, &d)){
        for(j = 0; j < n; j++){
            for(i = 0; i < n; i++) col[i] = 0.0;
            col[j] = 1.0;
            lubksb(a, n, indx, col);
            for(i = 0; i < n; i++) ai[i][j] = col[i];
        }
        for(i = 0; i < n; i++){
            for(j = 0; j < n; j++){
                a[i][j] = ai[i][j];
            }
        }
        ret = TRUE;
    }
    
    free_dmatrix(ai, 0, n - 1, 0);
    free_dvector(col, 0);
    free_ivector(indx, 0);
    
    return ret;
}

NUMERICS_EXPORT BOOL linsolve(double **m, double *b, int n, int method)
{
    int* indx;
    double d;
    BOOL ret = FALSE;
    
    if(method | LINSOLVE_LU){
        indx = ivector(0, n - 1);
        ret = ludcmp(m, n, indx, &d);
        if(ret){
            lubksb(m, n, indx, b);
        }
        free_ivector(indx, 0);
    }
    
    return ret;
}

NUMERICS_EXPORT void lubksb(double **m, int n, int *indx, double *b)
{
    int i, ii = -1, ip, j;
    double sum;
    
    for(i = 0; i < n; i++){
        ip = *(indx + i);
        sum = *(b + ip);
        *(b + ip) = *(b + i);
        if(ii > -1){
            for(j = ii; j <= i - 1; j++){
                sum -= m[i][j] * *(b + j);
            }
        } else if(sum){
            ii = i;
        }
        *(b + i) = sum;
    }
    for(i = n - 1; i >= 0; i--){
        sum = *(b + i);
        for(j = i + 1; j < n; j++){
            sum -= m[i][j] * *(b + j);
        }
        *(b + i) = sum / m[i][i];
    }
}

NUMERICS_EXPORT BOOL ludcmp(double **m, int n, int *indx, double *d)
{
    int i, imax, j, k;
    double big, dum, sum, temp;
    double* vv = dvector(0, n - 1);
    
    *d = 1.0;
    for(i = 0; i < n; i++){
        big = 0.0;
        for(j = 0; j < n; j++){
            if((temp = fabs(m[i][j])) > big){
                big = temp;
            }
        }
        if(big == 0.0){
            free_dvector(vv, 0);
            NUMERICS_ERROR("ludcmp", "Singular Matrix");
            return FALSE;
        }
        vv[i] = 1.0 / big;
    }
    for(j = 0; j < n; j++){
        for(i = 0; i < j; i++){
            sum = m[i][j];
            for(k = 0; k < i; k++){
                sum -= m[i][k] * m[k][j];
            }
            m[i][j] = sum;
        }
        big = 0.0;
        for(i = j; i < n; i++){
            sum = m[i][j];
            for(k = 0; k < j; k++){
                sum -= m[i][k] * m[k][j];
            }
            m[i][j] = sum;
            if((dum = vv[i] * fabs(sum)) >= big){
                big = dum;
                imax = i;
            }
        }
        if(j != imax){
            for(k = 0; k < n; k++){
                dum = m[imax][k];
                m[imax][k] = m[j][k];
                m[j][k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        *(indx + j) = imax;
        if(m[j][j] == 0.0){
            m[j][j] = NUMERICS_FLOAT_MIN;
        }
        if(j != n - 1){
            dum = 1.0 / (m[j][j]);
            for(i = j + 1; i < n; i++){
                m[i][j] *= dum;
            }
        }
    }
    
    free_dvector(vv, 0);
    
    return TRUE;
};

NUMERICS_EXPORT void matmat(double **a, int nra, int nca, double **b, int ncb, double **prod)
{
    int i, j, k;
    double sum;
    
    for(i = 0; i < nra; i++){
        for(j = 0; j < ncb; j++){
            sum = 0.0;
            for(k = 0; k < nca; k++){
                sum += a[i][k] * b[k][j];
            }
            prod[i][j] = sum;
        }
    }
}

NUMERICS_EXPORT void matvec(double **a, int nra, int nca, double *x, double *b)
{
    int i, j;
    double sum;
    
    for(i = 0; i < nra; i++){
        sum = 0.0;
        for(j = 0; j < nca; j++){
            sum += a[i][j] * x[j];
        }
        b[i] = sum;
    }
}

NUMERICS_EXPORT void transpose(double **a, int nr, int nc, double **at)
{
    int i, j;
    
    for(i = 0; i < nr; ++i){
        for(j = 0; j < nc; ++j){
            at[j][i] = a[i][j];
        }
    }
}

NUMERICS_EXPORT void vecmat(double *x, double **a, int nra, int nca, double *b)
{
    double** t = dmatrix(0, nca - 1, 0, nra - 1);
    
    transpose(a, nra, nca, t);
    matvec(t, nca, nra, x, b);
    
    free_dmatrix(t, 0, nca - 1, 0);
}

NUMERICS_EXPORT double vecvec(double *first1, double* last1, double* first2)
{
    double p = 1.0;
    
    while(first1 < last1){
        p += *first1 * *first2;
        ++first1;
        ++first2;
    }
    
    return p;
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
===================================================================*/

