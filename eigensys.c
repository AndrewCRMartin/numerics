/*====================================================================
 eigensys.c

 Version 1.1

 Written by:
   Brent Worden
   WordenWare
   email:  brent.worden@poboxes.com

 Copyright (c) 1998-1999 WordenWare

 Created:  August 28, 1998
 Revised:  May 11, 1999
====================================================================*/

#include <math.h>

#include "descript.h"
#include "eigensys.h"
#include "mathx.h"
#include "numerror.h"

NUMERICS_EXPORT void tridred(double** a, int n, double* d, double* e)
{
    double f, g, h;
    int i, j, k;
    double scale, tmp;

	for(i = 0; i < n; ++i){
		d[i] = a[i][n - 1];
		a[i][n - 1] = a[i][i];
	}
	for(i = n - 1; i >= 0; --i){
		h = 0.0;
		scale = 0.0;
		if(i < 1){
			e[i] = 0.0;
		} else {
			for(k = 0; k <= i; ++k){
				scale += fabs(d[k]);
			}
			if(scale == 0.0){
				for(j = 0; j < i; ++j){
					d[j] = a[j][i - 1];
					a[j][i - 1] = a[j][i];
					a[j][i] = 0.0;
				}
				e[i] = 0.0;
			} else {
				for(k = 0; k < i; ++k){
					tmp = (d[k] /= scale);
					h += tmp * tmp;
				}
				f = d[i - 1];
				g = (f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i] = scale * g;
				h -= f * g;
				d[i - 1] = f - g;
				if(i != 1){
					for(j = 0; j < i; ++j){
						e[j] = 0.0;
					}
					for(j = 0; j < i; ++j){
						f = d[j];
						g = e[j] + a[j][j] * f;
						k = j + 1;
						if(i > k){
							for(; k < i; ++k){
								g += a[j][k] * d[k];
								e[k] += a[j][k] * f;
							}
						}
						e[j] = g;
					}
					f = 0.0;
					for(j = 0; j < i; ++j){
						e[j] /= h;
						f += e[j] * d[j];
					}
					h = f / (h + h);
					for(j = 0; j < i; ++j){
						e[j] -= h * d[j];
					}
					for(j = 0; j < i; ++j){
						f = d[j];
						g = e[j];
						for(k = j; k < i; ++k){
							a[j][k] -= (f * e[k] + g * d[k]);
						}
					}
				}	
				for(j = 0; j < i; ++j){
					f = d[j];
					d[j] = a[j][i - 1];
					a[j][i - 1] = a[j][i];
					a[j][i] = f * scale;
				}
			}
		}
	}
}

NUMERICS_EXPORT void trideig(double *dfirst, double *dlast, double *e)
{
 	int n = length(dfirst, dlast);
	double b, c, f, g;
	int i, j, k, m;
	double p, r, s;
	double tmp;
	double *d = dfirst;

	if(n != 1) {
		for(i = 1; i < n; ++i){
			e[i - 1] = e[i];
		}
		e[n - 1] = 0.0;
		for(k = 0; k < n; ++k){
			j = 0;
			for(m = k + 1; m < n; ++m){
				tmp = fabs(d[m - 1]) + fabs(d[m]);
				if(tmp == fabs(e[m - 1]) + tmp){
					break;
				}
			}
			p = d[k];
			if(m != k + 1){
				if(j++ == NUMERICS_ITMAX){
					NUMERICS_ERROR("tridql", "Routine failed to converge.");
					return;
				}
				g = (d[k + 1] - p) / (e[k] * 2.0);
				r = pythag(g, 1.0);
				g = d[m - 1] - p + e[k] / (g + (g >= 0.0 ? r : -r));
				s = 1.0;
				c = 1.0;
				p = 0.0;
				for(i = m - 2; i >= 0; --i){
					tmp = e[i];
					f = s * tmp;
					b = c * tmp;
					r = pythag(f, g);
					e[i + 1] = r;
					if(r == 0.0){
						d[i + 1] -= p;
						e[m - 1] = 0.0;
						break;
					}
					s = f / r;
					c = g / r;
					g = d[i + 1] - p;
					r = (d[i] - g) * s + c * 2.0 * b;
					p = s * r;
					d[i + 1] = g + p;
					g = c * r - b;
				}
				if(i <= k + 1 || r != 0.0){
					d[k] -= p;
					e[k] = g;
					e[m - 1] = 0.0;
				}
			}
		}
	}
}

/*===================================================================
 Revision History

 Version 1.0 - 08/28/1998 - New.
 Version 1.1 - 05/11/1999 - Removed calcution of eigenvectors from
                            tridred and trideig
===================================================================*/
