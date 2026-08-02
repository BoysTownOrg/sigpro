// vector.c

#include "sigpro.h"

/************************ real-vector functions *************************************/

// vector add

FUNC(void) sp_vadd(float *x, float *y, float *z, int n)
{
    int     i;

    for (i = 0; i < n; i++) {
	z[i] = x[i] + y[i];
    }
}

// vector dot-product

FUNC(double) sp_vdot(float *x, float *y, int n)
{
    double  d;
    int     i;

    d = 0;
    for (i = 0; i < n; i++) {
	d += x[i] * y[i];
    }

    return (d);
}

// vector divide

FUNC(int) sp_vdiv(float *x, float *y, float *z, int n)
{
    int     i, err = 0;

    for (i = 0; i < n; i++) {
        if (y[i] == 0) {
            err++;
        } else {
            z[i] = x[i] / y[i];
        }
    }

    return (err);
}

// vector maximum

FUNC(int) sp_vmax(float *x, int n)
{
    int i, j = 0;

    for (i = 1; i < n; i++) {
	if (x[j] < x[i]) {
	    j = i;
	}
    }

    return (j);
}

// vector minimum

FUNC(int) sp_vmin(float *x, int n)
{
    int i, j = 0;

    for (i = 1; i < n; i++) {
	if (x[j] > x[i]) {
	    j = i;
	}
    }

    return (j);
}

// vector multiply

FUNC(void) sp_vmul(float *x, float *y, float *z, int n)
{
    int     i;

    for (i = 0; i < n; i++) {
	z[i] = x[i] * y[i];
    }
}

// vector subtract

FUNC(void) sp_vsub(float *x, float *y, float *z, int n)
{
    int     i;

    for (i = 0; i < n; i++) {
	z[i] = x[i] - y[i];
    }
}

/************************ complex-vector functions *************************************/

// complex-vector add

FUNC(void) sp_cvadd(float *x, float *y, float *z, int n) 
{
    int i;

    for (i = 0; i < (n * 2); i++) {
        z[i] = x[i] + y[i];
    }
}

// complex-vector divide

FUNC(int) sp_cvdiv(float *x, float *y, float *z, int n) 
{
    double xr, xi, yr, yi, ym;
    int i, ir, ii, err = 0;

    for (i = 0; i < n; i++) {
	ir = i * 2;
	ii = ir + 1;
	xr = x[ir];
	xi = x[ii];
	yr = y[ir];
	yi = y[ii];
        ym = yr * yr + yi * yi;
        if (ym == 0) {
            err++;
        } else {
            z[ir] = (float) ((xr * yr + xi * yi) / ym);
            z[ii] = (float) ((xi * yr - xr * yi) / ym);
        }
    }

    return (err);
}

// complex-vector multiply

FUNC(void) sp_cvmul(float *x, float *y, float *z, int n) 
{
    double xr, xi, yr, yi;
    int i, ir, ii;

    for (i = 0; i < n; i++) {
	ir = i * 2;
	ii = ir + 1;
	xr = x[ir];
	xi = x[ii];
	yr = y[ir];
	yi = y[ii];
        z[ir] = (float) (xr * yr - xi * yi);
        z[ii] = (float) (xi * yr + xr * yi);
    }
}

// complex magnitude squared

FUNC(void) sp_cmagsq(float *x, float *y, int n)
{
    int i, ir, ii;

    for (i = 0; i < n; i++) {
	ir = 2 * i;
	ii = ir + 1;
	y[ir] = x[ir] * x[ir] + x[ii] * x[ii];
	y[ii] = 0;
    }
}

// complex-vector subtract

FUNC(void) sp_cvsub(float *x, float *y, float *z, int n) 
{
    int i;

    for (i = 0; i < (n * 2); i++) {
        z[i] = x[i] - y[i];
    }
}
