// scalar.c

#include "sigpro.h"

FUNC(void) sp_sadd(float *x, float *y, int n, double a)
{
    int     i;

    for (i = 0; i < n; i++) {
	y[i] = x[i] + (float) a;
    }
}

FUNC(void) sp_smul(float *x, float *y, int n, double b)
{
    int     i;

    for (i = 0; i < n; i++) {
	y[i] = x[i] * (float) b;
    }
}


FUNC(void) sp_sma(float *x, float *y, int n, double b, double a)
{
    int     i;

    for (i = 0; i < n; i++) {
	y[i] = ((float) b) * x[i] + ((float) a);
    }
}

FUNC(void) sp_copy(float *x, float *y, int n)
{
    int     i;

    for (i = 0; i < n; i++) {
	y[i] = x[i];
    }
}

FUNC(void) sp_zero(float *y, int n)
{
    int     i;

    for (i = 0; i < n; i++) {
	y[i] = 0;
    }
}

