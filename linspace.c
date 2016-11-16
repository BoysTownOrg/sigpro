// linspace.c

#include "sigpro.h"

FUNC(void) sp_linspace(float *x, int n, double a, double b)
{
    int     i;
    float   d;

    d = (float) ((b - a) / (n - 1));
    for (i = 0; i < n; i++) {
        x[i] = ((float) a) + i * d;
    }
}

