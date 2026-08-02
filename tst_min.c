//tst_min.c - test fmins function

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sigpro.h"

#define NV 4

static float c[NV] = {11, 12, 13, 14};

static double 
test_variance(float *x)
{
    int     i;
    double  y, z = 0;

    for (i = 0; i < NV; i++) {
	y = (x[i] - c[i]);
	z += y * y;
    }
    return (sqrt(z));
}

static void 
report(float *x)
{
    printf(" x1=%4.1f x2=%4.1f x3=%4.1f x4=%4.1f\n",
        x[0], x[1], x[2], x[3]);
}

int 
main(int ac, char **av)
{
    static float x[4] = {3, 2, 1, 1};

    printf("  initial:");
    report(x);
    sp_fmins(x, 4, test_variance, NULL);
    printf("    final:");
    report(x);
    printf("should be:");
    report(c);

    return (0);
}

