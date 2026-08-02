// randflat

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sigpro.h"

FUNC(int) sp_randflat(float *x, int n)
{
    double  v;
    float  *p, *y;
    int     i, ir, ii, m;
    static double tpi = 2 * M_PI;

    m = n / 2;
    p = (float *) calloc(m, sizeof(float));
    y = (float *) calloc(n + 2, sizeof(float));
    sp_rand(p, m);
    y[0] = 1;
    for (i = 1; i < m; i++) {
	ir = i * 2;
	ii = ir + 1;
	v = tpi * p[i];
	y[ir] = (float) cos(v);
	y[ii] = (float) sin(v);
    }
    y[n] = 1;
    sp_crfft(y, n);
    memcpy(x, y, n * sizeof(float));
    free(p);
    free(y);
    return (0);
}
