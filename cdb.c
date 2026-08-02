// cdb.c

#include <math.h>
#include "sigpro.h"

#ifndef WIN32
#define _hypot  hypot
#endif

FUNC(void) sp_cdb(	// complex decibels (20*log10)
    float *x, 		// complex-input array
    float *y, 		// real-output array (dB)
    int n		// complex array size
)
{
    double  a;
    int     i, j, k;

    for (i = 0; i < n; i++) {
	j = i * 2;
	k = j + 1;
	a = _hypot(x[j], x[k]);
	if (a < 1e-40) {
	    y[i] = -800;
	} else {
	    y[i] = (float) (20 * log10(a));
	}
    }
}

FUNC(void)
sp_cgd(float *x, float *y, int n, double df)
{
    double  a, dp, ph, pr;
    int     i, j, k;
    static double tpi = 2 * M_PI;

    pr = 0;
    for (i = 0; i < n; i++) {
	j = i * 2;
	k = j + 1;
	a = _hypot(x[j], x[k]);
	if (a < 1e-40) {
	    ph = 0;
	} else {
	    ph = (float) (atan2(x[k], x[j]) / tpi);
	}
	dp = pr - ph;
	while (dp < -0.5) {
	    dp += 1;
	}
	while (dp > 0.5) {
	    dp -= 1;
	}
	y[i] = (float) (dp / df);
	pr = ph;
    }
}

FUNC(void)
sp_cph(float *x, float *y, int n)
{
    double  a;
    int     i, j, k;
    static double tpi = 2 * M_PI;

    for (i = 0; i < n; i++) {
	j = i * 2;
	k = j + 1;
	a = _hypot(x[j], x[k]);
	if (a < 1e-40) {
	    y[i] = 0;
	} else {
	    y[i] = (float) (atan2(x[k], x[j]) / tpi);
	}
    }
}

FUNC(void)
sp_unwrap(float *x, float *y, int n)
{
    int     i;

    y[0] = x[0];
    for (i = 1; i < n; i++) {
        y[i] = x[i] + (float) floor(y[i - 1] - x[i] + 0.5);
    }
}
