// filter.c

#include <stdlib.h>
#include <math.h>
#include "sigpro.h"

FUNC(int) sp_filterz(
    float *b, int nb, 
    float *a, int na,
    float *x, float *y, int n, 
    float *z
) {
    float   yyyy;
    int     i, j, k, nz;

    if (na < 0)
	return (1);
    if (nb <= 0)
	return (2);
    if (n <= 0)
	return (3);
    if (a == NULL)
	na = 0;
    // normalize coefficients
    if ((na > 0) && (a[0] != 1)) {
	for (i = 1; i < na; i++)
	    a[i] /= a[0];
	for (i = 0; i < nb; i++)
	    b[i] /= a[0];
	a[0] = 1;
    }
    nz = ((na > nb) ? na : nb) - 1;
    for (i = 0; i < n; i++) {
	// delay output
	yyyy = b[0] * x[i] + z[0];
	for (j = 0; j < nz; j++) {
	    k = j + 1;
	    if (k < nz)
		z[j] = z[k];
	    else
		z[j] = 0;
	    if (k < nb)
		z[j] += b[k] * x[i];
	    if (k < na)
		z[j] -= a[k] * yyyy;
	}
	y[i] = yyyy;
    }

    return (0);
}

FUNC(int) sp_filter(
    float *b, int nb, 
    float *a, int na,
    float *x, float *y, int n
) {
    float  *z;
    int     nz, err;

    nz = ((na > nb) ? na : nb) - 1;
    z = (float *) calloc(nz, sizeof(float));
    err = sp_filterz(b, nb, a, na, x, y, n, z);
    free(z);

    return (err);
}

FUNC(int) sp_transfer(
    float *stm, 
    float *rsp, 
    int np,
    float *H
) {
    float *x, *y, *z, xms;
    int i, ir, ii, nf, err = 0;
 
    nf = np /2 + 1; // number of frequencies
    x = (float *) calloc(np + 2, sizeof(float)); 
    y = (float *) calloc(np + 2, sizeof(float)); 
    z = H;
    for (i = 0; i < np; i++) {
        x[i] = stm[i];
        y[i] = rsp[i];
    }
    sp_rcfft(x, np);  // real, in-place FFT
    sp_rcfft(y, np);
    // calculate complex z = x / y;
    for (i = 0; i < nf; i++) {
        ir = 2 * i;     // index real
        ii = ir + 1;    // index imaginary
        xms = x[ir] * x[ir] + x[ii] * x[ii];  // magnitude squared 
        z[ir] = (x[ir] * y[ir] + x[ii] * y[ii]) / xms; 
        z[ii] = (x[ir] * y[ii] - x[ii] * y[ir]) / xms;
    }
    free(x);
    free(y);

    return (err);
}

FUNC(int) sp_freqz(
    float *b,
    int nb,
    float *a,
    int na,
    float *f, 
    float *H, 
    int nf,
    double fs
) {
    double br, bi, ar, ai, am, cr, ci, dp, ph;
    int i, j, ir, ii, nc, err = 0;

    nc = (nb > na) ? nb : na;
    for (i = 0; i < nf; i++) {
        dp = 2 * M_PI * f[i] / fs;
	br = bi = 0;
	ar = ai = 0;
	for (j = 0; j < nc; j++) {
	    ph = -j * dp;
	    cr = cos(ph);
	    ci = sin(ph);
	    if (j < nb) {
	        br += b[j] * cr;
	        bi += b[j] * ci;
	    }
	    if (j < na) {
	        ar += a[j] * cr;
	        ai += a[j] * ci;
	    }
	}
        if (na < 1) {
            ar = 1;
        }
	am = ar * ar + ai * ai;
        if (am == 0) {
            err = 1;
            break;
        }
	ir = i * 2;
	ii = ir + 1;
	H[ir] = (float) ((br * ar + bi * ai) / am);
	H[ii] = (float) ((bi * ar - br * ai) / am);
    }

    return (err);
}
