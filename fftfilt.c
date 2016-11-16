// fftfilt.c

#include <stdlib.h>
#include <math.h>
#include "sigpro.h"

/*****************************************************/

FUNC(int) sp_nxtpow2(
    int n
) {
    int     m;

    m = 1;
    while (m < n)
	m *= 2;
    return (m);
}

/*****************************************************/

FUNC(int) sp_fftfiltz(
    float *b, int nb, 
    float *x, float *y, int n, 
    float *z
) {
    float  *bb, *xx, xxr, xxi;
    int     nt, nn, nf, err, i, j, jr, ji, k;

    nt = sp_nxtpow2(nb) * 2;
    nn = nt - nb;
    nf = nt / 2 + 1;
    bb = (float *) calloc(nt + 2, sizeof(float));
    xx = (float *) calloc(nt + 2, sizeof(float));
    for (j = 0; j < nb; j++) {
	bb[j] = b[j];
    }
    err = sp_rcfft(bb, nt);
    for (i = 0; i < n; i += nn) {
	for (j = 0; j < nt; j++) {
	    xx[j] = ((j < nn) && ((i + j) < n)) ? x[i + j] : 0;
	}
        err = sp_rcfft(xx, nt);
	for (j = 0; j < nf; j++) {
	    jr = 2 * j;
	    ji = jr + 1;
	    xxr = xx[jr] * bb[jr] - xx[ji] * bb[ji];
	    xxi = xx[jr] * bb[ji] + xx[ji] * bb[jr];
	    xx[jr] = xxr;
	    xx[ji] = xxi;
	}
        err = sp_crfft(xx, nt);
	k = ((n - i) < nn) ? (n - i) : nn;
	for (j = 0; j < k; j++) {
	    y[i + j] = (j < nb) ? xx[j] + z[j] : xx[j];
	}
	for (j = 0; j < nb; j++) {
	    z[j] = xx[j + k];
	}
    }
    free(bb);
    free(xx);

    return (err);
}

FUNC(int) sp_fftfilt(
    float *b, int nb, 
    float *x, float *y, int n, 
    int wrap
) {
    float  *z;
    int     i, err;

    z = (float *) calloc(nb, sizeof(float));
    err = sp_fftfiltz(b, nb, x, y, n, z);
    if (!err && wrap) {
	for (i = 0; i < nb; i++) {
	    y[i % n] += z[i];
	}
    }
    free(z);

    return (err);
}
