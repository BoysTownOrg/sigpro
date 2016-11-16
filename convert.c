// convert.c - convert sampling rate

#include <stdlib.h>
#include <math.h>
#include "sigpro.h"

static float *hh;	    // low-pass filter impulse response
static int fltwid = 16;	    // low-pass filter width
static int mxupsr = 512;    // maximum up-sample ratio
static int o1, o2, wflg;
static int n_up, n_fw, upsr;

static int
rationalize(double rr, double ep)
{
    int i1, i2;
    double irat = 1, dif, mindif;

    i1 = i2 = o1 = o2 = 1;
    mindif = fabs(rr - 1);
    while((i1 <= mxupsr) && (i2 <= mxupsr) && mindif) {
	irat = (double) i2 / (double) i1;
	dif = fabs(rr - irat);
	if (mindif > dif) {
	    mindif = dif;
	    o1 = i1;
	    o2 = i2;
	}
	if (irat < rr) {
	    i2++;
	} else if(irat > rr) {
	    i1++;
	} else {
	    break;
	}
    }
    if (mindif < ep)
	return (1);
    o1 = 0;
    o2 = mxupsr;
    return (0);
}

static void
init_lpf()
{
    double  a, p, w;
    int     nn, i, j, k;
    float   s, sum;
    static double pi = M_PI;

    upsr = (o1 > o2) ? o1 : o2;
    n_up = o2;
    n_fw = (upsr * fltwid) / n_up;
    nn = n_up * n_fw + 1;
    hh = calloc(nn, sizeof(float));
    hh[0] = 1;
    for (i = 1; i < nn; i++) {
	a = (float) i * pi / upsr;
	p = (float) i * pi / (nn - 1);
	w = 0.42 + 0.5 * cos(p) + 0.08 * cos(2 * p);	/* Blackman window  */
        hh[i] = (float) (w * sin(a) / a);		/* LPF impulse response */
    }
    for (k = 0; k < n_up; k++) {			/* normalize */
	sum = 0;
	for (i = 0; i < n_fw; i++) {
	    j = i + 1;
	    sum += hh[i * n_up + k];
	    sum += hh[j * n_up - k];
	}
	s = 1 / sum;
	for (i = 0; i < n_fw; i++) {
	    j = i + 1;
	    hh[i * n_up + k] *= s;
	    hh[j * n_up - k] *= s;
	}
    }
}

static float
up(float *x1, int n1, int ii)
{
    float   x2_i, sum, x1_l, x1_r;
    int     jj, kk, i, j;

    if (wflg)
	ii %= o2 * n1;
    jj = ii / o2;
    kk = ii % o2;
    if (kk == 0) {
	x2_i = x1[jj];
    } else {
	if (hh == NULL)
	    init_lpf();
	sum = 0;
	for (i = 0; i < n_fw; i++) {
	    j = i + 1;
	    if (wflg) {
		x1_l = x1[(jj - i + n1) % n1];
		x1_r = x1[(jj + i + 1) % n1];
	    } else {
		if (((jj - i) < 0) || ((jj - i) >= n1))
		    x1_l = 0;
		else
		    x1_l = x1[jj - i];
		if (((jj + i +  1) < 0) || ((jj + i +  1) >= n1))
		    x1_r = 0;
		else
		    x1_r = x1[jj + i +  1];
	    }
	    sum += hh[i * n_up + kk] * x1_l;
	    sum += hh[j * n_up - kk] * x1_r;
	}
	x2_i = sum;
    }

    return (x2_i);
}

int
sp_convert(float *x1, int n1, float *x2, int n2,
	double rr, int wrap)
{
    double  k, ep;
    float   a;
    int     i, j;

    if ((n1 < 1) || (n2 < 1) || (rr < 0))
	return (1);
    if (rr == 0)
	rr = (double) n2 / (double) n1;
    ep = 1 / ((double) n1 * (double) mxupsr);
    wflg = wrap;
    hh = NULL;
    if (rr == 1) {				// copy
	if (wflg) {
	    for (i = 0; i < n2; i++) {
		x2[i] = x1[i % n1];
	    }
	} else {
	    for (i = 0; i < n2; i++) {
		x2[i] = (i < n1) ? x1[i] : 0;
	    }
	}
    } else if (rationalize(rr, ep)) {		// updn ratio = integer
	for (i = 0; i < n2; i++) {
	    x2[i] = up(x1, n1, i * o1);
	}
    } else {					// updn ratio = max
	for (i = 0; i < n2; i++) {
	    k = i * (o2 / rr);
	    j = (int) k;
	    a = (float) (k - j);
	    x2[i] = (1 - a) * up(x1, n1, j) + a * up(x1, n1, j + 1);
	}
    }
    if (hh != NULL)
	free(hh);
    return (0);
}
