/* crfft.c */

#include <stdlib.h>
#include <math.h>
#include "sigpro.h"

//-----------------------------------------------------------

static double tpi = 2 * M_PI;
static double sr2 = 1.414213562731;

static int
ilog2(int n)
{
    int     m;

    for (m = 1; m < 32; m++)
	if (n == (1 << m))
	    return (m);
    return (-1);
}

static int
bitrev(int ii, int m)
{
    register int jj;

    jj = ii & 1;
    --m;
    while (--m > 0) {
	ii >>= 1;
	jj <<= 1;
	jj |= ii & 1;
    }
    return (jj);
}

static void
rad2(int ii, float *x0, float *x1)
{
    int     k;
    float   t;

    for (k = 0; k < ii; k++) {
	t = x0[k] + x1[k];
	x1[k] = x0[k] - x1[k];
	x0[k] = t;
    }
}

static void
rad4(int jj, int nn,
    float *x0, float *x1, float *x2, float *x3,
    float *x4, float *x5, float *x6, float *x7)
{
    double  arg, tpiovn;
    float   c1, c2, c3, s1, s2, s3;
    float   t0, t1, t2, t3, t4, t5, t6, t7;
    int     ii, j, j0, ji, jr, jl, jlast, j4, k, k0, kl, m, n, ni;

    tpiovn = tpi / nn;
    ji = 3;
    jl = 2;
    jr = 2;
    n = nn / 4;
    for (m = 1; (1 << m) < n; m++)
	continue;
    ni = (n + 1) / 2;
    for (ii = 0; ii < ni; ii++) {
	if (ii == 0) {
	    for (k = 0; k < jj; k++) {
		t0 = x0[k] + x1[k];
		t1 = x0[k] - x1[k];
		t2 = x2[k] * 2;
		t3 = x3[k] * 2;
		x0[k] = t0 + t2;
		x2[k] = t0 - t2;
		x1[k] = t1 + t3;
		x3[k] = t1 - t3;
	    }
	    if (nn > 4) {
		k0 = jj * 4;
		kl = k0 + jj;
		for (k = k0; k < kl; k++) {
		    t2 = x0[k] - x2[k];
		    t3 = x1[k] + x3[k];
		    x0[k] = (x0[k] + x2[k]) * 2;
		    x2[k] = (x3[k] - x1[k]) * 2;
		    x1[k] = (float) ((t2 + t3) * sr2);
		    x3[k] = (float) ((t3 - t2) * sr2);
		}
	    }
	} else {
	    arg = tpiovn * bitrev(ii, m);
	    c1 = (float) cos(arg);
	    s1 = (float) -sin(arg);
	    c2 = c1 * c1 - s1 * s1;
	    s2 = c1 * s1 + c1 * s1;
	    c3 = c1 * c2 - s1 * s2;
	    s3 = c2 * s1 + s2 * c1;
	    j4 = jj * 4;
	    j0 = jr * j4;
	    k0 = ji * j4;
	    jlast = j0 + jj;
	    for (j = j0; j < jlast; j++) {
		k = k0 + j - j0;
		t0 = x0[j] + x6[k];
		t1 = x7[k] - x1[j];
		t2 = x0[j] - x6[k];
		t3 = x7[k] + x1[j];
		t4 = x2[j] + x4[k];
		t5 = x5[k] - x3[j];
		t6 = x5[k] + x3[j];
		t7 = x4[k] - x2[j];
		x0[j] = t0 + t4;
		x4[k] = t1 + t5;
		x1[j] = (t2 + t6) * c1 - (t3 + t7) * s1;
		x5[k] = (t2 + t6) * s1 + (t3 + t7) * c1;
		x2[j] = (t0 - t4) * c2 - (t1 - t5) * s2;
		x6[k] = (t0 - t4) * s2 + (t1 - t5) * c2;
		x3[j] = (t2 - t6) * c3 - (t3 - t7) * s3;
		x7[k] = (t2 - t6) * s3 + (t3 - t7) * c3;
	    }
	    jr += 2;
	    ji -= 2;
	    if (ji <= jl) {
		ji = 2 * jr - 1;
		jl = jr;
	    }
	}
    }
}

static void
reorder1(int m, float *x)
{
    int     j, k, kl, n;
    float   t;

    k = 4;
    kl = 2;
    n = 1 << m;
    for (j = 4; j <= n; j += 2) {
	if (k > j) {
	    t = x[j - 1];
	    x[j - 1] = x[k - 1];
	    x[k - 1] = t;
	}
	k -= 2;
	if (k <= kl) {
	    k = 2 * j;
	    kl = j;
	}
    }
}

static void
reorder2(int m, float *x)
{
    int     ji, ij, n;
    float   t;

    n = 1 << m;
    for (ij = 0; ij <= (n - 2); ij += 2) {
	ji = bitrev(ij >> 1, m) << 1;
	if (ij < ji) {
	    t = x[ij];
	    x[ij] = x[ji];
	    x[ji] = t;
	    t = x[ij + 1];
	    x[ij + 1] = x[ji + 1];
	    x[ji + 1] = t;
	}
    }
}

//-----------------------------------------------------------

static int
crfft2(float *x, int m)
{
    int     n, i, it, nn, jj, m2;

    n = 1 << m;
    x[1] = x[n];
    m2 = m / 2;

// re-order

    for (i = 3; i < n; i += 2)
	x[i] = -x[i];
    reorder2(m, x);
    reorder1(m, x);

// radix 4

    if (m2 != 0) {
	nn = 4 * n;
	for (it = 0; it < m2; it++) {
	    nn = nn / 4;
	    jj = n / nn;
	    rad4(jj, nn, x, x + jj, x + 2 * jj, x + 3 * jj,
		x, x + jj, x + 2 * jj, x + 3 * jj);
	}
    }

// radix 2

    if (m > m2 * 2) {
	jj = n / 2;
	rad2(jj, x, x + jj);
    }

    return (0);
}

//-----------------------------------------------------------

static int
crdft(float *x, int n)
{
    double a, ur, ui, vr, vi, wr, wi, xr, xi, yr, zz;
    float *x0, *x1, *y;
    int i, j, jj, k;
    static double tpi = 2 * M_PI;

    k = n / 2;
    y = (float *) calloc(n, sizeof(float));
    x0 = x;
    x1 = x0 + 1;
    a = tpi / n;
    ur = cos(a);
    ui = sin(a);
    for (i = 0; i < n; i++) {
	if (i == 0) {
	    vr = 1;
	    vi = 0;
	} else {
	    zz = ur * vr - ui * vi;
	    vi = ur * vi + ui * vr;
	    vr = zz;
	}
        xr = x0[0];
        xi = x1[0];
	yr = xr * vr - xi * vi;
	wr = vr;
	wi = vi;
	for (j = 1; j <= k; j++) {
	    zz = vr * wr - vi * wi;
	    wi = vr * wi + vi * wr;
	    wr = zz;
	    jj = j * 2;
	    xr = x0[jj];
	    xi = x1[jj];
	    yr += xr * wr - xi * wi;
	    if (jj < n)
		yr += xr * wr - xi * wi;
	}
	y[i] = (float) yr;
    }
    for (i = 0; i < n; i++) {
	x[i] = y[i];
    }
    free(y);

    return (0);
}

//-----------------------------------------------------------

// sp_crfft - complex-to-real inverse FFT

FUNC(int) sp_crfft(float *x, int n)
{
    int m, err, i;

    if (n < 2)
	return (0);

    m = ilog2(n);
    if (m > 0)
	err = crfft2(x, m);
    else
        err = crdft(x, n);

// scale inverse by 1/n

    for (i = 0; i < n; i++) {
	x[i] /= n;
    }

    return (err);
}
