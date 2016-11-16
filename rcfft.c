// rcfft.c 

#include <stdlib.h>
#include <math.h>
#include "sigpro.h"

//-----------------------------------------------------------

static double tpi = 2 * M_PI;
static double srh = 0.70710678118655;

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
rad4(int ii, int nn,
    float *x0, float *x1, float *x2, float *x3,
    float *x4, float *x5, float *x6, float *x7)
{
    double  arg, tpiovn;
    float   c1, c2, c3, s1, s2, s3, pr, pi, r1, r5;
    float   t0, t1, t2, t3, t4, t5, t6, t7;
    int     i0, i4, j, j0, ji, jl, jr, jlast, k, k0, kl, m, n, ni;

    n = nn / 4;
    for (m = 1; (1 << m) < n; m++)
	continue;
    tpiovn = tpi / nn;
    ji = 3;
    jl = 2;
    jr = 2;
    ni = (n + 1) / 2;
    for (i0 = 0; i0 < ni; i0++) {
	if (i0 == 0) {
	    for (k = 0; k < ii; k++) {
		t0 = x0[k] + x2[k];
		t1 = x1[k] + x3[k];
		x2[k] = x0[k] - x2[k];
		x3[k] = x1[k] - x3[k];
		x0[k] = t0 + t1;
		x1[k] = t0 - t1;
	    }
	    if (nn > 4) {
		k0 = ii * 4;
		kl = k0 + ii;
		for (k = k0; k < kl; k++) {
		    pr = (float) (srh * (x1[k] - x3[k]));
		    pi = (float) (srh * (x1[k] + x3[k]));
		    x3[k] = x2[k] + pi;
		    x1[k] = pi - x2[k];
		    x2[k] = x0[k] - pr;
		    x0[k] += pr;
		}
	    }
	} else {
	    arg = tpiovn * bitrev(i0, m);
	    c1 = (float) cos(arg);
	    s1 = (float) sin(arg);
	    c2 = c1 * c1 - s1 * s1;
	    s2 = c1 * s1 + c1 * s1;
	    c3 = c1 * c2 - s1 * s2;
	    s3 = c2 * s1 + s2 * c1;
	    i4 = ii * 4;
	    j0 = jr * i4;
	    k0 = ji * i4;
	    jlast = j0 + ii;
	    for (j = j0; j < jlast; j++) {
		k = k0 + j - j0;
		r1 = x1[j] * c1 - x5[k] * s1;
		r5 = x1[j] * s1 + x5[k] * c1;
		t2 = x2[j] * c2 - x6[k] * s2;
		t6 = x2[j] * s2 + x6[k] * c2;
		t3 = x3[j] * c3 - x7[k] * s3;
		t7 = x3[j] * s3 + x7[k] * c3;
		t0 = x0[j] + t2;
		t4 = x4[k] + t6;
		t2 = x0[j] - t2;
		t6 = x4[k] - t6;
		t1 = r1 + t3;
		t5 = r5 + t7;
		t3 = r1 - t3;
		t7 = r5 - t7;
		x0[j] = t0 + t1;
		x7[k] = t4 + t5;
		x6[k] = t0 - t1;
		x1[j] = t5 - t4;
		x2[j] = t2 - t7;
		x5[k] = t6 + t3;
		x4[k] = t2 + t7;
		x3[j] = t3 - t6;
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
rcfft2(float *x, int m)
{
    int     ii, nn, m2, it, n;

    n = 1 << m;;
    m2 = m / 2;

// radix 2

    if (m <= m2 * 2) {
	nn = 1;
    } else {
	nn = 2;
	ii = n / nn;
	rad2(ii, x, x + ii);
    }

// radix 4

    if (m2 != 0) {
	for (it = 0; it < m2; it++) {
	    nn = nn * 4;
	    ii = n / nn;
	    rad4(ii, nn, x, x + ii, x + 2 * ii, x + 3 * ii,
		x, x + ii, x + 2 * ii, x + 3 * ii);
	}
    }

// re-order

    reorder1(m, x);
    reorder2(m, x);
    for (it = 3; it < n; it += 2)
	x[it] = -x[it];
    x[n] = x[1];
    x[1] = 0.0;
    x[n + 1] = 0.0;

    return (0);
}

//-----------------------------------------------------------

static int
rcdft(float *x, int n)
{
    double a, ur, ui, vr, vi, wr, wi, xr, yr, yi, zz;
    float *y;
    int i, ir, ii, j, m, k;
    static double tpi = 2 * M_PI;

    k = n / 2;
    m = (k + 1) * 2;
    y = (float *) calloc(m, sizeof(float));
    a = -tpi / n;
    ur = cos(a);
    ui = sin(a);
    vr = 1;
    vi = 0;
    for (i = 0; i <= k; i++) {
        xr = x[0];
	yr = xr * vr;
	yi = xr * vi;
	wr = vr * vr - vi * vi;
	wi = vr * vi + vi * vr;
	for (j = 1; j < n; j++) {
	    xr = x[j];
	    yr += xr * wr;
	    yi += xr * wi;
	    zz = vr * wr - vi * wi;
	    wi = vr * wi + vi * wr;
	    wr = zz;
	}
	ir = 2 * i;
	ii = ir + 1;
	y[ir] = (float) yr;
	y[ii] = (float) yi;
	zz = ur * vr - ui * vi;
	vi = ur * vi + ui * vr;
	vr = zz;
    }
    for (i = 0; i < m; i++) {
	x[i] = y[i];
    }
    free(y);

    return (0);
}


//-----------------------------------------------------------

// sp_rcfft -  real-to-complex FFT

FUNC(int)
sp_rcfft(float *x, int n)
{
    int m, err;

    if (n < 2)
	return (0);

    m = ilog2(n);
    if (m > 0) {
	err = rcfft2(x, m);
    } else {
        err = rcdft(x, n);
    }

    return (err);
}
