// anafilt.c

#include <stdlib.h>
#include <math.h>
#include "sigpro.h"

static void lp2bax(float *, int, float *, float *, float *, int);

// poly_mul - multiply two polynomials
static void
poly_mul(float *p1, int n1, float *p2, int n2, float *p3)
{
    float *pp;
    int i, j, n3;

    // allocate temp array for polynomial product
    n3 = n1 + n2 - 1;
    pp = (float *) calloc(n3, sizeof(float));
    // multiply input polynomials
    for (i = 0; i < n3; i++) {
	pp[i] = 0;
    }
    for (i = 0; i < n1; i++) {
	for (j = 0; j < n2; j++) {
	    pp[i + j] += p1[i] * p2[j];
	}
    }
    // copy product polynomial to output
    for (i = 0; i < n3; i++) {
	p3[i] = pp[i];
    }
    free(pp);
}

// poly_add - add two polynomials
static void
poly_add(float *p1, float *p2, float *p3, int n)
{
    int i;

    for (i = 0; i < n; i++) {
	p3[i] = p1[i] + p2[i];
    }
}

// sp2ba - single-pole bilinear transformation
static void
sp2ba(float *p, float *b, float *a)
{
    double  aa;
    
    aa = p[0];
    a[0] = (float) (aa - 1);
    a[1] = (float) (aa + 1);
    b[0] = (float) aa;
    b[1] = (float) aa;
}

// dp2ba - double-pole bilinear transformation
static void
dp2ba(float *p, float *b, float *a)
{
    double  aa, bb;
    
    aa = 2 * p[0];
    bb = p[0] * p[0] + p[1] * p[1];
    a[0] = (float) (1 - aa + bb);
    a[1] = (float) (2 * (bb - 1));
    a[2] = (float) (1 + aa + bb);
    b[0] = (float) bb;
    b[1] = (float) (2 * bb);
    b[2] = (float) bb;
}

// warp - all-pole low-pass frequency warp
static void
warp(float *p, int np, double wn)
{
    double  rr;
    int i;
    
    rr = tan(M_PI * wn / 2);
    for (i = 0; i < np; i++) {
	p[i] = (float) (p[i] * rr);
    }
}

// lp2pa - all-pole low-pass bilinear transformation
static void
lp2ba(float *p, int np, float *b, float *a, float *wn, int ft)
{
    float aa[4], bb[4];
    double  cc, s, suma, sumb;
    int i, nc;
    
    if (ft < 0 || ft > 1) {
	lp2bax(p, np, b, a, wn, ft);
	return;
    }
    warp(p, np, *wn);
    nc = 1;
    a[0] = b[0] = 1;
    for (i = 1; i < (np + 1); i++) {
	a[i] = b[i] = 0;
    }
    for (i = 0; i < (np - 1); i += 2) {
        dp2ba(p + i, bb, aa);
	if (ft) bb[1] = -bb[1];
	poly_mul(b, nc, bb, 3, b);
	poly_mul(a, nc, aa, 3, a);
	nc +=  2;
    }    
    if (np % 2) {
        sp2ba(p + np - 1, bb, aa);
	if (ft) bb[1] = -bb[1];
	poly_mul(b, nc, bb, 2, b);
	poly_mul(a, nc, aa, 2, a);
	nc += 1;
    }
    // normalize so that a[0] = 1
    nc = np + 1;
    cc = a[0];
    s = 1;
    suma = sumb = 0;
    for (i = 0; i < nc; i++) {
        a[i] = (float) (a[i] / cc);
        b[i] = (float) (b[i] / cc);
        suma += a[i] * s;
        sumb += b[i] * s;
        if (ft) s = -s;
    }
    cc = sumb / suma;
    for (i = 0; i < nc; i++) {
        b[i] = (float) (b[i] / cc);
    }
}

// lp2pax - lp2ba extension to do band-pass & band-stop
static void
lp2bax(float *p, int np, float *b, float *a, float *wn, int ft)
{
    int i, nc, ns;
    float *p1, *b1, *a1, *w1, *c1;
    float *p2, *b2, *a2, *w2, *c2;

    if (ft < 2 || ft > 3) {
	return;
    }
    nc = 1 + np;
    ns = 1 + np * 2;
    p1 = p;
    p2 = (float *) calloc(nc, sizeof(float));
    a1 = (float *) calloc(nc, sizeof(float));
    a2 = (float *) calloc(nc, sizeof(float));
    b1 = (float *) calloc(nc, sizeof(float));
    b2 = (float *) calloc(nc, sizeof(float));
    c1 = (float *) calloc(ns, sizeof(float));
    c2 = (float *) calloc(ns, sizeof(float));
    w1 = wn;
    w2 = wn + 1;
    for (i = 0; i < np; i++) {
	p2[i] = p1[i];
    }
    if (ft == 2) {		// band-pass
	lp2ba(p1, np, b1, a1, w1, 1);
	lp2ba(p2, np, b2, a2, w2, 0);
	poly_mul(a1, nc, a2, nc, a);
	poly_mul(b1, nc, b2, nc, b);
    } else {			// band-stop
	lp2ba(p1, np, b1, a1, w1, 0);
	lp2ba(p2, np, b2, a2, w2, 1);
	poly_mul(a1, nc, a2, nc, a);
	poly_mul(b1, nc, a2, nc, c1);
	poly_mul(a1, nc, b2, nc, c2);
	poly_add(c1, c2, b, ns);
    }
    free(p2);
    free(a1);
    free(a2);
    free(b1);
    free(b2);
    free(c1);
    free(c2);
}

/**********************************************************/

static void
butterp(float *p, int n)
{
    double aa;
    int i, ir, ii, m;
    
    if (n < 1) {
	return;
    }
    m = n / 2;
    for (i = 0; i < m; i++) {
	ir = 2 * i;
	ii = ir + 1;
	aa = ii * M_PI / (2 * n);
	p[ir] = (float) (-sin(aa));
	p[ii] = (float) (cos(aa));
    }
    if (n % 2) {
	p[n - 1] = -1;
    }
}

static void
besselp(float *p, int n)
{
    int i;
    static double p01[] = {-1.000000000};
    static double p02[] = {-0.866025404, 0.500000000};
    static double p03[] = {-0.745640386, 0.711366625,-0.941600027};
    static double p04[] = {-0.657211172, 0.830161435,-0.904758797, 0.270918733};
    static double p05[] = {-0.590575945, 0.907206756,-0.851553619, 0.442717464,
	-0.926442077};
    static double p06[] = {-0.538552682, 0.961687688,-0.799654186, 0.562171735,
	-0.909390683, 0.185696440};
    static double p07[] = {-0.496691726, 1.002508508,-0.752735543, 0.650469631,
	-0.880002934, 0.321665276,-0.919487156};
    static double p08[] = {-0.462174041, 1.034388681,-0.711138181, 0.718651731,
	-0.847325080, 0.425901754,-0.909683155, 0.141243798};
    static double p09[] = {-0.433141556, 1.060073670,-0.674362269, 0.773054621,
	-0.814802111, 0.508581569,-0.891121702, 0.252658093,-0.915495780};
    static double p10[] = {-0.408322073, 1.081274843,-0.641751387, 0.817583617,
	-0.783769441, 0.575914754,-0.868845964, 0.343000823,
	-0.909134732, 0.113958314};
    static double p11[] = {-0.386814951, 1.099117467,-0.612687155, 0.854781389,
	-0.754693893, 0.631915005,-0.845304401, 0.417869692,
	-0.896365671, 0.208048038,-0.912906724};
    static double p12[] = {-0.367964009, 1.114373576,-0.586636932, 0.886377275,
	-0.727668162, 0.679296118,-0.821729694, 0.481021212,
	-0.880253434, 0.287177950,-0.908447823, 0.095506365};
    static double p13[] = {-0.351279232, 1.127591548,-0.563155984, 0.913590034,
	-0.702623468, 0.719961189,-0.798746069, 0.535075212,
	-0.862509420, 0.354741373,-0.899131467, 0.176834296,
	-0.911091467};
    static double p14[] = {-0.336386822, 1.139172298,-0.541876678, 0.937304368,
	-0.679425643, 0.755285731,-0.776659139, 0.581917068,
	-0.844119916, 0.413165383,-0.886950667, 0.247007918,
	-0.907793214, 0.082196399};
    static double p15[] = {-0.322996306, 1.149416155,-0.522495407, 0.958178726,
	-0.657919659, 0.786289550,-0.755602717, 0.622939636,
	-0.825663145, 0.464234875,-0.873126462, 0.308235247,
	-0.900698169, 0.153768120,-0.909748236};
    static double p16[] = {-0.310878276, 1.158552841,-0.504760644, 0.976713748,
	-0.637950251, 0.813745354,-0.735616630, 0.659195088,
	-0.807479029, 0.509293375,-0.858426423, 0.362169727,
	-0.891172307, 0.216708966,-0.907209960, 0.072142113};
    static double *pp[] = {p01, p02, p03, p04, p05, p06, p07, p08,
	p09, p10, p11, p12, p13, p14, p15, p16};

    if (n < 1 || n > 16) {
	return;
    }
    for (i = 0; i < n; i++) {
	p[i] = (float) ((pp[n-1])[i]);
    }
}

static void
chebyp(float *p, int n, double rip)
{
    double aa, eps, mu;
    int i, ir, ii, m;
    
    if (n < 1) {
	return;
    }
    m = n / 2;
    eps = pow(10.0, rip / 20) - 1;
    mu = log((1 + sqrt(1 + eps * eps)) / eps) / n;   // =asinh(1/eps)/n
    for (i = 0; i < m; i++) {
	ir = 2 * i;
	ii = ir + 1;
	aa = ii * M_PI / (2 * n);
	p[ir] = (float) (-sin(aa) * sinh(mu));
	p[ii] = (float) (cos(aa) * cosh(mu));
    }
    if (n % 2) {
	p[n - 1] = -1;
    }
}

/**********************************************************/

FUNC(void) sp_bessel(	// Bessel filter design
    float *b, 		// input coeffcients
    float *a, 		// output coeffcients
    int n,		// filter order
    float *wn,		// cutoff frequency
    int ft		// filter type 
)
{
    float *p;
    
    p = (float *) calloc(n, sizeof(float));
    besselp(p, n);
    lp2ba(p, n, b, a, wn, ft);
    free(p);
}

FUNC(void) sp_butter(	// Butterworth filter design
    float *b, 		// input coeffcients
    float *a, 		// output coeffcients
    int n,		// filter order
    float *wn,		// cutoff frequency
    int ft		// filter type 
)
{
    float *p;
    
    p = (float *) calloc(n, sizeof(float));
    butterp(p, n);
    lp2ba(p, n, b, a, wn, ft);
    free(p);
}


FUNC(void) sp_cheby(	// Chebyshev filter design
    float *b, 		// input coeffcients
    float *a, 		// output coeffcients
    int n,		// filter order
    float *wn,		// cutoff frequency
    int ft,		// filter type 
    double rip		// pass-band ripple
)
{
    float *p;
    
    p = (float *) calloc(n, sizeof(float));
    chebyp(p, n, rip);
    lp2ba(p, n, b, a, wn, ft);
    free(p);
}
