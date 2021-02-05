// anafilt.c

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sigpro.h"

#define dcopy(x,y,n)    memcpy(x,y,(n)*sizeof(double))
#define dzero(x,n)      memset(x,0,(n)*sizeof(double))

/***********************************************************/

// bilinear_pole - transform analog pole to IIR pole
static void
bilinear_pole(float* p, double* ap, double wp)
{
    double aa, bb, c1, c2, c3, p1, p2;

    aa = 2 * ap[0];
    bb = ap[0] * ap[0] + ap[1] * ap[1];
    c1 = 1 - aa + bb;
    c2 = 2 * (bb - 1);
    c3 = 1 + aa + bb;
    p1 = -c2 / (2 * c1);
    if (ap[1] == 0) {
        p[0] = (float)p1;
        p[1] = 0;
    }
    else {
        p2 = sqrt(c3 / c1 - p1 * p1);
        p[0] = (float)p1;
        p[1] = (float)p2;
        p[2] = (float)p1;
        p[3] = (float)-p2;
    }
}

static double
gain(float* z, float* p, int nz, double w)
{
    double f[2], x[2], y[2], xr, xi, yr, yi, xm, ym, temp;
    int j, jr, ji;

    f[0] = cos(M_PI * w);
    f[1] = sin(M_PI * w);
    x[0] = y[0] = 1;
    x[1] = y[1] = 0;
    for (j = 0; j < nz; j++) {
        jr = j * 2;
        ji = jr + 1;
        xr = f[0] - z[jr];
        xi = f[1] - z[ji];
        yr = f[0] - p[jr];
        yi = f[1] - p[ji];
        temp = x[0] * xr - x[1] * xi;
        x[1] = x[0] * xi + x[1] * xr;
        x[0] = temp;
        temp = y[0] * yr - y[1] * yi;
        y[1] = y[0] * yi + y[1] * yr;
        y[0] = temp;
    }
    xm = sqrt(x[0] * x[0] + x[1] * x[1]);
    ym = sqrt(y[0] * y[0] + y[1] * y[1]);
    temp = ym / xm;
    return (temp);
}

// pole2zp - transform analog pole to IIR zeros, poles, and gain
static void
pole2zp(float* z, float* p, float* g, double* ap, int np, float* wn, int ft)
{
    double bw, u0, u1, wc, wp, Q, M, A, zz[2], pp[2], M1[2], M2[2], N[2];
    float  p1[4], p2[4], z1[4];
    int j, m = 4;

    if ((ft == 0) || (ft == 1)) {
        u0 = tan(M_PI * wn[0] / 2);
        pp[0] = ap[0] * u0;
        pp[1] = ap[1] * u0;
        wp = (ft == 0) ? 1 : -1;
        bilinear_pole(p, pp, wp);
        z[0] = (float)-wp;
        z[1] = 0;
        if (ap[1] == 0) {
            g[0] = (float)fabs(wp - p[0]) / 2;
        }
        else {
            z[2] = z[0];
            z[3] = z[1];
            g[0] = (float)((wp - p[0]) * (wp - p[0]) + p[1] * p[1]) / 4;
        }
    }
    else {
        u0 = tan(M_PI * wn[0] / 2);
        u1 = tan(M_PI * wn[1] / 2);
        bw = u1 - u0;        // bandwidth
        wc = sqrt(u1 * u0);  // center frequency
        if (ft == 2) {
            wp = 1;
            z1[0] = 1;
            z1[1] = 0;
            z1[2] = -1;
            z1[3] = 0;
        }
        else {
            wp = -1;
            zz[0] = 0;
            zz[1] = wc;
            bilinear_pole(z1, zz, wp);
        }
        Q = wc / bw;
        M1[0] = (ap[0] / Q) / 2;
        M1[1] = (ap[1] / Q) / 2;
        N[0] = M1[0] * M1[0] - M1[1] * M1[1] - 1;
        N[1] = M1[0] * M1[1] + M1[1] * M1[0];
        M = sqrt(sqrt(N[0] * N[0] + N[1] * N[1]));
        A = atan2(N[1], N[0]) / 2;
        M2[0] = M * cos(A);
        M2[1] = M * sin(A);
        u0 = tan(M_PI * wn[0] / 2);
        pp[0] = (M1[0] + M2[0]) * wc;
        pp[1] = (M1[1] + M2[1]) * wc;
        bilinear_pole(p1, pp, wp);
        pp[0] = (M1[0] - M2[0]) * wc;
        pp[1] = (M1[1] - M2[1]) * wc;
        bilinear_pole(p2, pp, wp);
        if (ap[1] == 0) {
            for (j = 0; j < m; j++) {
                z[j] = z1[j];
                p[j] = p1[j];
            }
        }
        else {
            for (j = 0; j < m; j++) {
                z[j] = z1[j];
                p[j] = p1[j];
                z[j + m] = z1[j];
                p[j + m] = p2[j];
            }
        }
        wc = (ft == 2) ? sqrt(wn[0] * wn[1]) : 0;
        g[0] = (float)gain(z, p, m, wc);
    }
}

// ap2zp - transform analog prototype to IIR zeros, poles, and gain
static void
ap2zp(float* z, float* p, float* g, double* ap, int np, float* wn, int ft)
{
    double gg;
    int j, jm, m;

    m = ((ft == 0) || (ft == 1)) ? 1 : 2;
    gg = 1;
    for (j = 0; j < np; j += 2) {
        jm = j * m * 2;
        pole2zp(z + jm, p + jm, g, ap + j * 2, np, wn, ft);
        gg *= g[0];
    }
    *g = (float)gg;
}

/***********************************************************/

// transform polynomial roots to coefficients
static void
root2poly(float* r, double* p, int n)
{
    double* pp, * qq;
    int i, ir, ii, j, jr, ji;

    pp = (double*)calloc((n + 1) * 2, sizeof(double));
    qq = (double*)calloc((n + 1) * 2, sizeof(double));
    dzero(pp, (n + 1) * 2);
    dzero(qq, (n + 1) * 2);
    pp[0] = qq[0] = 1;
    for (i = 0; i < n; i++) {
        ir = i * 2;
        ii = i * 2 + 1;
        qq[2] = pp[2] - r[ir];
        qq[3] = pp[3] - r[ii];
        for (j = 0; j < i; j++) {
            jr = j * 2;
            ji = j * 2 + 1;
            qq[jr + 4] = pp[jr + 4] - (pp[jr + 2] * r[ir] - pp[ji + 2] * r[ii]);
            qq[ji + 4] = pp[ji + 4] - (pp[ji + 2] * r[ir] + pp[jr + 2] * r[ii]);
        }
        dcopy(pp, qq, (n + 1) * 2);
    }
    // return real part of product-polynomial coefficients
    for (i = 0; i < (n + 1); i++) {
        p[i] = pp[i * 2];
    }
    free(pp);
    free(qq);
}

// transform filterbank poles and zeros to IIR coefficients
static void
zp2ba(float* z, float* p, int nz, float* b, float* a)
{
    double *bd, *ad;
    int i;

    bd = (double*)calloc(nz + 1, sizeof(double));
    ad = (double*)calloc(nz + 1, sizeof(double));
    root2poly(z, bd, nz);
    root2poly(p, ad, nz);
    for (i = 0; i <= nz; i++) {
        b[i] = (float)bd[i];
        a[i] = (float)ad[i];
    }
    free(bd);
    free(ad);
}

// lp2bx - lp2bx band-pass or band-stop
static void
lp2bx(float* lp, int np, float* b, float* a, float* wn, int ft)
{
    double *ap;
    float *z, *p, g[1];
    int i, ir, ii, m;

    // allocate local arrays
    ap = (double*)calloc(np * 2, sizeof(double));
    z = (float*)calloc(np * 4, sizeof(float));
    p = (float*)calloc(np * 4, sizeof(float));
    // copy low-pass prototype to analog prototype
    m = np / 2;
    for (i = 0; i < m; i++) {
        ir = 2 * i;
        ii = ir + 1;
        ap[ir] = lp[ir];
        ap[ii] = lp[ii];
        ap[2 * m + ir] = lp[ir];
        ap[2 * m + ii] = -lp[ii];
    }
    if (np % 2) { // copy real pole, if any
        ap[2 * np - 2] = lp[np - 1];
        ap[2 * np - 1] = 0;
    }
    // transform analog prototype to zeros and poles
    ap2zp(z, p, g, ap, np, wn, ft);
    // transform poles & zeros to IIR coeficients
    if (ft > 1) np *= 2;
    zp2ba(z, p, np, b, a);
    for (i = 0; i <= np; i++) {
        b[i] *= g[0];
    }
    // free local arrays
    free(ap);
    free(z);
    free(p);
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

FUNC(void) sp_bessel(   // Bessel filter design
    float *b,           // input coeffcients
    float *a,           // output coeffcients
    int n,              // filter order
    float *wn,          // cutoff frequency
    int ft              // filter type 
)
{
    float *p;
    
    p = (float *) calloc(n, sizeof(float));
    besselp(p, n);
    lp2bx(p, n, b, a, wn, ft);
    free(p);
}

FUNC(void) sp_butter(	// Butterworth filter design
    float *b,           // input coeffcients
    float *a,           // output coeffcients
    int n,              // filter order
    float *wn,          // cutoff frequency
    int ft              // filter type 
)
{
    float *p;
    
    p = (float *) calloc(n, sizeof(float));
    butterp(p, n);
    lp2bx(p, n, b, a, wn, ft);
    free(p);
}


FUNC(void) sp_cheby(    // Chebyshev filter design
    float *b,           // input coeffcients
    float *a,           // output coeffcients
    int n,              // filter order
    float *wn,          // cutoff frequency
    int ft,             // filter type 
    double rip          // pass-band ripple
)
{
    float *p;
    
    p = (float *) calloc(n, sizeof(float));
    chebyp(p, n, rip);
    lp2bx(p, n, b, a, wn, ft);
    free(p);
}
