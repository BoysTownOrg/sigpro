// tst_afd.c - test Butterworth band-pass filter

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sigpro.h"

/**********************************************************/

void
impl_resp(float* y, int np, char* fn)
{
    int i;
    FILE* fp;

    fp = fopen(fn, "w");
    if (!fp) {
	printf("can't open %s\n", fn);
	return;
    }
    for (i = 0; i < np; i++) {
	fprintf(fp, "%6d %9.6f\n", i, y[i]);
    }
    fclose(fp);
}

void
freq_resp(float* y, int np, char* fn)
{
    double f, df;
    float* db, * ph, * gd;
    int i, nf;
    FILE* fp;

    db = (float*)calloc(np + 2, sizeof(float));
    ph = (float*)calloc(np + 2, sizeof(float));
    gd = (float*)calloc(np + 2, sizeof(float));
    df = 1.0 / np;
    sp_rcfft(y, np);		// frequency response
    nf = np / 2;
    sp_cdb(y, db, nf);		// decibels
    sp_cph(y, ph, nf);		// phase
    sp_unwrap(ph, ph, nf);	// unwrap phase
    while (ph[nf / 2] > 0) {
	for (i = 0; i < nf; i++) {
	    ph[i] -= 1;
	}
    }
    sp_cgd(y, gd, nf, df);	// delay

    fp = fopen(fn, "w");
    if (!fp) {
	printf("can't open %s\n", fn);
	return;
    }
    fprintf(fp, "; %s - analog-filter design test\n\n", fn);
    fprintf(fp, ";  f        db       ph       gd\n");
    for (i = 0; i < nf; i++) {
	f = (double)i / nf;
	fprintf(fp, "%8.3f %8.3f %8.3f %8.3f\n", f, db[i], ph[i], gd[i]);
    }
    fclose(fp);
    free(db);
    free(ph);
    free(gd);
}

static void
show_coef(char* s, float* c, int n, double g)
{
    int i;

    printf("%s =", s);
    for (i = 0; i < n; i++) printf(" %8.4f", c[i] / g);
    if (g < 1) {
	printf(" * %g\n", g);
    }
    else {
	printf("\n");
    }
}

int
main()
{
    double gn = 1;
    float *a, *b, *x, *y;
    int no, np, nc, ft;
    static float wn[2] = { 0.029462f, 0.058925f };

    ft = 2;	    // filter type: 0=low_pass, 1=high_pass, 2=band_pass
    no = 2;	    // filter order
    np = 1024;	    // number of samples
    nc = 2 * no + 1;

    a = (float*)calloc(nc, sizeof(float));
    b = (float*)calloc(nc, sizeof(float));
    x = (float*)calloc(np + 2, sizeof(float));
    y = (float*)calloc(np + 2, sizeof(float));
    x[0] = 1;

    // Butterworth bandpass filter
    sp_butter(b, a, no, wn, ft);        // filter coefficients
    sp_filter(b, nc, a, nc, x, y, np);  // impulse response
    freq_resp(y, np, "tst_bbf.txt");

    // report coefficients
    show_coef("b", b, nc, gn);
    show_coef("a", a, nc, 1);

    free(a);
    free(b);
    free(x);
    free(y);

    return (0);
}

