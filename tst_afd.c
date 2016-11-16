// tst_afd.c - test analog-filter design

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sigpro.h"

/**********************************************************/

void
impl_resp(float *y, int np, char *fn)
{
    int i;
    FILE *fp;
    
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
freq_resp(float *y, int np, char *fn)
{
    double f, df;
    float *db, *ph, *gd;
    int i, nf;
    FILE *fp;
    
    db = (float *) calloc(np+2, sizeof(float));
    ph = (float *) calloc(np+2, sizeof(float));
    gd = (float *) calloc(np+2, sizeof(float));
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
	f = (double) i / nf;
	fprintf(fp, "%8.3f %8.3f %8.3f %8.3f\n", f, db[i], ph[i], gd[i]);
    }
    fclose(fp);
    free(db);
    free(ph);
    free(gd);
}

int
main()
{
    double rip;
    float wn[2];
    float *a, *b, *x, *y;
    int no, np, nc, ft;

    ft = 0;	    // filter type: 0=low_pass, 1=high_pass
    no = 3;	    // filter order
    np = 1024;	    // number of samples
    nc = no + 1;
    wn[0] = 0.1f;   // cutoff_frequency / half_sample_rate

    a = (float *) calloc(nc, sizeof(float));
    b = (float *) calloc(nc, sizeof(float));
    x = (float *) calloc(np+2, sizeof(float));
    y = (float *) calloc(np+2, sizeof(float));
    x[0] = 1;
    
    // Butterworth filter
    sp_butter(b, a, no, wn, ft);	// filter coefficients
    sp_filter(b, nc, a, nc, x, y, np);	// impulse response
    impl_resp(y, np, "butteri.txt");
    freq_resp(y, np, "butterf.txt");
    
    // Bessel filter
    sp_bessel(b, a, no, wn, ft);	// filter coefficients
    sp_filter(b, nc, a, nc, x, y, np);	// impulse response
    impl_resp(y, np, "besseli.txt");
    freq_resp(y, np, "besself.txt");
    
    // Chebyshev filter
    rip = 3;				// pass-band ripple
    sp_cheby(b, a, no, wn, ft, rip);	// filter coefficients
    sp_filter(b, nc, a, nc, x, y, np);	// impulse response
    impl_resp(y, np, "chebyi.txt");
    freq_resp(y, np, "chebyf.txt");
    
    free(a);
    free(b);
    free(x);
    free(y);

    return (0);
}

