// tst_xfr.c - test transfer-function computaton

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sigpro.h"

/**********************************************************/

void
freq_resp(float *H, int nf, double df, char *fn)
{
    double f;
    float *db, *ph, *gd;
    int i;
    FILE *fp;
    
    db = (float *) calloc(nf, sizeof(float));
    ph = (float *) calloc(nf, sizeof(float));
    gd = (float *) calloc(nf, sizeof(float));
    sp_cdb(H, db, nf);          // decibels
    sp_cph(H, ph, nf);          // phase
    sp_unwrap(ph, ph, nf);      // unwrap phase
    sp_cgd(H, gd, nf, df);      // delay
    
    fp = fopen(fn, "w");
    if (!fp) {
	printf("can't open %s\n", fn);
	return;
    }
    fprintf(fp, "; %s - transfer-function test test\n\n", fn);
    fprintf(fp, ";  f        dB       ph       gd\n");
    for (i = 0; i < nf; i++) {
	f = i * df;
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
    double df, fc, sr;
    float wn[2];
    float *a, *b, *x, *y, *z;
    int ft, no, np, nc, nf;

    sr = 48000;         // sampling rate
    fc = 8000;	        // low-pass cut-off frequency
    ft = 0;	        // filter type: 0=low_pass, 1=high_pass
    no = 4;	        // filter order
    np = 8192;	        // number of samples = buffer length
    nc = no + 1;        // number of filter coefficients
    nf = np / 2 + 1;    // number of frequencies

    a = (float *) calloc(nc, sizeof(float));
    b = (float *) calloc(nc, sizeof(float));
    x = (float *) calloc(np, sizeof(float));
    y = (float *) calloc(np, sizeof(float));
    z = (float *) calloc(nf * 2, sizeof(float));

    // compute white noise
    sp_randflat(x, np);                 // noise stimulus in x

    // Butterworth low-pass filter
    wn[0] = (float)(2 * fc / sr);       // cutoff_frequency / Nyquist_rate
    sp_butter(b, a, no, wn, ft);        // Butterworth filter coefficients
    sp_filter(b, nc, a, nc, x, y, np);  // filter response in y

    // compute transfer function
    sp_transfer(x, y, np, z);           // transfer function in z

    // write results
    df = (sr / np) / 1000;              // frequency resolution (kHz)
    freq_resp(z, nf, df, "transfer.txt");
    
    free(a);
    free(b);
    free(x);
    free(y);
    free(z);

    return (0);
}
