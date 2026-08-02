// tst_shp.c - test frequency-shaping functions

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sigpro.h"

#define NT 6
#define NF 250
#define NP 2040

void
wr_spec(float *x, int n, double fs, char *fn)
{
    float  *f, *y;
    int     i, n1;
    FILE   *fp;

    n1 = n / 2 + 1;
    f = (float *) calloc(n + 2, sizeof(float));
    y = (float *) calloc(n + 2, sizeof(float));
    sp_linspace(f, n1, 0, fs / 2);
    sp_copy(x, y, n);
    sp_rcfft(y, n);
    sp_cdb(y, y, n1);

    fp = fopen(fn, "wt");
    fprintf(fp, "; %s\n", fn);
    for (i = 0; i < n1; i++) {
	fprintf(fp, "%12.4g %12.4g\n", f[i], y[i]);
    }
    fclose(fp);

    free(f);
    free(y);
}

int
main(int ac, char **av)
{
    static double fs = 20000;
    static float x[NP+2];
    static float fr[NT] = {0, 2000, 4000, 6000, 8000, 10000};
    static float at[NT] = {0,    0,   40,   40,    0,     0};

    // frequency-shape FIR

    sp_firdb(x, NF, fs, fr, at, NT); 
    wr_spec(x, NF, fs, "spec1.txt");

    // white noise wave form

    sp_randflat(x, NP); 
    wr_spec(x, NP, fs, "spec2.txt"); 

    // filtered waveform

    sp_frqshp(x, x, NP, NF, fs, fr, at, NT, 1); 
    wr_spec(x, NP, fs, "spec3.txt");

    return(0);
}
