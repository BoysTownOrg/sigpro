// tst_tb.c - test tone-burst waveform calculation

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sigpro.h"

#define round(x) ((int)floor((x)+0.5))

void
print_tb(char *fn, float *x, int np, double sr)
{
    double tms;
    int i;
    FILE *fp;

    fp = fopen(fn, "wt");
    fprintf(fp, "; %s\n", fn);
    fprintf(fp, "rate=%.0f : size=%d\n", sr, np);
    for (i = 0; i < np; i++) {
	tms = i * 1000 / sr;
	fprintf(fp, "%9.3f %9.6f\n", tms, x[i]);
    }
    fclose(fp);
}

/***********************************************************/

void
gen_tb(float *x, int n, double sr,
    double ft, double at, double td, double ph, int wt)
{
    int i, m;
    double p, o;
    float s, a;

    m = round(td * sr / 1000);
    a = (float) pow(10, -at / 20);
    sp_window(x, m, wt);
    p = 2 * M_PI * ft / sr;
    o = 2 * M_PI * ph - p * m / 2; 
    for (i = 0; i < m; i++) {
	s = (float) cos(i * p + o);
        x[i] *= a * s;
    }
    for (i = m; i < n; i++) {
        x[i] = 0;
    }
}

/***********************************************************/

int
main(int ac, char *av[])
{
    float *x;
    static int np = 1024;	// number of points
    static int wt = 4;		// window type (4=Blackman)
    static double ft = 1000; 	// tone frequency (Hz)
    static double at = 0; 	// tone attenuation (dB)
    static double ph = 0; 	// tone center phase (cyc)
    static double td = 8;	// tone duration (ms)
    static double sr = 32000;	// sampling rate (Hz)

    x = (float *) calloc(np, sizeof(float));
    gen_tb(x, np, sr, ft, at, td, ph, wt); 
    print_tb("tst_tb.txt", x, np, sr); 
    free(x);
    printf("%s\n", sp_version());

    return (0);
}

