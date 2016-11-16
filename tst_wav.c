// tst_wav.c - test WAV file read & write

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sigpro.h"
#ifdef WIN32
#include <io.h>
#else
#define _strdup strdup
#endif

void
print_wav(char *fn, VAR *vl, float *fs)
{
    float *f;
    int i, j, nr, nc;
    FILE *fp;

    if (vl[0].data == NULL) {
	return;
    }
    f = (float *) vl[0].data;
    nr = vl[0].rows;
    nc = vl[0].cols;
    fp = fopen(fn, "w");
    fprintf(fp, "; %s: nchn=%d\n", fn, nc);
    fprintf(fp, "rate=%.0f : size=%d\n", *fs, nr);
    for (i = 0; i < nr; i++) {
	for (j = 0; j < nc; j++) {
	    fprintf(fp, " %9.6f", f[i * nc + j]);
	}
	fprintf(fp, "\n");
    }
    fclose(fp);
}

/***********************************************************/

static VAR *
tone_gen(double f, int n, float *fs)
{
    int i;
    double p;
    float *s;
    VAR *vl;
    static double tpi = 2 * M_PI;

    s = (float *) calloc(n, sizeof(float));
    p = tpi * f / *fs;
    for (i = 0; i < n; i++) {
	s[i] = (float) sin(i * p);
    }
    vl = sp_var_alloc(1);
    vl[0].name = _strdup("tone");
    vl[0].data = s;
    vl[0].rows = n;
    vl[0].cols = 1;
    vl[0].dtyp = SP_DTYP_F4;
    vl[0].last = 1;

    return (vl);
}

/***********************************************************/

int
main(int ac, char *av[])
{
    float fs, ft;
    int nb, np;
    VAR *vl;
    static char *ifn = "test/grab.wav";
    static char *ofn = "tst_wav.txt";
    static char *tfn = "tone.wav";

    // test wav_read
    vl = sp_wav_read(ifn, 0, 0, &fs);
    if (vl == NULL) {
	fprintf(stderr, "can't open %s\n", ifn);
	return (1);
    }
    print_wav(ofn, vl, &fs); 

    // test wav_write
    fs = 32000;	        // sampling rate
    ft = 1000;	        // tone frequency
    np = 1024;	        // number of points
    nb = 16;	        // number of bits
    vl = tone_gen(ft, np, &fs);
    sp_wav_write(tfn, vl, &fs, nb);

    sp_var_clear_all();
    return (0);
}

