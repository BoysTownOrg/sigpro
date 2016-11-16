// sptest.c - test sigpro library functions

#include <stdio.h>
#include <math.h>
#include "sigpro.h"

#define NP 1024

static void
linspace(float *x, int n, double a, double b)
{
    int     i;
    double  d;

    d = (b - a) / (n - 1);
    for (i = 0; i < n; i++) {
        x[i] = a + i * d;
    }
}

static void
cdb(float *x, float *y, int n)
{
    double  a;
    int     i, j, k, n2;

    for (i = 0; i < (n / 2); i++) {
	j = i * 2;
	k = j + 1;
	a = hypot(x[j], x[k]);
	if (a < 1e-20) {
	    y[i] = -400;
	} else {
	    y[i] = 20 * log10(a);
	}
    }
}

void
sp_copy(float *x, float *y, int n)
{
    int     i;

    for (i = 0; i < n; i++) {
	y[i] = x[i];
    }
}

void
wr_spec(float *x, int n, double fs, char *fn)
{
    float   f[n + 1], y[n + 1];
    int     i, n1, n2;
    FILE   *fp;

    sp_copy(x, y, n);
    sp_rcfft(y, n);
    n1 = n + 1;
    linspace(f, n1, 0, fs);
    n2 = n / 2;
    cdb(y, y, n2);

    fp = fopen(fn, "wt");
    fprintf(fp, "; %s\n", fn);
    for (i = 0; i < n2; i++) {
	fprintf(fp, "%12.4g %12.4g\n", f[i], y[i]);
    }
    fclose(fp);
}

int
main(int ac, char **av)
{
    float   x[NP];
    int     i;
    static double fs = 20000;

    sp_randn(x, NP); 
    wr_spec(x, NP, fs, "spec1.txt"); 
    wr_spec(x, NP, fs, "spec2.txt");

    return(0);
}
