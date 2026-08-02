// tst_src.c - test sampling-rate conversion

#include <stdio.h>
#include <math.h>
#include "sigpro.h"

#define N0 550
#define N1 15
#define N2 91
#define NC 5

static double tpi = 2 * M_PI;
static float s0[N0], s1[N1], s2[N2];

void
wrtdat(float *s, int n, char *fn)
{
    double dt;
    int i;
    FILE *fp;

    fp = fopen(fn, "wt");
    dt = 1.0 / n;
    fprintf(fp, "; %s\n", fn);
    for (i = 0; i < n; i++)
	fprintf(fp, "%9.3f %9.3f\n", i * dt, s[i]);
    fclose(fp);
}

int
main(int ac, char **av)
{
    int i;

    for (i = 0; i < N0; i++)
	s0[i] = (float) sin((tpi * NC * i) / N0);
    sp_convert(s0, N0, s1, N1, 0, 1);
    sp_convert(s1, N1, s2, N2, 0, 1);
    wrtdat(s0, N0, "s0.txt");
    wrtdat(s1, N1, "s1.txt");
    wrtdat(s2, N2, "s2.txt");
    return(0);
}
