// tst_fft.c - test real & complex FFT and inverse FFT

#include <stdio.h>
#include "sigpro.h"

#define NP 32

int
main()
{
    int     i, j, m, n, c;
    float   x[NP], b[NP+2];
    FILE   *fp;

    c = 4;  // cols
    fp = fopen("tst_fft.txt", "wt");

// store random numbers  in b

    sp_rand(x, NP);
    for (i = 0; i < NP; i++) {;
	b[i] = x[i];
    }

    n = NP;
    m = n / 2 + 1;

    fprintf(fp, "\n\nReal FFT Test (n=%d):\n", n);
    fprintf(fp, "\nreal input sequence\n");
    for (i = 0; i < n; i = i + c) {
	for (j = i; j < i + c; j++)
	    fprintf(fp, "%17.8f", b[j]);
	fprintf(fp, "\n");
    }
    sp_rcfft(b, n);
    fprintf(fp, "\nreal components of transform\n");
    for (i = 0; i < m; i = i + c) {
	for (j = i; j < i + c && j < m; j++)
	    fprintf(fp, "%17.8f", b[2 * j]);
	fprintf(fp, "\n");
    }
    fprintf(fp, "\nimag components of transform\n");
    for (i = 0; i < m; i = i + c) {
	for (j = i; j < i + c && j < m; j++)
	    fprintf(fp, "%17.8f", b[2 * j + 1]);
	fprintf(fp, "\n");
    }
    sp_crfft(b, n);
    fprintf(fp, "\nreal inverse transform\n");
    for (i = 0; i < n; i = i + c) {
	for (j = i; j < i + c; j++)
	    fprintf(fp, "%17.8f", b[j]);
	fprintf(fp, "\n");
    }


    n = NP / 2;

    fprintf(fp, "\n\nComplex FFT Test (n=%d):\n", n);
    fprintf(fp, "\nreal input sequence\n");
    for (i = 0; i < n; i = i + c) {
	for (j = i; j < i + c && j < n; j++)
	    fprintf(fp, "%17.8f", b[2 * j]);
	fprintf(fp, "\n");
    }
    fprintf(fp, "\nimag input sequence\n");
    for (i = 0; i < n; i = i + c) {
	for (j = i; j < i + c && j < n; j++)
	    fprintf(fp, "%17.8f", b[2 * j + 1]);
	fprintf(fp, "\n");
    }
    sp_fft(b, n);
    fprintf(fp, "\nreal components of transform\n");
    for (i = 0; i < n; i = i + c) {
	for (j = i; j < i + c && j < n; j++)
	    fprintf(fp, "%17.8f", b[2 * j]);
	fprintf(fp, "\n");
    }
    fprintf(fp, "\nimag components of transform\n");
    for (i = 0; i < n; i = i + c) {
	for (j = i; j < i + c && j < n; j++)
	    fprintf(fp, "%17.8f", b[2 * j + 1]);
	fprintf(fp, "\n");
    }
    sp_ifft(b, n);
    fprintf(fp, "\nreal inverse transform\n");
    for (i = 0; i < m; i = i + c) {
	for (j = i; j < i + c && j < m; j++)
	    fprintf(fp, "%17.8f", b[2 * j]);
	fprintf(fp, "\n");
    }
    fprintf(fp, "\nimag inverse transform\n");
    for (i = 0; i < m; i = i + c) {
	for (j = i; j < i + c && j < m; j++)
	    fprintf(fp, "%17.8f", b[2 * j + 1]);
	fprintf(fp, "\n");
    }

    fclose(fp);
    return (0);
}
