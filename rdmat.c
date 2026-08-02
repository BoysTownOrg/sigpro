/* rdmat.c  - reads MAT data files and prints contents */

#include <stdio.h>
#include <stdlib.h>
#ifdef WIN32
#include <io.h>
#else
#include <unistd.h>
#define _access access
#endif
#include "sigpro.h"

static int bflg = 0;
static int ntxt = 50;
static int nprn = 4;
static int prnv = 0;
static int pflg = 0;

/* rdmat - read MAT file */

void
rdmat(char *fn)
{
    char *vn, *dt, c;
    float *f4;
    int i, j, nv, nd, np, nr, nc, nx;
    VAR *vl, *vc;
    static char mod;
    static short nrc[2] = {8, 8};

    if (_access(fn, 0)) {
	printf("%s: can't open\n", fn);
	return;
    }
    nv = sp_mat_size(fn);
    if (nv < 1) {
	printf("%s: can't read\n", fn);
	return;
    }
    printf("%s (v%d): %d variable%s\n", fn, sp_mat_version(fn),
	nv, (nv == 1) ? "" : "s");
    if (bflg) {		// brief list ?
	return;
    }
    vl = sp_mat_whos(fn);
    if (prnv) {
	printf("     name   rows   cols type        data\n");
    } else {
	printf("     name   rows   cols type\n");
    }
    for (i = 0; i < nv; i++) {
	vn = vl[i].name;
	nr = vl[i].rows;
	nc = vl[i].cols;
	nx = vl[i].cmpx ? 2 : 1;
	dt = sp_var_dattyp(vl[i].dtyp);
	if (vl[i].dtyp == 15) {
	    printf("  %s: %d bytes compressed\n", vn, nr);
	    continue;
	}
	if (vl[i].text) {
	    mod = 't';
	} else if (vl[i].cmpx) {
	    mod = 'c';
	} else {
	    mod = ' ';
	}
	printf("%9s %6d %6d  %s%c", vn, nr, nc, dt, mod);
	if (prnv) {
	    vc = sp_mat_fetch(fn, vn, NULL, nrc);
	    sp_var_float(vc);
	    nd = nr * nc * nx;
	    f4 = (float *) vc[0].data;
	    if (vl[i].text) {
		np = (nd < ntxt) ? nd : ntxt;
		printf(" ");
		for (j = 0; j < np; j++) {
                    c = ((f4[j] < 32) || (f4[j] > 255)) ? '.' : (char) f4[j];
		    printf("%c", c);
		}
		printf("%s\n", np >= nd ? "" : "...");
	    } else {
		np = (nd < nprn) ? nd : nprn;
		for (j = 0; j < np; j++) {
		    printf("%12.5g", f4[j]);
		}
		printf("%s\n", np >= nd ? "" : " ...");
	    }
	    sp_var_clear(vc);
	} else {
	    printf("\n");
	}
    }
    sp_var_clear_all();
}

int
main(int ac, char **av)
{
    if (ac < 2) {
	printf("usage: rdmat [options] file.mat ...\n");
	printf("options:\n");
	printf("  -b    brief list\n");
	printf("  -p    wait for Enter\n");
	printf("  -v    display variable values\n");
	_exit(0);
    }
    while (ac > 1) {
        if (av[1][0] == '-') {
	    if (av[1][1] == 'b') {
		bflg = 1;
	    } else if (av[1][1] == 'p') {
		pflg = 1;
	    } else if (av[1][1] == 'v') {
		prnv = 1;
	        if (av[1][2] >= '1'
		    && av[1][2] <= '9') {
		    nprn = av[1][2] - '0';
	        }
	    }
	} else {
	    rdmat(av[1]);
	}
	ac--;
	av++;
    }
    if (pflg) {
	printf("\n[Press Enter]\n");
	getchar();
    }
    return (0);
}
