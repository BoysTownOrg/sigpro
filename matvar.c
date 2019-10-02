// matvar.c - MATLAB variable and MAT file functions

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "sigpro.h"
#ifdef WIN32
#include <io.h>
#else
#include <unistd.h>
#define _strdup strdup
#endif

#define MAXNCH    256
#define free_null(x)  if(x){free(x);x=NULL;}
#ifndef O_BINARY
#define O_BINARY  0
#endif
#define PMODE       (S_IREAD|S_IWRITE)
#define AFLAG       (O_RDWR|O_BINARY)
#define OFLAG       (O_RDWR|O_CREAT|O_BINARY)
#define OPMOD       (O_RDONLY|O_BINARY)

/********************************* MAT **************************************/

// v4 data types
// 0:	float*8
// 1:	float*4
// 2:	int*4
// 3:	int*2
// 4:	uint*4

// v5 data types
// 0:	???
// 1:	int*1
// 2:	uint*1
// 3:	int*2
// 4:	uint*2
// 5:	int*4
// 6:	uint*4
// 7:	float*4
// 8:	???
// 9:	float*8
//10:	???
//11:	???
//12:	int*8
//13:	uint*8
//14:	matrix
//15:	compressed
//16:	utf*1
//17:	utf*2
//18:	utf*4

static int bs4[6] = {8, 4, 4, 2, 4, 2};
static int bs5[20] = {
    0, 1, 1, 2, 2, 4, 4, 4, 0, 8,
    0, 0, 8, 8, 0, 0, 1, 2, 4, 0
};
static int dt4[20] = {
    -1, -1, -1,  3,  3,  2,  4,  1, -1,  0,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
};
static int dt5[6] = {9, 7, 5, 3, 6, 4};
static unsigned char *uncmp = NULL;

/*............................ MAT READ  ................................*/

typedef struct {
    int32_t data, next, nbyt, ofst, imag, rows, cols;
    char type, text, cmpx, vers, bswp;
    char name[MAXNCH];
} MATHDR;

/* mswab - swap bytes of m-byte sets */

static void
mswab(char *a, int m, int n)
{
    char *b, *c, t;

    while (n-- > 0) {
        b = a;
	c = a + m - 1;
	while (c > b) {
	    t = *b;
	    *b++ = *c;
	    *c-- = t;
	}
        a += m;
    }
}

/* cpu_lsb - check processor byte order */

static int
cpu_lsb(int i)
{
    int16_t   a = 1;
    char   *b = (char *) (&a);

    return (b[i]);
}

/* ismat5 - test if file is MAT v5 format */

static int
ismat5(char *s)
{
    int     status = 0;

    if (strncmp(s, "MATLAB", 6) == 0
	&& s[124] == 0
	&& s[125] == 1
	&& s[126] == 'I'
	&& s[127] == 'M') {
	status = 1;
    }

    return (status);
}

/* ismat4 - test if file is MAT v4 format */

static int
ismat4(char *s, int bs)
{
    char    typ[4];
    int     status = 0;
    uint32_t   *hdr = (uint32_t *)s;
    
    if (bs)
	mswab((char *) hdr, 4, 5);
    if (hdr[0] > 9999)
	return (0);
    typ[0] = (char) (hdr[0] / 1000) % 10;  /* M */
    typ[1] = (char) (hdr[0] / 100) % 10;   /* O */
    typ[2] = (char) (hdr[0] / 10) % 10;    /* P */
    typ[3] = (char) (hdr[0] / 1) % 10;     /* T */
    if (   (typ[0] == 0 || typ[0] == 1)
	&& (typ[1] == 0 || typ[1] == 1)
	&& (typ[2] >= 0 && typ[2] <= 4)
	&& (typ[3] == 0 || typ[3] == 1)
	&& (hdr[0] == 1 || hdr[1] > 0 || hdr[2] > 0)
	&& (hdr[3] == 0 || hdr[3] == 1)
	&& (hdr[4] > 0 && hdr[4] <= MAXNCH))
	status = 1;

    return (status);
}

/* rdhdr5 - read MAT v5 header */

static int
rdhdr5(FILE *fp, MATHDR *mh)
{
    int32_t nc, nbd, ofst, gap;
    uint32_t h[10];
    uint16_t w[2];

    if (fread(h, 4, 2, fp) < 2)
	return (0);
    mh->next += h[1] + 8;
    if (h[0] == 15) {	// compressed data ?
#ifndef ZLIB
	strcpy(mh->name, "");
        mh->type = h[0];
        mh->rows = h[1];
	mh->data = mh->next + 8;
        mh->vers = 6;
#else
#undef FAR
#define FAR
#include <zlib.h>
	Bytef *compr;
	uLongf nbc, nbu;
	int status;

        mh->data = ftell(fp);
	nbc = h[1];
	nbu = nbc * 10 + 100; // ???
	compr = (Bytef *) malloc(nbc);
	uncmp = (Bytef *) malloc(nbu);
        if (fread(compr, 1, nbc, fp) < (int) nbc) {
	    strcpy(mh->name, "read error");
	    return (0);
	}
	status = uncompress(uncmp, &nbu, compr, nbc);
	if (status == Z_OK) {
	    strcpy(mh->name, "");
	} else if (status == Z_MEM_ERROR) {
	    strcpy(mh->name, "mem error");
	} else if (status == Z_BUF_ERROR) {
	    strcpy(mh->name, "buf error");
	} else if (status == Z_DATA_ERROR) {
	    strcpy(mh->name, "data error");
	} else {
	    strcpy(mh->name, "<error>");
	}
	if (status != Z_OK) {
            mh->type = (char) h[0];
            mh->rows = nbc;
	    return (1);
	}
	ofst = 0;
	memcpy(h, uncmp + ofst, 40);
	ofst += 40;
	memcpy(w, uncmp + ofst, 4);
	ofst += 4;
        if (w[1]) {
	    nc = w[1];
        } else {
	    memcpy(&nc, uncmp + ofst, 4);
	    ofst += 4;
        }
        memcpy(mh->name, uncmp + ofst, nc);
	ofst += nc;
        mh->name[nc] = '\0';
	ofst += 8 - (ofst % 8);
	memcpy(w, uncmp + ofst, 4);
	ofst += 4;
        if (w[1]) {
	    nbd = w[1];
        } else {
	    memcpy(&nbd, uncmp + ofst, 4);
	    ofst += 4;
        }
	free(compr);
        mh->rows = h[8];
        mh->cols = h[9];
        mh->type = (char) w[0];
        mh->text = (h[4] == 4);
        mh->cmpx = (w[0] == 0);
	mh->nbyt = nbd;
	mh->ofst = ofst;
        mh->vers = 6;
#endif
    } else {
        if (fread(h + 2, 4, 8, fp) < 8)
	    return (0);
        if (fread(w, 2, 2, fp) < 2)
	    return (0);
        if (w[1]) {
	    nc = w[1];
        } else {
	    if (fread(&nc, 4, 1, fp) < 1)
	        return (0);
        }
        if (fread(mh->name, nc, 1, fp) < 1)
	    return (0);
        mh->name[nc] = '\0';
        ofst = ftell(fp);
        if (ofst % 8)
	    ofst = fseek(fp, 8 - (ofst % 8), 1);
        if (fread(w, 2, 2, fp) < 2)
	    return (0);
        if (w[1]) {
	    nbd = w[1];
        } else {
	    if (fread(&nbd, 4, 1, fp) < 1)
	        return (0);
        }
	gap = (w[1] & 4) ? 4 : 8;
        mh->rows = h[8];
        mh->cols = h[9];
        mh->type = (char) w[0];
        mh->text = (h[4] == 4);
        mh->cmpx = (h[4] >> 11) & 1;
        mh->data = ftell(fp);
	mh->nbyt = nbd;
	mh->imag = mh->cmpx ? mh->data + nbd + gap : 0;
    }

    return (1);
}

/* rdhdr4 - read MAT v4 header */

static int
rdhdr4(FILE *fp, MATHDR *mh)
{
    int P, T, nbyt, ndat, nchr = 0;
    uint32_t hdr[12];
    /* P= [0,1,2,3,4]=>[R*8,R*4,I*4,I*2,U*4] */

    if (fread(hdr, 4, 5, fp) < 5)
	return (0);
    if (mh->bswp)
	mswab((char *) hdr, 4, 5);
    /* M = (hdr[0] / 1000) % 10; */
    /* O = (hdr[0] / 100) % 10;  */
    P = (hdr[0] / 10) % 10;
    T = (hdr[0] / 1) % 10;
    mh->rows = hdr[1];
    mh->cols = hdr[2];
    mh->cmpx = (char) hdr[3];
    nchr = hdr[4];
    if (nchr <= 0 || nchr > MAXNCH)
        return (0);
    if (fread(mh->name, nchr, 1, fp) < 1)
        return (0);
    nbyt = bs4[P];
    ndat = mh->rows * mh->cols;
    mh->data = ftell(fp);
    if (mh->cmpx) {
        mh->imag = mh->data + ndat * nbyt;
        mh->next = mh->data + ndat * nbyt *2;
    } else {
        mh->imag = 0;
        mh->next = mh->data + ndat * nbyt;
    }
    mh->type = dt5[P];
    mh->text = T;
    mh->nbyt = ndat * nbyt;

    return (1);
}

/* rdhdr - read MAT variable header */

static int
rdhdr(FILE *fp, MATHDR *mh)
{
    free_null(uncmp);
    fseek(fp, mh->next, 0);
    if (mh->vers > 4)
	return (rdhdr5(fp, mh));
    return (rdhdr4(fp, mh));
}

/* ismat - check MAT version */

static int
ismat(char *s, MATHDR *mh)
{
    if (ismat5(s)) {
	mh->next = 128;
	mh->vers = 5;
	mh->bswp = ((int16_t *)s)[62] & 0xFF;
    } else if (ismat4(s, 0)) {
	mh->next = 0;
	mh->vers = 4;
	mh->bswp = 0;
    } else if (ismat4(s, 1)) {
	mh->next = 0;
	mh->vers = 4;
	mh->bswp = 1;
    } else {
	mh->next = 0;
	mh->vers = 0;
	mh->bswp = 0;
    }
 
    return (mh->vers);
}

/* mat_version - return MAT-file version */

static int
mat_version(const char *fn)
{
    char s[128];
    int nv = 0;
    FILE *fp;
    MATHDR mh;

    if ((fp = fopen(fn, "rb")) == NULL) {
	return (0);
    }
    fread(s, 1, 128, fp);
    if (ismat(s, &mh)) {
	for(;;) {
	    if (!rdhdr(fp, &mh)) {
		break;
	    }
	    nv++;
	}
    }
    fclose(fp);

    return (mh.vers);
}

/* mat_size - count number of variables in MAT file */

static int
mat_size(const char *fn)
{
    char s[128];
    int nv = 0;
    FILE *fp;
    MATHDR mh;

    if ((fp = fopen(fn, "rb")) == NULL) {
	return (0);
    }
    fread(s, 1, 128, fp);
    if (ismat(s, &mh)) {
	for(;;) {
	    if (!rdhdr(fp, &mh)) {
		break;
	    }
	    nv++;
	}
    }
    fclose(fp);

    return (nv);
}

/* mat_name - read variable names from MAT file */

static int
mat_name(const char *fn, VAR *vl, int nv)
{
    char s[128];
    int i;
    FILE *fp;
    MATHDR mh;

    if ((fp = fopen(fn, "rb")) == NULL) {
	return (0);
    }
    fread(s, 1, 128, fp);
    if (ismat(s, &mh)) {
	for(i = 0; i < nv; i++) {
	    if (!rdhdr(fp, &mh)) {
		break;
	    }
	    vl[i].name = _strdup(mh.name);
	    vl[i].rows = mh.rows;
	    vl[i].cols = mh.cols;
	    vl[i].dtyp = mh.type;
	    vl[i].text = mh.text;
	    vl[i].cmpx = mh.cmpx;
	    vl[i].data = NULL;
	}
    }
    fclose(fp);

    return (nv);
}

static int
rd_var(FILE *fp, MATHDR *mh, VAR *vl, int16_t *irc, int16_t *nrc)
{
    char *v, *w;
    int b, c, m, n, o, r, s, nr, nc;

    o = 0;
    n = mh->nbyt;
    c = mh->cmpx ? 2 : 1;
    s = mh->rows * mh->cols;
    if (s <= 0) {
        return (0);
    }
    b = mh->nbyt / s;
    if (irc) {		// specify initial row & column
	o = (irc[0] * mh->cols + irc[1]) * b;
	n = (n > o) ? (n - o) : 0;
    }
    if (nrc) {		// specify number of rows & columns
	nr = nrc[0];
	nc = nrc[1];
	m = nr * nc;
	if (s > m) {
	    s = m;
	    n = s * b;
	} else {
	    nr = 1;
	    nc = s;
	}
    }
    if (n <= 0) {
        return (0);
    }
    v = vl->data = malloc(n * c);
    if (nrc) {
	vl->rows = nr;
	vl->cols = nc;
    } else {
        vl->rows = mh->rows;
        vl->cols = mh->cols;
    }
    if (uncmp) {
        memcpy(v, uncmp + mh->ofst + o, n);
    } else {
        fseek(fp, mh->data + o, 0);
	r = fread(v, 1, n, fp);
	if (mh->imag) {
            fseek(fp, mh->imag + o, 0);
            r = fread(v + n, 1, n, fp);
	}
        fseek(fp, mh->next, 0);
        if (r < 0) {
            return (1);
        }
    }
    if (mh->cmpx) {
        w = malloc(n * c);
        n = b * s;
        for (r = 0; r < s; r++) {
            o = r * b;
	    memcpy(w + o * 2 + 0, v + o + 0, b);
	    memcpy(w + o * 2 + b, v + o + n, b);
        }
        vl->data = w;
        free(v);
    }
    if (mh->bswp) {
	mswab((char *) vl->data, b, s * c);
    }

    return (0);
}

/* mat_fetch - read specified variable from MAT file */

static int
mat_fetch(const char *fn, VAR *vl, char *vn, int16_t *irc, int16_t *nrc)
{
    char s[128];
    int i;
    FILE *fp;
    MATHDR mh;

    if ((fp = fopen(fn, "rb")) == NULL) {
	return (0);
    }
    fread(s, 1, 128, fp);
    if (ismat(s, &mh)) {
	for(;;) {
	    if (!rdhdr(fp, &mh)) {
		break;
	    }
	    if (strcmp(vn, mh.name) == 0) {
		break;
	    }
	}
        if (strcmp(vn, mh.name) == 0) {
	    i = 0;
	    rd_var(fp, &mh, vl + i, irc, nrc);
	    vl[i].name = _strdup(mh.name);
	    vl[i].dtyp = mh.type;
	    vl[i].text = mh.text;
	    vl[i].cmpx = mh.cmpx;
	    vl[i].last = 1;
	}
    }
    free_null(uncmp);
    fclose(fp);

    return (1);
}

/* mat_read - read variables from MAT file */

static int
mat_read(const char *fn, VAR *vl, int nv)
{
    char s[128];
    int i;
    FILE *fp;
    MATHDR mh;

    if ((fp = fopen(fn, "rb")) == NULL) {
	return (0);
    }
    fread(s, 1, 128, fp);
    if (ismat(s, &mh)) {
	for(i = 0; i < nv; i++) {
	    if (!rdhdr(fp, &mh)) {
		break;
	    }
	    vl[i].name = _strdup(mh.name);
	    vl[i].rows = mh.rows;
	    vl[i].cols = mh.cols;
	    vl[i].dtyp = mh.type;
	    vl[i].text = mh.text;
	    vl[i].cmpx = mh.cmpx;
	    if (rd_var(fp, &mh, vl + i, NULL, NULL)) {
	        break;
	    }	
	}
    }
    free_null(uncmp);
    fclose(fp);

    return (nv);
}

/*............................ MAT WRITE ................................*/

/*
 *  matlab header -
 *	type = MOPT where:
 *	M=0 for pc; 1 for sun
 *	O=0 col or 1 row;
 *	P= [0,1,2,3,4]=>[r*8,r*4,I*4,I*2,I*4 unsigned integers]
 *	T= 0 for matrix or 1 for text (stored as r*4 numbers 0<i<255).
 */
static int32_t
encode_mopt(int p, int t)
{
    int     a = 1;
    char   *b = (char *) (&a);

    return (b[1] * 1000 + p * 10 + t);    
}

static void
mat_wr(FILE *fp, char *nam, void *d, int r, int c, int dtyp, int txt, int cx)
{
    int     sz, sw, ne;
    int32_t    ns, hdr[5];
    static int bo = 1; // desired byte order: 0=BE, 1=LE 

    ns = strlen(nam) + 1;
    hdr[0] = encode_mopt(dtyp, txt);
    hdr[1] = r;			/* rows */
    hdr[2] = c;			/* cols */
    hdr[3] = cx;		/* 0 for real, 1 for complex */
    hdr[4] = ns;	        /* name length (including null byte) */
    sz = bs4[dtyp];
    ne = r * c;			/* number of elements = rows * cols */
    if (cx) {			/* complex ? */
        ne *= 2;
    }
    sw = cpu_lsb(bo);		/* swap bytes ? */
    if (sw) {
        mswab((char *) hdr, 4, 5);
        mswab((char *) d, sz, ne);
    }
    fwrite(hdr, 4, 5, fp);
    fwrite(nam, 1, ns, fp);
    fwrite(d, sz, ne, fp);
    if (sw) {
        mswab((char *) hdr, 4, 5);
        mswab((char *) d, sz, ne);
    }
}

/*..........................................................................*/

typedef struct {
    int lstsiz;
    VAR *varlst;
} VRL;

static int nvl = 0;
static int nvlsiz = 256;
static VRL *vrl = NULL;

static void
var_clear(int n)
{
    int i, nv;
    VAR *vl;

    vl = vrl[n].varlst;
    nv = vrl[n].lstsiz;
    for (i = 0; i < nv; i++) {
	free_null(vl[i].name);
	free_null(vl[i].data);
    }
    free_null(vl);
}

static int
var_size(VAR *vl)
{
    int i;

    for (i = 0; i < SP_MAXVAR; i++) {
	if (vl[i].last) {
	    break;
	}
    }

    return (i + 1);
}

static int
var_verify(VAR *vl)
{
    int i;
    VAR *v1, *v2;

    for (i = 0; i < nvl; i++) {
	v1 = vrl[i].varlst;
	v2 = vrl[i].lstsiz + v1;
	if (vl >= v1 && vl < v2) {
	    return (1);
	}
    }

    return (0);
}

static void
var_string(VAR *vl, const char *name, void *strptr)
{
    char *s;
    float *data;
    int i, n;

    s = (char *) strptr;
    n = strlen(s) + 1;
    data = (float *) calloc(n, sizeof(float));
    for (i = 0; i < n; i++) {
        data[i] = s[i];
    }
    sp_var_set(vl, name, data, 1, n, "f4");
    vl[0].text = 1;
    free(data);
}

/*..........................................................................*/

FUNC(VAR *) sp_mat_fetch(
    const char  *fn,
    char  *vn,
    int16_t *irc,
    int16_t *nrc
)
{
    VAR *vl;

    vl = sp_var_alloc(1);
    mat_fetch(fn, vl, vn, irc, nrc);

    return (vl);
}

FUNC(VAR *) sp_mat_load(
    const char *fn
)
{
    int nv;
    VAR *vl = NULL;

    nv = mat_size(fn);
    if (nv > 0) {
        vl = sp_var_alloc(nv);
        mat_read(fn, vl, nv);
    }

    return (vl);
}

FUNC(VAR *) sp_mat_whos(
    const char *fn
)
{
    int nv;
    VAR *vl;

    nv = mat_size(fn);
    vl = sp_var_alloc(nv);
    mat_name(fn, vl, nv);

    return (vl);
}

FUNC(int) sp_mat_append(
    const char *fn,
    VAR *vl
)
{
    int  i, dt;
    FILE *fp;

    fp = fopen(fn, "ab");		// open the file
    if (fp == NULL) {
	return (-1);
    }
    fseek(fp, 0L, SEEK_END);
    for (i = 0; i < SP_MAXVAR; i++) {
	dt = dt4[(int)vl[i].dtyp];
	if (vl[i].name && vl[i].data && dt >= 0
	    && vl[i].rows > 0 && vl[i].cols > 0) {
	    mat_wr(fp, vl[i].name, vl[i].data,
		vl[i].rows, vl[i].cols,
	        dt, vl[i].text, vl[i].cmpx);
	}
	if (vl[i].last) {
	    break;
	}
    }
    fclose(fp);

    return (0);
}

FUNC(int) sp_mat_save(
    const char *fn,
    VAR *vl
)
{
    int  i, dt;
    FILE *fp;

    fp = fopen(fn, "wb");		// open the file
    if (fp == NULL) {
	return (-1);
    }
    for (i = 0; i < SP_MAXVAR; i++) {
	dt = dt4[(int)vl[i].dtyp];
	if (vl[i].name && vl[i].data && dt >= 0
	    && vl[i].rows > 0 && vl[i].cols > 0) {
	    mat_wr(fp, vl[i].name, vl[i].data,
		vl[i].rows, vl[i].cols,
	        dt, vl[i].text, vl[i].cmpx);
	}
	if (vl[i].last) {
	    break;
	}
    }
    fclose(fp);

    return (0);
}

FUNC(int) sp_mat_version(
    const char *fn
)
{
    return (mat_version(fn));
}

FUNC(int) sp_mat_size(
    const char *fn
)
{
    return (mat_size(fn));
}

FUNC(VAR *) sp_var_alloc(
    int nvar
)
{
    VAR *vl;

    vl = (VAR *) calloc(nvar, sizeof(VAR));
    vl[nvar - 1].last = 1;
    if (!vrl) {
	vrl = (VRL *) calloc(nvlsiz, sizeof(VRL));
	nvl = 0;
    } else if (nvl == nvlsiz) {
	nvlsiz *= 2;
	vrl = (VRL *) realloc(vrl, nvlsiz * sizeof(VRL));
    }
    vrl[nvl].lstsiz = nvar;
    vrl[nvl].varlst = vl;
    nvl++;

    return (vl);
}

FUNC(void) sp_var_clear(
    VAR *vl
)
{
    int i, j;

    for (i = 0; i < nvl; i++) {
	if (vrl[i].varlst == vl) {
	    var_clear(i);
	    nvl--;
	    for (j = i; j < nvl; j++) {
		vrl[j].varlst = vrl[j+1].varlst;
	    }
	    vrl[j].varlst = NULL;
	    break;
	}
    }
}

FUNC(void) sp_var_clear_all(
)
{
    int i;

    for (i = 0; i < nvl; i++) {
	var_clear(i);
    }
    nvl = 0;
}

FUNC(VAR *) sp_var_copy(
    VAR *vl
)
{
    int i, nv, nd, nb;
    VAR *vc;

    nv = var_size(vl);
    vc = sp_var_alloc(nv);
    for(i = 0; i < nv; i++) {
	vc[i].name = _strdup(vl[i].name);
	vc[i].rows = vl[i].rows;
	vc[i].cols = vl[i].cols;
	vc[i].dtyp = vl[i].dtyp;
	vc[i].cmpx = vl[i].cmpx;
	vc[i].text = vl[i].text;
	vc[i].last = vl[i].last;
	nb = bs5[(int)vl[i].dtyp];
	nd = vl[i].rows * vl[i].cols;
	nd = vl[i].cmpx ? 2 * nd : nd;
	vc[i].data = calloc(nd, nb);
	memcpy(vc[i].data, vl[i].data, nd * nb);
    }
    return (vc);
}

FUNC(void) sp_var_float(
    VAR *vl
)
{
    char *i1;
    double *f8;
    float *f4;
    int i, j, n;
    int32_t *i4;
    int16_t *i2;
    unsigned char *u1;
    uint32_t *u4;
    uint16_t *u2;

    for(i = 0; i < SP_MAXVAR; i++) {
	n = vl[i].rows * vl[i].cols;
	n = vl[i].cmpx ? 2 * n : n;
	switch (vl[i].dtyp) {
	case 1:		// int*1
	    i1 = (char *) vl[i].data;
	    f4 = (float *) calloc(n, sizeof(float));
	    for (j = 0; j < n; j++) {
		f4[j] = (float) i1[j];
	    }
	    free(i1);
	    vl[i].data = f4;
	    break;
	case 2:		// uint*1
	    u1 = (unsigned char *) vl[i].data;
	    f4 = (float *) calloc(n, sizeof(float));
	    for (j = 0; j < n; j++) {
		f4[j] = (float) u1[j];
	    }
	    free(u1);
	    vl[i].data = f4;
	    break;
	case 3:		// int*2
	    i2 = (int16_t *) vl[i].data;
	    f4 = (float *) calloc(n, sizeof(float));
	    for (j = 0; j < n; j++) {
		f4[j] = (float) i2[j];
	    }
	    free(i2);
	    vl[i].data = f4;
	    break;
	case 4:		// uint*2
	    u2 = (uint16_t *) vl[i].data;
	    f4 = (float *) calloc(n, sizeof(float));
	    for (j = 0; j < n; j++) {
		f4[j] = (float) u2[j];
	    }
	    free(u2);
	    vl[i].data = f4;
	    break;
	case 5:		// int*4
	    i4 = (int32_t *) vl[i].data;
	    f4 = (float *) calloc(n, sizeof(float));
	    for (j = 0; j < n; j++) {
		f4[j] = (float) i4[j];
	    }
	    free(i4);
	    vl[i].data = f4;
	    break;
	case 6:		// uint*4
	    u4 = (uint32_t *) vl[i].data;
	    f4 = (float *) calloc(n, sizeof(float));
	    for (j = 0; j < n; j++) {
		f4[j] = (float) u4[j];
	    }
	    free(u4);
	    vl[i].data = f4;
	    break;
	case 0:		// text as float*8
	case 9:		// float*8
	    f8 = (double *) vl[i].data;
	    f4 = (float *) calloc(n, sizeof(float));
	    for (j = 0; j < n; j++) {
		f4[j] = (float) f8[j];
	    }
	    free(f8);
	    vl[i].data = f4;
	    break;
	}
	if (vl[i].last) {
	    break;
	}
    }
}

FUNC(void) sp_var_set(
    VAR *vl,
    const char *name,
    void *data,
    int32_t rows,
    int32_t cols,
    const char *frmt
)
{
    char *s;
    double *d, *d1, *d2;
    float *f1, *f2;
    int   i, ir, ii, jr, ji, t, c, ndat, nimg, nbyt, nsiz;

    if (!var_verify(vl)) {
	return;
    }
    if (strncmp(frmt, "f4str", 5) == 0) {
        var_string(vl, name, data);
	return;
    }
    t = strchr(frmt, 'T') || strchr(frmt, 't');
    c = strchr(frmt, 'C') || strchr(frmt, 'c');
    vl->name = _strdup(name);
    vl->rows = rows;
    vl->cols = cols;
    vl->text = (char) t;
    vl->cmpx = (char) c;
    nbyt = 0;
    if (vl->text) {
	vl->cmpx = 0;
	vl->dtyp = dt5[0];
	nbyt = bs4[0];
    } else if (strchr(frmt, 'F') || strchr(frmt, 'f')) {
	if (strchr(frmt, '8')) {
	    vl->dtyp = dt5[0];
	    nbyt = bs4[0];
	} else {
	    vl->dtyp = dt5[1];
	    nbyt = bs4[1];
	}
    } else if (strchr(frmt, 'I') || strchr(frmt, 'i')) {
	if (strchr(frmt, '4')) {
	    vl->dtyp = dt5[2];
	    nbyt = bs4[2];
	} else {
	    vl->dtyp = dt5[3];
	    nbyt = bs4[3];
	}
    } else if (strchr(frmt, 'U') || strchr(frmt, 'u')) {
	if (strchr(frmt, '1')) {
	    vl->dtyp = SP_DTYP_U1;
	    nbyt = rows * cols;
	} else if (strchr(frmt, '4')) {
	    vl->dtyp = dt5[4];
	    nbyt = bs4[4];
	} else {
	    vl->dtyp = dt5[5];
	    nbyt = bs4[5];
	}
    }
    if (nbyt == 0) {
	return;
    }
    ndat = rows * cols;
    nimg = vl->cmpx ? ndat : 0;
    nsiz = ndat + nimg;
    vl->data = calloc(nsiz, nbyt);
    if (vl->text) {
	s = (char *) data;
	d = (double *) vl->data;
	for (i = 0; i < ndat; i++) {
	    d[i] = (double) s[i];
	}
    } else if (vl->cmpx) {
        if (vl->dtyp == dt5[0]) {               // f8c
	    d1 = (double *) data;
	    d2 = (double *) vl->data;
	    for (i = 0; i < ndat; i++) {
                ir = 2 * i;
                ii = ir + 1;
                jr = i;
                ji = i + ndat;
	        d2[jr] = d1[ir];
	        d2[ji] = d1[ii];
	    }
        } else if (vl->dtyp == dt5[1]) {        // f4c
	    f1 = (float *) data;
	    f2 = (float *) vl->data;
	    for (i = 0; i < ndat; i++) {
                ir = 2 * i;
                ii = ir + 1;
                jr = i;
                ji = i + ndat;
	        f2[jr] = f1[ir];
	        f2[ji] = f1[ii];
	    }
        } else {
	    memcpy(vl->data, data, nsiz * nbyt);
        }
    } else {
	memcpy(vl->data, data, nsiz * nbyt);
    }
}

FUNC(int) sp_var_idx(
    VAR *vl 
)
{
    int i;

    for (i = 0; i < SP_MAXVAR; i++) {
	if (vl[i].name == NULL) {
            break;
        }
	if (vl[i].last) {
            i = -1; // variable list full
	    break;
	}
    }
    return (i);
}

FUNC(void) sp_var_add(
    VAR *vl,
    const char *name,
    void *data,
    int32_t rows,
    int32_t cols,
    const char *frmt
)
{
    int idx;
 
    idx = sp_var_idx(vl);
    sp_var_set(vl + idx, name, data, rows, cols, frmt);
}

FUNC(int) sp_var_size(	    // count variables in list
    VAR *vl 		    // pointer to variable
)
{
    return (var_size(vl));
}

/*..........................................................................*/

FUNC(int)
sp_var_find(VAR *v, const char *vn)
{
    int i, n;

    n = var_size(v);
    for (i = 0; i < n; i++) {
        if (strcmp(v[i].name, vn) == 0) {
            break;
        }
    }
    if (i == n) {
        i = -1;
    }

    return (i);
}

FUNC(float)
sp_var_f4(VAR *v, const char *vn)
{
    int i;

    i = sp_var_find(v, vn);
    return *((float *) v[i].data);
}

FUNC(double)
sp_var_f8(VAR *v, const char *vn)
{
    int i;

    i = sp_var_find(v, vn);
    return *((double *) v[i].data);
}

FUNC(int16_t)
sp_var_i2(VAR *v, const char *vn)
{
    int i;

    i = sp_var_find(v, vn);
    return *((int16_t *) v[i].data);
}

FUNC(int32_t)
sp_var_i4(VAR *v, const char *vn)
{
    int i;

    i = sp_var_find(v, vn);
    return *((int32_t *) v[i].data);
}

FUNC(char *)
sp_var_dattyp(int dt)
{
    static char *data_type[] = {
	"??", "I1", "U1", "I2", "U2", "I4", "U4", "F4",
	"??", "F8", "??", "??", "I8", "U8", "MX", "CM",
	"T1", "T2", "T4"
    };
    static int ndt = sizeof(data_type) / sizeof(data_type[0]); 

    if ((dt < 0) || (dt >= ndt)) {
	dt = 0;
    }
    return (data_type[dt]);
}

/*..........................................................................*/
