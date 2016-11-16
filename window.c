// window.c 

#include <math.h>
#include "sigpro.h"

/*
 * window type: 
 *	0-Rect,    1-Bartlet
 *      2-Hanning, 3=Hamming, 4=Blackman, 5=Nuttall
 */
FUNC(int) sp_window(
    float   *x,	    // output array
    int      m,	    // array size
    int      wt	    // window type
) 
{
    double  dp, ph, w, mc, *aa;
    int     i, k, c, nc;
    static double tpi = 2 * M_PI;
    static double Hanning[3] = {0.5, 0.5};
    static double Hamming[3] = {0.53836, 0.46164};
    static double Blackman[3] = {0.42, 0.5, 0.08};
    static double Nuttall[4] = {0.355768, 0.487396, 0.144232, 0.012604};

    if (wt < 0 || wt > 4)
	return (1);	    // invalid window type

    switch (wt) {

    case SP_WT_RECT:
        for (i = 0; i < m; i++) {
    	    x[i] = 1;
	}
	break;
    case SP_WT_BARTLET:
	c = (m + 1) / 2;
        mc = m / 2.0;
        for (i = 0; i < c; i++) {
    	    x[i] = (float) (i / mc);
	}
        for (i = c; i < m; i++) {
    	    x[i] = (float) (2 - i / mc);
	}
	break;
    case SP_WT_HANNING:
        dp = tpi / m;
        mc = m / 2.0;
	nc = 2;
	aa = Hanning;
        for (i = 0; i < m; i++) {
    	    ph = (i - mc) * dp;
    	    w = aa[0];
    	    for (k = 1; k < nc; k++)
    	        w = w + aa[k] * cos(ph * k);
    	    x[i] = (float) w;
        }
	break;
    case SP_WT_HAMMING:
        dp = tpi / m;
        mc = m / 2.0;
	nc = 2;
	aa = Hamming;
        for (i = 0; i < m; i++) {
    	    ph = (i - mc) * dp;
    	    w = aa[0];
    	    for (k = 1; k < nc; k++)
    	        w = w + aa[k] * cos(ph * k);
    	    x[i] = (float) w;
        }
	break;
    case SP_WT_BLACKMAN:
        dp = tpi / m;
        mc = m / 2.0;
	nc = 3;
	aa = Blackman;
        for (i = 0; i < m; i++) {
    	    ph = (i - mc) * dp;
    	    w = aa[0];
    	    for (k = 1; k < nc; k++)
    	        w = w + aa[k] * cos(ph * k);
    	    x[i] = (float) w;
        }
	break;
    case SP_WT_NUTTALL:
        dp = tpi / m;
        mc = m / 2.0;
	nc = 4;
	aa = Nuttall;
        for (i = 0; i < m; i++) {
    	    ph = (i - mc) * dp;
    	    w = aa[0];
    	    for (k = 1; k < nc; k++)
    	        w = w + aa[k] * cos(ph * k);
    	    x[i] = (float) w;
        }
	break;
    }

    return (0);
}
