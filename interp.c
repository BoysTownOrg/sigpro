// interp.c

#include "sigpro.h"

FUNC(int) sp_interp(
    float *x1, float *y1, int n1,
    float *x2, float *y2, int n2)
{
    float   a;
    int     i, j;

    if (n1 < 2)
	return (1);		// table too small
    for (j = 1; j < n1; j++) {
	if (x1[j - 1] > x1[j])
	    return (2);		// table nonmontonic
    }
    j = 0;
    for (i = 0; i < n2; i++) {
	if ((x2[i] < x1[0]) || (x2[i] > x1[n1 - 1])) {
	    return (3);		// outside table range
	}
	while ((j < (n2 - 1)) && (x2[i] >= x1[j + 1])) {
	    j++;
	}
	if (x2[i] == x1[j]) {
	    y2[i] = y1[j];
	} else if (x2[i] == x1[j + 1]) {
	    y2[i] = y1[j + 1];
	} else {
	    a = (x1[j + 1] - x2[i]) / (x1[j + 1] - x1[j]);
	    y2[i] = a * y1[j] + (1 - a) * y1[j + 1];
	}
    }

    return (0);
}

