// chirp.c 

#include <math.h>
#include "sigpro.h"

FUNC(void) sp_chirp(	// frequency-sweep signal
    float *x, 		// output array
    int n		// complex array size
)
{
    double df, ph;
    int i;

    df = M_PI / n;
    ph = 0;
    for (i = 0; i < n; i++) {
	ph += i * df;
	x[i] = (float) sin(ph);
    }
}
