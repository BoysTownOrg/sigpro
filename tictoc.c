// tictoc.c - stopwatch-timer functions

#include <time.h>
#include "sigpro.h"

static double tic_time;

/*..........................................................................*/

FUNC(double) sp_tic(void)
{
    return (tic_time / CLOCKS_PER_SEC);
}

FUNC(double) sp_toc(void)
{
    return ((clock() - tic_time) / CLOCKS_PER_SEC);
}

/*..........................................................................*/
