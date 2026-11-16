#include "dw_time.h"

void
dw_time_get(dw_time * t)
{
#ifdef WINDOWS
    timespec_get(t, TIME_UTC); // since C11
#else
    clock_gettime(CLOCK_REALTIME, t);
#endif
    return;
}

double
dw_time_delta(dw_time * end, dw_time * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}
