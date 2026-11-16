#pragma once
#include <time.h>

typedef struct timespec dw_time;

// Get the current time
void dw_time_get(dw_time *);

// Difference between two timepoints
double dw_time_delta(dw_time *, dw_time *);

// MATLAB like macros for timings
#define tictoc dw_time tictoc_start, tictoc_end;
#define tic dw_time_get(&tictoc_start);
#define toc(X) {                                                      \
        dw_time_get(&tictoc_end);                                     \
        printf(#X);                                                   \
        printf(" %f s\n", dw_time_delta(&tictoc_end, &tictoc_start)); \
        fflush(stdout);                                               \
        }
