#ifndef __lanczos_h__
#define __lanczos_h__

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

/*
 * Lanczos-3 interpolation of v at floating point indices
 * assumes symmetry around 0
 * assumes that x >= 0 and that x+3 < nV
 */
double lanczos3(const double * v, size_t nV, double x);
double lanczos5(const double * v, size_t nV, double x);


#endif
