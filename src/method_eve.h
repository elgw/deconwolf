#ifndef __method_eve_h__
#define __method_eve_h__

#include "dw.h"

/* Exponential Vector Extrapolation (EVE) */
float biggs_alpha_eve(const afloat * restrict Xk,
                      const afloat * restrict Xkm1,
                      const afloat * restrict Ukm1,
                      const afloat * restrict Ukm2,
                      const size_t wMNP);

/* Exponential Vector Extrapolation (EVE) */
float * deconvolve_eve(afloat * restrict im, const int64_t M, const int64_t N, const int64_t P,
                       afloat * restrict psf, const int64_t pM, const int64_t pN, const int64_t pP,
                       dw_opts * s);


#endif
