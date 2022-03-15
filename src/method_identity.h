#ifndef __method_identity_h__
#define __method_identity_h__

#include "dw.h"

/* Exponential Vector Extrapolation (EVE) */
float * deconvolve_identity(afloat * restrict im, const int64_t M, const int64_t N, const int64_t P,
                       afloat * restrict psf, const int64_t pM, const int64_t pN, const int64_t pP,
                       dw_opts * s);


#endif
