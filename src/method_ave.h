#ifndef __method_ave_h__
#define __method_ave_h__

#include "dw.h"

/* Additive Vector Extrapolation (AVE) */

float * deconvolve_ave(afloat * restrict im,
                       const int64_t M, const int64_t N, const int64_t P,
                       afloat * restrict psf,
                       const int64_t pM, const int64_t pN, const int64_t pP,
                       dw_opts * s);

float alpha_ave(const afloat * restrict g,
                const afloat * restrict gm,
                const size_t wMNP, int mode);
#endif
