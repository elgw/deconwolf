#ifndef __method_ave_h__
#define __method_ave_h__

#include "dw.h"

/* Additive Vector Extrapolation (AVE) */


float biggs_alpha(const afloat * restrict g,
                  const afloat * restrict gm,
                  const size_t wMNP, int mode);

float * deconvolve_ave(afloat * restrict im,
                       const int64_t M, const int64_t N, const int64_t P,
                       afloat * restrict psf,
                       const int64_t pM, const int64_t pN, const int64_t pP,
                       dw_opts * s);

float iter(
    afloat ** xp, // Output, f_(t+1)
    const float * restrict im, // Input image
    fftwf_complex * restrict cK, // fft(psf)
    afloat * restrict f, // Current guess
    afloat * restrict W, // Bertero Weights
    const int64_t wM, const int64_t wN, const int64_t wP, // expanded size
    const int64_t M, const int64_t N, const int64_t P, // input image size
    __attribute__((unused)) const dw_opts * s);

#endif
