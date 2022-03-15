#ifndef __method_rl_h__
#define __method_rl_h__

#include "dw.h"

/* Richardson-Lucy */
float * deconvolve_rl(afloat * restrict im, const int64_t M, const int64_t N, const int64_t P,
                       afloat * restrict psf, const int64_t pM, const int64_t pN, const int64_t pP,
                       dw_opts * s);

#endif
