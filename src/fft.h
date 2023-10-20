/*    Copyright (C) 2020 Erik L. G. Wernersson
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/* Provides some abstraction over fftw. Please note that this module
 * is not thread safe since it has some global data
 */

#pragma once

#include <fftw3.h>
#include <stdint.h>
#ifndef __APPLE__
#include <malloc.h>
#endif

#include "dw_util.h"

/*
 * Initialization commands.
 *
 */

/* @brief select FFTW planning type
 *
 * Set the plan to use, call before myfftw_start
 * can be FFTW_ESTIMATE, FFTW_MEASURE etc...
 * see the FFTW3 documentation.
 * TODO: Integrate with myfftw_start
*/
void fft_set_plan(unsigned int plan);

/* @brief To enable inplace-FFTs
 *
 * Set to use in-place transformations when possible
 * call before myfftw_start
 * TODO: integrate with myfftw_start
*/
void fft_set_inplace(int use_inplace);

/** @brief Required initialization routines
 *
 * Initialize, run before any other commands futher down.
 *log can be NULL */
void myfftw_start(int nThreads, int verbosity, FILE * log);

/** @brief Generate FFTW plans for the specified size
 *
 * Will generate both in-place and out-of place
 * has to be called before using the fft.
 */
void fft_train(size_t M, size_t N, size_t P,
               int verbosity, int nThreads,
               FILE * log);


/* @brief Free allocated memory
 *
 * Call this when you are done.
*/
void myfftw_stop(void);

void dim3_real_float_inverse(fftwf_complex * in, float * out,
                             const int n1, const int n2, const int n3);

void dim3_real_float(float * in, fftwf_complex* out,
                     const int n1, const int n2, const int n3);


/* Return the FFT of X. X is also freed (or re-used for in-place
 * transformations)
 */
fftwf_complex * fft_and_free(float * restrict X,
                             const int n1, const int n2, const int n3);

float * ifft_and_free(fftwf_complex * F,
                      const size_t n1, const size_t n2, const size_t n3);

/* Fast Fourier Transform, out of place */
fftwf_complex * fft(const float * in, int n1, int n2, int n3);

/* In-Place fft can be faster for small problems but is slower for
 * large ones The time it takes to pad the data can be neglected.
 * After return the input pointer should not be used and does not need
 * to be freed.
 */
fftwf_complex * fft_inplace(float * X,
                            const size_t M, const size_t N, const size_t P);


float * ifft(const fftwf_complex * fX,
             size_t M, size_t N, size_t P);
float * ifft_inplace(fftwf_complex * fX,
                     const size_t M, const size_t N, const size_t P);

void fft_mul(fftwf_complex * restrict C,
             fftwf_complex * restrict A,
             fftwf_complex * restrict B,
             const size_t n1, const size_t n2, const size_t n3);

/* Y = ifft(A*B) */
float *
fft_convolve_cc(fftwf_complex * A, fftwf_complex * B,
                int M, int N, int P);

/* Y = ifff(conj(A)*B)) */
float *
fft_convolve_cc_conj(fftwf_complex * A, fftwf_complex * B,
                     int M, int N, int P);

/* Highly specialised versions where the second argument is freed */
float * fft_convolve_cc_f2(fftwf_complex * A, fftwf_complex * B, int M, int N, int P);
float * fft_convolve_cc_conj_f2(fftwf_complex * A, fftwf_complex * B, int M, int N, int P);



void fft_ut(void);

/* Benchmark 1D ffts of size from, from+1, ... to
 * return time for each size
 */
double * fft_bench_1d(int64_t from, int64_t to, int niter);
