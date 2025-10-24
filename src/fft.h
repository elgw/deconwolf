#pragma once

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

typedef struct {
    unsigned int FFTW3_PLANNING;
    int n_thread;
    int verbose;
    int use_inplace;
    fftwf_plan plan_r2c;
    fftwf_plan plan_c2r;;
    fftwf_plan plan_r2c_inplace;
    fftwf_plan plan_c2r_inplace;
    /* Real size */
    int64_t M;
    int64_t N;
    int64_t P;
    /* Complex size */
    int64_t cM;
    int64_t CN;
    int64_t cP;
    FILE * log; // Can be NULL
    // Will be set by fft_init
    char * wisdom_file;
} dw_fft;

/* @brief select FFTW planning type
 *
 * Set the plan to use, call before myfftw_start
 * can be FFTW_ESTIMATE, FFTW_MEASURE etc...
 * see the FFTW3 documentation.
 * TODO: Integrate with myfftw_start
 */
void fft_set_plan(dw_fft *, unsigned int plan);

/* @brief To enable inplace-FFTs
 *
 * Set to use in-place transformations when possible
 * call before myfftw_start
 * TODO: integrate with myfftw_start
 */
void fft_set_inplace(dw_fft *, int use_inplace);

/** @brief Required initialization routines
 *
 * Initialize, run before any other commands futher down.
 * log can be NULL
 * This function will also perform fftw planning
 * planner flags can be:
 * FFTW_MEASURE, FFTW_ESTIMATE, FFTW_PATIENT, FFTW_EXHAUSTIVE
 */
dw_fft * dw_fft_new(int nThreads, int verbosity, FILE * log,
                    int64_t M, int64_t N, int64_t P,
                    int planner_flags);


/* @brief Free allocated memory
 *
 * Call this when you are done.
 */
void dw_fft_destroy(dw_fft *);

void dim3_real_float_inverse(fftwf_complex * in, float * out,
                             const int n1, const int n2, const int n3);

void dim3_real_float(float * in, fftwf_complex* out,
                     const int n1, const int n2, const int n3);


/* Return the FFT of X. X is also freed (or re-used for in-place
 * transformations)
 */
fftwf_complex * fft_and_free(dw_fft *, float * restrict X,
                             const int n1, const int n2, const int n3);

float * ifft_and_free(dw_fft *, fftwf_complex * F,
                      const size_t n1, const size_t n2, const size_t n3);

/* Fast Fourier Transform, out of place */
fftwf_complex * fft(dw_fft *, const float * in, int n1, int n2, int n3);

/* In-Place fft can be faster for small problems but is slower for
 * large ones The time it takes to pad the data can be neglected.
 * After return the input pointer should not be used and does not need
 * to be freed.
 */
fftwf_complex * fft_inplace(dw_fft *, float * X,
                            const size_t M, const size_t N, const size_t P);


float * ifft(dw_fft *, const fftwf_complex * fX,
             size_t M, size_t N, size_t P);

float * ifft_inplace(dw_fft *, fftwf_complex * fX,
                     const size_t M, const size_t N, const size_t P);

void fft_mul(dw_fft *, fftwf_complex * restrict C,
             fftwf_complex * restrict A,
             fftwf_complex * restrict B,
             const size_t n1, const size_t n2, const size_t n3);

/* Y = ifft(A*B) */
float *
fft_convolve_cc(dw_fft *, fftwf_complex * A, fftwf_complex * B,
                int M, int N, int P);

/* Y = ifff(conj(A)*B)) */
float *
fft_convolve_cc_conj(dw_fft *, fftwf_complex * A, fftwf_complex * B,
                     int M, int N, int P);

/**
 * @brief like fft_convolve_cc but 2nd argument freed
 *
 * @param B is freed during the call and should be set to NULL afterwards
 */
float * fft_convolve_cc_f2(dw_fft *, fftwf_complex * A,
                           fftwf_complex * B,
                           int M, int N, int P);

/**
 * @brief like fft_convolve_cc_conj but 2nd argument freed
 *
 * @param B is freed during the call and should be set to NULL afterwards
 */
float * fft_convolve_cc_conj_f2(dw_fft *, fftwf_complex * A,
                                fftwf_complex * B,
                                int M, int N, int P);


/**
 * @brief run unit tests
 */
void fft_ut(void);

/* Benchmark 1D ffts of size from, from+1, ... to
 * return time for each size
 */
double * fft_bench_1d(int64_t from, int64_t to, int niter);
