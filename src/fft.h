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

#ifndef fft_h
#define fft_h

#include <fftw3.h>
#include <stdint.h>

/* Call this before any other commands
 * log can be NULL */
void myfftw_start(int nThreads, int verbosity, FILE * log);
/* Call this when you are done. */
void myfftw_stop(void);

void dim3_real_float_inverse(fftwf_complex * in, float * out,
    const int n1, const int n2, const int n3);

void dim3_real_float(float * in, fftwf_complex* out,
    const int n1, const int n2, const int n3);

fftwf_complex * fft(float * in, int n1, int n2, int n3);

void fft_mul(fftwf_complex * restrict C,
    fftwf_complex * restrict A,
    fftwf_complex * restrict B,
    const size_t n1, const size_t n2, const size_t n3);

/* Y = ifft(A*B) */
float * fft_convolve_cc(fftwf_complex * A, fftwf_complex * B, int M, int N, int P);
/* Y = ifff(conj(A)*B)) */
float * fft_convolve_cc_conj(fftwf_complex * A, fftwf_complex * B, int M, int N, int P);

/* Highly specialised versions where the second argument is freed */
float * fft_convolve_cc_f2(fftwf_complex * A, fftwf_complex * B, int M, int N, int P);
float * fft_convolve_cc_conj_f2(fftwf_complex * A, fftwf_complex * B, int M, int N, int P);

/* Set the plan to use, should only be called once before fft_train */
void fft_set_plan(unsigned int plan);

/* Generate FFTW plans for the specified size */
void fft_train(size_t M, size_t N, size_t P,
               int verbosity, int nThreads,
               FILE * log);

void fft_ut(void);

/* Benchmark 1D ffts of size from, from+1, ... to
 * return time for each size
 */
double * fft_bench_1d(int64_t from, int64_t to, int niter);

#endif
