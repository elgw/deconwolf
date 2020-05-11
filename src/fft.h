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

#include <fftw3.h>

void myfftw_start(int nThreads);
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

float * fft_convolve_cc(fftwf_complex * A, fftwf_complex * B, int M, int N, int P);

float * fft_convolve_cc_conj(fftwf_complex * A, fftwf_complex * B, int M, int N, int P);


void fft_train(size_t M, size_t N, size_t P , int, int nThreads);
