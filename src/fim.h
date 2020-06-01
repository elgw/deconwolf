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

#ifndef fim_h
#define fim_h

#define INLINED inline __attribute__((always_inline))

/* fim : operations on 3D floating point images 
 * all allocations are done with fftw3f_malloc (for alignment)
 *
 * functions ending with `_ref` are reference implementations to be compared 
 * tweaked or alternative versions.
 * */


float fim_min(const float * A, size_t N);
float fim_mean(const float * A, size_t N);
float fim_max(const float * A, size_t N);
float fim_sum(const float * restrict A, size_t N);
void fim_minus(float * restrict  A, 
    const float * restrict B, 
    const float * restrict C, 
    const size_t N);
  // A = B - C

void fim_invert(float * restrict A, const size_t N);

void fim_set_min_to_zero(float * I, size_t N);

int fim_maxAtOrigo(const float * restrict V, const int M, const int N, const int P);
  /* Check that the MAX of the fim is in the middle
   * returns 1 on success.
   * Returns 0 if any of the image dimensions are even
   */

void fim_stats(const float * A, size_t N);
// Print some info about A to stdout

void fim_flipall(float * restrict T, const float * restrict A, const int a1, const int a2, const int a3);
  /* 
   * MATLAB:
   * T = flip(flip(flip(A,1),2),3)*/


void fim_insert(float * restrict T, const int t1, const int t2, const int t3, 
    const float * restrict F, const int f1, const int f2, const int f3);
  /* Insert F [f1xf2xf3] into T [t1xt2xt3] in the "upper left" corner 
   * MATLAB:
   * T(1:size(F,1), 1:size(F,2), 1:sizes(F,3) = F;
   * */
void fim_insert_ref(float * T, int t1, int t2, int t3, 
    float * F, int f1, int f2, int f3);

float * fim_get_cuboid(float * restrict A, const int M, const int N, const int P,
    const int m0, const int m1, const int n0, const int n1, const int p0, const int p1);
/* MATLAB:
 * Y = A(m0:m1, n0:n1, p0:p1)
 */

float * fim_subregion(float * restrict A, const int M, const int N, const int P, const int m, const int n, const int p);
/* MATLAB:
 * Y = A(1:m, 1:n, 1:p);
 */

float * fim_subregion_ref(float * A, int M, int N, int P, int m, int n, int p);

void fim_normalize_sum1(float * psf, int M, int N, int P);


float * fim_copy(const float * restrict V, const size_t N);
  // Return a newly allocated copy of V

float * fim_zeros(const size_t N);
  // Allocate and return an array of N floats

float * fim_constant(const size_t N, const float value);
  // Allocate and return an array of N floats sets to a constant value

void fim_circshift(float * restrict A, 
    const int M, const int N, const int P, 
    const int sm, const int sn, const int sp);
  /* Shift the image A [MxNxP] by sm, sn, sp in each dimension */

float * fim_expand(const float * restrict in, 
    const int pM, const int pN, const int pP, 
    const int M, const int N, const int P);
  /* "expand an image" by making it larger 
   * pM, ... current size
   * M, Nm ... new size
   * */

float fim_mse(float * A, float * B, size_t N);
  /* mean( (A(:)-B(:)).^(1/2) )
   */

void shift_vector(float * restrict V, 
    const int S, 
    const int N,
    const int k);
  /* Circular shift of a vector of length N with stride S by step k */

void shift_vector_buf(float * restrict V, 
    const int S, 
    const int N,
    int k, float * restrict buffer);

void fim_mult_scalar(float * I, size_t N, float x);

void fim_ut(void);

#endif
