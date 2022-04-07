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

#ifndef _fim_h_
#define _fim_h_

#include <assert.h>
#include <inttypes.h>
#include <fftw3.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#include "fft.h"
typedef float afloat;
#include "fim_tiff.h"
#include "ftab.h"

#ifdef _OPENMP // turned on with -fopenmp
#include <omp.h>
#endif


#define INLINED inline __attribute__((always_inline))


/* fim : operations on 3D floating point images
 * all allocations are done with fftw3f_malloc (for alignment)
 *
 * functions ending with `_ref` are reference implementations to be compared
 * tweaked or alternative versions.
 * */

typedef struct{
    float * V;
    size_t M;
    size_t N;
    size_t P;
} fim_image_t;

/* Create a new object with a pointer to V */
fim_image_t * fim_image_from_array(const float * V, size_t M, size_t N, size_t P);

float fim_min(const float * A, size_t N);
float fim_mean(const float * A, size_t N);
float fim_max(const float * A, size_t N);
float fim_sum(const float * restrict A, size_t N);

/* Standard deviation, normalizing by (N-1) */
float fim_std(const float * V, size_t N);

float * fim_maxproj(const float * A, size_t M, size_t N, size_t P);

float * fim_sumproj(const float * A, size_t M, size_t N, size_t P);

/* Cumulative sum along dimension dim
 * Only supports 2D images
*/
void fim_cumsum2(float *, size_t M, size_t N, int dim);

/* A = B - C */
void fim_minus(float * restrict  A,
    const float * restrict B,
    const float * restrict C,
    const size_t N);

/* A = B / C */
void fim_div(float * restrict  A,
               const float * restrict B,
               const float * restrict C,
               const size_t N);

/* A[kk] += B[kk] */
void fim_add(float * restrict A,
             const float * restrict B,
             size_t N);

void fim_invert(float * restrict A, const size_t N);

void fim_set_min_to_zero(float * , size_t N);

int fim_maxAtOrigo(const float * restrict V, const int64_t M, const int64_t N, const int64_t P);
  /* Check that the MAX of the fim is in the middle
   * returns 1 on success.
   * Returns 0 if any of the image dimensions are even
   */

void fim_stats(const float * A, size_t N);
// Print some info about A to stdout

void fim_flipall(float * restrict T, const float * restrict A, const int64_t a1, const int64_t a2, const int64_t a3);
  /*
   * MATLAB:
   * T = flip(flip(flip(A,1),2),3)*/


void fim_insert(float * restrict T,
                const int64_t t1, const int64_t t2, const int64_t t3,
                const float * restrict F,
                const int64_t f1, const int64_t f2, const int64_t f3);
  /* Insert F [f1xf2xf3] into T [t1xt2xt3] in the "upper left" corner
   * MATLAB:
   * T(1:size(F,1), 1:size(F,2), 1:sizes(F,3) = F;
   * */

void fim_insert_ref(float * T, int64_t t1, int64_t t2, int64_t t3,
    float * F, int64_t f1, int64_t f2, int64_t f3);

float * fim_get_cuboid(float * restrict A, const int64_t M, const int64_t N, const int64_t P,
    const int64_t m0, const int64_t m1, const int64_t n0, const int64_t n1, const int64_t p0, const int64_t p1);
/* MATLAB:
 * Y = A(m0:m1, n0:n1, p0:p1)
 */

float * fim_subregion(const float * restrict A, const int64_t M, const int64_t N, const int64_t P, const int64_t m, const int64_t n, const int64_t p);
/* MATLAB:
 * Y = A(1:m, 1:n, 1:p);
 */

float * fim_subregion_ref(float * A, int64_t M, int64_t N, int64_t P, int64_t m, int64_t n, int64_t p);

/* Normalize an image to have the sum 1.0 */
void fim_normalize_sum1(float * psf, int64_t M, int64_t N, int64_t P);

/* Return a newly allocated copy of V */
float * fim_copy(const float * restrict V, const size_t N);


/* Allocate and return an array of N floats */
float * fim_zeros(const size_t N);


/* Allocate and return an array of N floats sets to a constant value */
float * fim_constant(const size_t N, const float value);


/* Shift the image A [MxNxP] by sm, sn, sp in each dimension */
void fim_circshift(float * restrict A,
    const int64_t M, const int64_t N, const int64_t P,
    const int64_t sm, const int64_t sn, const int64_t sp);

/* Shift the image A [MxNxP] by dm, dn, dp in each dimension,
 * What is outside of the image is interpreted as zero */
void fim_shift(float * restrict A,
                   const int64_t M, const int64_t N, const int64_t P,
                   const float dm, const float dn, const float dp);


float * fim_expand(const float * restrict in,
    const int64_t pM, const int64_t pN, const int64_t pP,
    const int64_t M, const int64_t N, const int64_t P);
  /* "expand an image" by making it larger
   * pM, ... current size
   * M, Nm ... new size
   * */

float fim_mse(float * A, float * B, size_t N);
  /* mean( (A(:)-B(:)).^(1/2) )
   */

void shift_vector(float * restrict V,
    const int64_t S,
    const int64_t N,
    const int64_t k);
  /* Circular shift of a vector of length N with stride S by step k */

void shift_vector_buf(float * restrict V,
    const int64_t S,
    const int64_t N,
    int64_t k, float * restrict buffer);

/* Shift vector by interpolation */
void shift_vector_float_buf(afloat * restrict V, // data
                            const int64_t S, // stride
                            const int64_t N, // elements
                            int n, // integer shift
                            afloat * restrict kernel, // centered kernel used for sub pixels shift
                            const int nkernel, // kernel size (odd!)
                            afloat * restrict buffer);

/* Multiply a float array of size N by x */
void fim_mult_scalar(float * fim, size_t N, float x);

void fim_ut(void);


/* Arg max: find coordinates of largest element
* If more than one maxima, return the first found
*/
void fim_argmax(const float * fim,
                size_t M, size_t N, size_t P,
                int64_t * _aM, int64_t *_aN, int64_t *_aP);




/*  */
float * fim_local_sum(const float * A, size_t M, size_t N, size_t pM, size_t pN);

/* Cumulative sum along dimension dim */
void fim_cumsum(float * A, const size_t M, const size_t N, const int dim);



/* Normalized cross correlation between T and A
 * See MATLAB's normxcorr2
 * This function requires that T and A have the same size
 * The returned image is of size [2xM-1, 2xN-1] and the values
 * should be in the range [0,1]
*/
float * fim_xcorr2(const float * T, const float * A,
                   const size_t M, const size_t N);

/* Image binarization with Otsu's method */
float * fim_otsu(float * Im, size_t M, size_t N);


typedef struct{
    double * C; /* Counts */
    float left; /* Left edge of first bin */
    float right; /* Right edge of last bin */
    size_t nbin; /* Number of bins */
} fim_histogram_t;

/* Return a histogram, the number of bins as well as the bin edges are
 * set automatically. The smallest end largest value should end up in
 * the middle of the first and last bin. */
fim_histogram_t * fim_histogram(const float * Im, size_t N);

/* Return a threshold that separates H at percentile p */
float fim_histogram_percentile(const fim_histogram_t * H, float p);

/* Return a global threshold by Otsu's method */
float fim_histogram_otsu(fim_histogram_t * H);

/* Replace count C[kk] = log(1+C[kk]) */
void fim_histogram_log(fim_histogram_t * H);

void fim_histogram_free(fim_histogram_t * H);

/* 2D connected components using 6-connectivity. TODO: see
 * wu2005optimizing for something smarter. */
int * fim_conncomp6(float * Im, size_t M, size_t N);


/* Find local maxima in I */
//ftab_t * fim_lmax(const float * I, size_t M, size_t N, size_t P);
ftab_t * fim_lmax(const float * Im, size_t M, size_t N, size_t P);
/* Sort with largest value first */
void ftab_sort(ftab_t * T, int col);

/* Spatial convolution */

/* Convolution of a single vector
 * In MATLAB that would be
 * Y = convn(V, K, 'same') / convn(ones(size(V)), K, 'same')
 * With normalized == 1 it would be
 * Y = convn(V, K, 'same') / convn(ones(size(V)), K, 'same')
 * That is only useful for gaussians
 */
void fim_conv1_vector(float * restrict V, int stride, float * restrict W,
                  const size_t nV,
                      const float * restrict K, const size_t nKu,
    const int normalized);

/* In-place convolution of a 3D volume, V, with a 1D filter, K along
 * dimension dim (0,1 or 2). The value of normalized is passed on to
 * fim_conv1_vector
 */
int fim_convn1(float * restrict V, size_t M, size_t N, size_t P,
               float * K, size_t nK,
               int dim, const int normalized);


/* Gaussian smoothing, normalized at edges */
void fim_gsmooth(float * restrict V, size_t M, size_t N, size_t P, float sigma);

/* Gaussian smoothing, normalized at edges, separate values for lateral and axial
 * filter */
void fim_gsmooth_aniso(float * restrict V,
                       size_t M, size_t N, size_t P,
                       float lsigma, float asigma);

/* Laplacian of Gaussian filter */
float * fim_LoG(const float * V, size_t M, size_t N, size_t P,
                float sigmaxy, float sigmaz);

/* Different implementation, using shiftdim */
float * fim_LoG_S(const float * V, size_t M, size_t N, size_t P,
                float sigmaxy, float sigmaz);


/* Extract a line centered at (x, y, z) with nPix pixels along dimension dim */
double * fim_get_line_double(fim_image_t * Im,
                          int x, int y, int z,
                          int dim, int nPix);

/* Similar to MATLABs shiftfim, [M,N,P] -> [N,P,M] */
fim_image_t * fim_shiftdim(fim_image_t *);

#endif /* _fim_h_ */
