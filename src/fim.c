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

#include <math.h>
#include <fftw3.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "fim.h"

typedef float afloat __attribute__ ((__aligned__(16)));


int fim_maxAtOrigo(const afloat * restrict V, const int64_t M, const int64_t N, const int64_t P)
  /* Check that the MAX of the fim is in the middle
   * returns 1 on success.
   * Returns 0 if any of the image dimensions are even
   */
{
  if( M % 2 == 0)
  { return 0; }

  if( N % 2 == 0)
  { return 0; }

  if (P % 2 == 0)
  { return 0; }

  int64_t mM = (M-1)/2;
  int64_t mN = (N-1)/2;
  int64_t mP = (P-1)/2;

  float maxV = 0;
  int64_t m = 0, n = 0, p = 0;

  for(int64_t pp = 0; pp < P; pp++) {
    for(int64_t nn = 0; nn < N; nn++) {
      for(int64_t mm = 0; mm < M; mm++) {
        size_t idx = mm + nn*M + pp*M*N;
        if(V[idx] > maxV)
        {
          maxV = V[idx];
          m = mm; n = nn; p = pp;
        }
      }
    }
  }

  float midValue = V[mM + mN*M + mP*M*N];


  if(maxV > midValue)
  { 
    printf("max I(%ld, %ld, %ld)=%f > mid I(%ld, %ld, %ld)=%f\n", 
        m, n, p, maxV, mM, mN, mP, midValue);
    return 0; 
  }

  return 1;
}

float fim_sum(const afloat * restrict A, size_t N)
{
  double sum = 0;
#pragma omp parallel for reduction(+:sum)
  for(size_t kk = 0; kk<N; kk++)
    sum+=(double) A[kk];
  return (float) sum;
}

float fim_mean(const afloat * A, size_t N)
{
  return fim_sum(A, N)/(float) N;
}

float fim_min(const afloat * A, size_t N)
{
  float amin = INFINITY;
  for(size_t kk = 0; kk<N; kk++)
  {
    if(A[kk] < amin)
      amin = A[kk];
  }
  return amin;
}

void fim_minus(afloat * restrict  A, 
    const afloat * restrict B, 
    const afloat * restrict C, 
    const size_t N)
  // A = B - C
{
  size_t kk = 0;

#pragma omp parallel for
  for(kk = 0; kk<N; kk++)
  {
    A[kk] = B[kk] - C[kk];
  }
  return;
}

float fim_max(const afloat * A, size_t N)
{
  float amax = -INFINITY;
  for(size_t kk = 0; kk<N; kk++)
  {
    if(A[kk] > amax)
      amax = A[kk];
  }
  return amax;
}


void fim_stats(const afloat * A, const size_t N)
{
  printf("min: %f mean: %f, max: %f\n",
      fim_min(A, N),
      fim_mean(A, N),
      fim_max(A, N));
  return;
}

float fim_mse(afloat * A, afloat * B, size_t N)
  /* mean( (A(:)-B(:)).^(1/2) )
  */
{
  double mse = 0;
  for(size_t kk = 0; kk<N; kk++)
  {
    mse += pow(A[kk]-B[kk], 2);
  }
  return mse/N;
}

void fim_flipall(afloat * restrict T, const afloat * restrict A, const int64_t a1, const int64_t a2, const int64_t a3)
  /* Equivalent to T = flip(flip(flip(A,1),2),3) in matlab */
{
  for(int64_t aa = 0; aa<a1; aa++){
    for(int64_t bb = 0; bb<a2; bb++){
      for(int64_t cc = 0; cc<a3; cc++){
        int64_t idxTo = aa + bb*a1 + cc*a1*a2;
        int64_t idxFrom = (a1-aa-1) + (a2-bb-1)*a1 + (a3-cc-1)*a1*a2;

        T[idxTo] = A[idxFrom];
      }
    }
  }
}


void fim_insert(afloat * restrict T, const int64_t t1, const int64_t t2, const int64_t t3, 
    const afloat * restrict F, const int64_t f1, const int64_t f2, const int64_t f3)
  /* Insert F [f1xf2xf3] into T [t1xt2xt3] in the "upper left" corner */
{
  for(int64_t pp = 0; pp<f3; pp++)
  {
    for(int64_t nn = 0; nn<f2; nn++)
    {
      for(int64_t mm = 0; mm<f1; mm++)
      {
        T[mm + nn*t1 + pp*t1*t2] = F[mm + nn*f1 + pp*f1*f2];
      }
    }
  }
  return;
}

void fim_insert_ref(afloat * T, int64_t t1, int64_t t2, int64_t t3, 
    afloat * F, int64_t f1, int64_t f2, int64_t f3)
  /* Insert F [f1xf2xf3] into T [t1xt2xt3] in the "upper left" corner */
{
  for(int64_t mm = 0; mm<f1; mm++)
  {
    for(int64_t nn = 0; nn<f2; nn++)
    {
      for(int64_t pp = 0; pp<f3; pp++)
      {
        afloat x = F[mm + nn*f1 + pp*f1*f2];
        T[mm + nn*t1 + pp*t1*t2] = x;
      }
    }
  }
  return;
}


afloat * fim_get_cuboid(afloat * restrict A, const int64_t M, const int64_t N, const int64_t P,
    const int64_t m0, const int64_t m1, const int64_t n0, const int64_t n1, const int64_t p0, const int64_t p1)
{
  /* Create a new array from V using [m0, m1]x[n0, n1]x[p0, p1] */
  int64_t m = m1-m0+1;
  int64_t n = n1-n0+1;
  int64_t p = p1-p0+1;

  afloat * C = fftwf_malloc(m*n*p*sizeof(float));
  assert(C != NULL);

  for(int64_t aa = m0; aa <= m1; aa++)
  {
    for(int64_t bb = n0; bb <= n1; bb++)
    {
      for(int64_t cc = p0; cc <= p1; cc++)
      {
        // printf("aa:%d, bb:%d, cc:%d\n", aa, bb, cc);
        size_t Aidx = aa + bb*M + cc*M*N;
        assert(Aidx < M*N*P);
        // New coordinates are offset ...
        size_t Cidx = (aa - m0) + 
          (bb - n0)*m + 
          (cc - p0)*m*n;
        assert(Cidx < m*n*p);
        C[Cidx] = A[Aidx];
      }
    }
  }
  return C;
}

afloat * fim_subregion(afloat * restrict A, const int64_t M, const int64_t N, const int64_t P, const int64_t m, const int64_t n, const int64_t p)
{
  /* Extract sub region starting at (0,0,0) */
  afloat * S = fftwf_malloc(m*n*p*sizeof(float));
  assert(S != NULL);
  for(int64_t pp = 0; pp<p; pp++)
  {
    for(int64_t nn = 0; nn<n; nn++)
    {
      for(int64_t mm = 0; mm<m; mm++)
      {
        size_t Aidx = mm + nn*M + pp*M*N;
        size_t Sidx = mm + nn*m + pp*m*n;
        assert(Aidx < M*N*P);
        assert(Sidx < m*n*p);
        S[Sidx] = A[Aidx];
      }
    }
  }
  return S;
}

afloat * fim_subregion_ref(afloat * A, int64_t M, int64_t N, int64_t P, int64_t m, int64_t n, int64_t p)
{
  afloat * S = fftwf_malloc(m*n*p*sizeof(float));
  assert(S != NULL);
  for(int64_t mm = 0; mm<m; mm++)
  {
    for(int64_t nn = 0; nn<n; nn++)
    {
      for(int64_t pp = 0; pp<p; pp++)
      {
        size_t Aidx = mm + nn*M + pp*M*N;
        size_t Sidx = mm + nn*m + pp*m*n;
        assert(Aidx < M*N*P);
        assert(Sidx < m*n*p);
        S[Sidx] = A[Aidx];
      }
    }
  }
  return S;
}

void fim_set_min_to_zero(afloat * I, size_t N)
{
  float min = fim_min(I, N);
  for(size_t kk = 0; kk<N; kk++)
  {
    I[kk] -= min;
  }
}

void fim_mult_scalar(afloat * I, size_t N, float x)
{
  for(size_t kk = 0; kk < N ; kk++)
  {
    I[kk]*=x;
  }
}

void fim_normalize_sum1(afloat * psf, int64_t M, int64_t N, int64_t P)
  /* 
   * MATLAB:
   * Y = X/max(X(:))
   */
{
  size_t pMNP = M*N*P;;
  double psf_sum = 0;
  for(size_t kk = 0; kk<pMNP; kk++)
  { psf_sum += psf[kk]; }
  //  printf("psf_sum: %f\n", psf_sum);
  for(size_t kk = 0; kk<pMNP; kk++) 
  { psf[kk]/=psf_sum; }
}

afloat * fim_copy(const afloat * restrict V, const size_t N)
  // Return a newly allocated copy of V
{
  afloat * C = fftwf_malloc(N*sizeof(float));
  memcpy(C, V, N*sizeof(float));
  return C;
}

afloat * fim_zeros(const size_t N)
  // Allocate and return an array of N floats
{
  afloat * A = fftwf_malloc(N*sizeof(float));
  memset(A, 0, N*sizeof(float));
  return A;
}

afloat * fim_constant(const size_t N, const float value)
  // Allocate and return an array of N floats sets to a constant value
{
  afloat * A = fftwf_malloc(N*sizeof(float));
  for(size_t kk = 0; kk<N; kk++)
  {
    A[kk] = value;
  }
  return A;
}

void fim_circshift(afloat * restrict A, 
    const int64_t M, const int64_t N, const int64_t P, 
    const int64_t sm, const int64_t sn, const int64_t sp)
  /* Shift the image A [MxNxP] by sm, sn, sp in each dimension */
{

  const size_t bsize = fmax(fmax(M, N), P);
#pragma omp parallel
  {
    afloat * restrict buf = malloc(bsize*sizeof(float));

    // Dimension 1
#pragma omp for
    for(int64_t cc = 0; cc<P; cc++)
    {
      for(int64_t bb = 0; bb<N; bb++)
      {    
        //shift_vector(A + bb*M + cc*M*N, 1, M, sm);
        shift_vector_buf(A + bb*M + cc*M*N, 1, M, sm, buf);
      }
    }

    // Dimension 2
#pragma omp for
    for(int64_t cc = 0; cc<P; cc++)
    {
      for(int64_t aa = 0; aa<M; aa++)
      {    
        //shift_vector(A + aa+cc*M*N, M, N, sn);
        shift_vector_buf(A + aa+cc*M*N, M, N, sn, buf);
      }
    }

    // Dimension 3
#pragma omp for
    for(int64_t bb = 0; bb<N; bb++)
    {
      for(int64_t aa = 0; aa<M; aa++)
      {  
        //shift_vector(A + aa+bb*M, M*N, P, sp);
        shift_vector_buf(A + aa+bb*M, M*N, P, sp, buf);
      }
    }

    free(buf);
  }
  return;
}


INLINED static int64_t mod_int(const int64_t a, const int64_t b)
{
  int64_t r = a % b;
  return r < 0 ? r + b : r;
}

void shift_vector_buf(afloat * restrict V, 
    const int64_t S, 
    const int64_t N,
    int64_t k, afloat * restrict buffer)
  /* Circular shift of a vector of length N with stride S by step k */
{

  k = -k;
  for(size_t pp = 0; pp<N; pp++)
  {
    buffer[pp] = V[pp*S];
  }
  for(size_t pp = 0; pp<N; pp++)
  {
    V[pp*S] = buffer[mod_int(pp+k, N)];
  }
  return;
}

void shift_vector(afloat * restrict V, 
    const int64_t S, 
    const int64_t N,
    const int64_t k)
  /* Circular shift of a vector of length N with stride S by step k */
{

  afloat * buffer = malloc(N*sizeof(float));
  shift_vector_buf(V, S, N, k, buffer);
  free(buffer);
  return;
}


afloat * fim_expand(const afloat * restrict in, 
    const int64_t pM, const int64_t pN, const int64_t pP, 
    const int64_t M, const int64_t N, const int64_t P)
  /* "expand an image" by making it larger 
   * pM, ... current size
   * M, Nm ... new size
   * */
{
  assert(pM<=M);
  assert(pN<=N);
  assert(pP<=P);

  afloat * out = fftwf_malloc(M*N*P*sizeof(float));
  assert(in != NULL);
  assert(out != NULL);
  for(size_t kk = 0; kk<M*N*P; kk++)
    out[kk] = 0;
  fim_insert(out, M, N, P, in, pM, pN, pP);
  return out;
}

void fim_flipall_ut()
{

  float * a = fftwf_malloc(3*3*3*sizeof(float));
  float * b = fftwf_malloc(3*3*3*sizeof(float));
  float * c = fftwf_malloc(3*3*3*sizeof(float));

  for(int64_t kk = 0; kk<27; kk++)
  {
    a[kk] = 0;
  }

  a[13] = 1;
  fim_flipall(b, a, 3, 3, 3);
  assert(b[13] == 1);

  for(int64_t kk = 0; kk<27; kk++)
  {
    a[kk] = rand();
  }

  fim_flipall(b, a, 3, 3, 3);
  fim_flipall(c, b, 3, 3, 3);
  for(int64_t kk = 0; kk<27; kk++)
  {
    assert(a[kk] == c[kk]);
  }

  fim_flipall(b, a, 4, 3, 2);
  fim_flipall(c, b, 4, 3, 2);
  for(int64_t kk = 0; kk<24; kk++)
    assert(a[kk] == c[kk]);

  fim_flipall(b, a, 2, 3, 4);
  fim_flipall(c, b, 2, 3, 4);
  for(int64_t kk = 0; kk<24; kk++)
    assert(a[kk] == c[kk]);

  free(a); free(b); free(c);
  return;
}

void shift_vector_ut()
{
  int64_t N = 5;
  int64_t S = 1; // stride
  float * V = fftwf_malloc(N*sizeof(float));

  for(int64_t k = -7; k<7; k++)
  {
    for(int64_t kk = 0; kk<N; kk++)
    {V[kk] = kk;}
    printf("shift: %ld -> ", k);
    shift_vector(V,S,N,k);
    for(int64_t kk =0; kk<N; kk++)
    { printf("%.0f ", V[kk]);}
    printf("\n");
  }
}

void fim_invert(float * restrict A, size_t N)
{
  for(size_t kk = 0; kk < N; kk++)
  {
    A[kk] = 1-A[kk];
  }
return;
}

void fim_ut()
{
  fim_flipall_ut();
  shift_vector_ut();

}
