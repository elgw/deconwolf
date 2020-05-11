#include <math.h>
#include <fftw3.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fim.h"


void shift_vector_ut()
{
  int N = 5;
  int S = 1; // stride
  float * V = malloc(N*sizeof(float));

  for(int k = -7; k<7; k++)
  {
    for(int kk = 0; kk<N; kk++)
    {V[kk] = kk;}
    printf("shift: % d -> ", k);
    shift_vector(V,S,N,k);
    for(int kk =0; kk<N; kk++)
    { printf("%.0f ", V[kk]);}
    printf("\n");
  }
}


int fim_maxAtOrigo(const float * restrict V, const int M, const int N, const int P)
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

  int mM = (M-1)/2;
  int mN = (N-1)/2;
  int mP = (P-1)/2;
  float maxV = 0;
  for(size_t kk = 0; kk< M*N*P; kk++)
  {
    V[kk] > maxV ? maxV = V[kk] : 0;
  }

  if(maxV < V[mM + mN*M + mP*M*N])
  { return 0; }

  return 1;
}

void fim_stats(float * A, size_t N)
{
  float amax = 0;
  for(size_t kk = 0; kk<N; kk++)
  {
    if(A[kk] > amax)
      amax = A[kk];
  }
  printf("max: %f\n", amax);
}


void fim_flipall(float * restrict T, const float * restrict A, const int a1, const int a2, const int a3)
  /* Equivalent to T = flip(flip(flip(A,1),2),3) in matlab */
{
  for(int aa = 0; aa<a1; aa++){
    for(int bb = 0; bb<a2; bb++){
      for(int cc = 0; cc<a3; cc++){
        int idxTo = aa + bb*a1 + cc*a1*a2;
        int idxFrom = (a1-aa-1) + (a2-bb-1)*a1 + (a3-cc-1)*a1*a2;

        T[idxTo] = A[idxFrom];
      }
    }
  }
}


void fim_insert(float * restrict T, const int t1, const int t2, const int t3, 
    const float * restrict F, const int f1, const int f2, const int f3)
  /* Insert F [f1xf2xf3] into T [t1xt2xt3] in the "upper left" corner */
{
  for(int pp = 0; pp<f3; pp++)
  {
    for(int nn = 0; nn<f2; nn++)
    {
      for(int mm = 0; mm<f1; mm++)
      {
        float x = F[mm + nn*f1 + pp*f1*f2];
        T[mm + nn*t1 + pp*t1*t2] = x;
      }
    }
  }
  return;
}

void fim_insert_ref(float * T, int t1, int t2, int t3, 
    float * F, int f1, int f2, int f3)
  /* Insert F [f1xf2xf3] into T [t1xt2xt3] in the "upper left" corner */
{
  for(int mm = 0; mm<f1; mm++)
  {
    for(int nn = 0; nn<f2; nn++)
    {
      for(int pp = 0; pp<f3; pp++)
      {
        float x = F[mm + nn*f1 + pp*f1*f2];
        T[mm + nn*t1 + pp*t1*t2] = x;
      }
    }
  }
  return;
}


float * fim_get_cuboid(float * restrict A, const int M, const int N, const int P,
    const int m0, const int m1, const int n0, const int n1, const int p0, const int p1)
{
  /* Create a new array from V using [m0, m1]x[n0, n1]x[p0, p1] */
  int m = m1-m0+1;
  int n = n1-n0+1;
  int p = p1-p0+1;

  float * C = fftwf_malloc(m*n*p*sizeof(float));
  assert(C != NULL);

  for(int aa = m0; aa <= m1; aa++)
  {
    for(int bb = n0; bb <= n1; bb++)
    {
      for(int cc = p0; cc <= p1; cc++)
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

float * fim_subregion(float * restrict A, const int M, const int N, const int P, const int m, const int n, const int p)
{
  /* Extract sub region starting at (0,0,0) */
  float * S = fftwf_malloc(m*n*p*sizeof(float));
  assert(S != NULL);
  for(int pp = 0; pp<p; pp++)
  {
    for(int nn = 0; nn<n; nn++)
    {
      for(int mm = 0; mm<m; mm++)
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

float * fim_subregion_ref(float * A, int M, int N, int P, int m, int n, int p)
{
  float * S = fftwf_malloc(m*n*p*sizeof(float));
  assert(S != NULL);
  for(int mm = 0; mm<m; mm++)
  {
    for(int nn = 0; nn<n; nn++)
    {
      for(int pp = 0; pp<p; pp++)
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

void fim_normalize_max1(float * psf, int M, int N, int P)
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

float * fim_copy(const float * restrict V, const size_t N)
  // Return a newly allocated copy of V
{
  float * C = fftwf_malloc(N*sizeof(float));
  memcpy(C, V, N*sizeof(float));
  return C;
}


float * fim_zeros(const size_t N)
  // Allocate and return an array of N floats
{
  float * A = fftwf_malloc(N*sizeof(float));
  memset(A, 0, N*sizeof(float));
  return A;
}

float * fim_constant(const size_t N, const float value)
  // Allocate and return an array of N floats sets to a constant value
{
  float * A = fftwf_malloc(N*sizeof(float));
  for(size_t kk = 0; kk<N; kk++)
  {
    A[kk] = value;
  }
  return A;
}

void fim_circshift(float * restrict A, 
    const int M, const int N, const int P, 
    const int sm, const int sn, const int sp)
  /* Shift the image A [MxNxP] by sm, sn, sp in each dimension */
{

  const size_t bsize = fmax(fmax(M, N), P);
  float * restrict buf = malloc(bsize*sizeof(float));

  // Dimension 1
  for(int cc = 0; cc<P; cc++)
  {
    for(int bb = 0; bb<N; bb++)
    {    
      //shift_vector(A + bb*M + cc*M*N, 1, M, sm);
      shift_vector_buf(A + bb*M + cc*M*N, 1, M, sm, buf);
    }
  }

  // Dimension 2
  for(int cc = 0; cc<P; cc++)
  {
    for(int aa = 0; aa<M; aa++)
    {    
      //shift_vector(A + aa+cc*M*N, M, N, sn);
      shift_vector_buf(A + aa+cc*M*N, M, N, sn, buf);
    }
  }

  // Dimension 3
  for(int bb = 0; bb<N; bb++)
  {
    for(int aa = 0; aa<M; aa++)
    {  
      //shift_vector(A + aa+bb*M, M*N, P, sp);
      shift_vector_buf(A + aa+bb*M, M*N, P, sp, buf);
    }
  }

  free(buf);
  return;
}


INLINED static int mod_int(const int a, const int b)
{
  int r = a % b;
  return r < 0 ? r + b : r;
}

void shift_vector_buf(float * restrict V, 
    const int S, 
    const int N,
    int k, float * restrict buffer)
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

void shift_vector(float * restrict V, 
    const int S, 
    const int N,
    const int k)
  /* Circular shift of a vector of length N with stride S by step k */
{

  float * buffer = malloc(N*sizeof(float));
  shift_vector_buf(V, S, N, k, buffer);
  free(buffer);
  return;
}


float * fim_expand(const float * restrict in, 
    const int pM, const int pN, const int pP, 
    const int M, const int N, const int P)
  /* "expand an image" by making it larger 
   * pM, ... current size
   * M, Nm ... new size
   * */
{
  assert(pM<=M);
  assert(pN<=N);
  assert(pP<=P);

  float * out = fftwf_malloc(M*N*P*sizeof(float));
  assert(in != NULL);
  assert(out != NULL);
  for(size_t kk = 0; kk<M*N*P; kk++)
    out[kk] = 0;
  fim_insert(out, M, N, P, in, pM, pN, pP);
  return out;
}

void fim_flipall_ut()
{

  float * a = malloc(3*3*3*sizeof(float));
  float * b = malloc(3*3*3*sizeof(float));
  float * c = malloc(3*3*3*sizeof(float));

  for(int kk = 0; kk<27; kk++)
  {
    a[kk] = 0;
  }

  a[13] = 1;
  fim_flipall(b, a, 3, 3, 3);
  assert(b[13] == 1);

  for(int kk = 0; kk<27; kk++)
  {
    a[kk] = rand();
  }

  fim_flipall(b, a, 3, 3, 3);
  fim_flipall(c, b, 3, 3, 3);
  for(int kk = 0; kk<27; kk++)
  {
    assert(a[kk] == c[kk]);
  }

  fim_flipall(b, a, 4, 3, 2);
  fim_flipall(c, b, 4, 3, 2);
  for(int kk = 0; kk<24; kk++)
    assert(a[kk] == c[kk]);

  fim_flipall(b, a, 2, 3, 4);
  fim_flipall(c, b, 2, 3, 4);
  for(int kk = 0; kk<24; kk++)
    assert(a[kk] == c[kk]);

  free(a); free(b); free(c);
  return;
}


