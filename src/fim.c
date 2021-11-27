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

#include "fim.h"

//typedef float afloat __attribute__ ((__aligned__(16)));
typedef float afloat;

static float * gaussian_kernel(float sigma, size_t * nK);


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
    printf("max I(%" PRId64 ", %" PRId64 ", %" PRId64 ")=%f > mid I(%" PRId64 ", %" PRId64 ", %" PRId64 ")=%f\n",
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
     if(f3 > t3 || f2 > t2 || f1 > t1)
    {
        printf("Error: can't insert image of size %" PRId64
               " x %" PRIu64
               " x %" PRIu64
               " into target of size "
               "%" PRIu64
               " x %" PRIu64
               " x %" PRIu64
               "\n",
               f1, f2, f3, t1, t2, t3);
        exit(-1);
    }
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
    if(f3 > t3 || f2 > t2 || f1 > t1)
    {
        printf("Error: can't insert image of size %" PRId64
               " x %" PRIu64
               " x %" PRIu64
               " into target of size "
               "%" PRIu64
               " x %" PRIu64
               " x %" PRIu64
               "\n",
               f1, f2, f3, t1, t2, t3);
        exit(-1);
    }
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


    ((void) P);

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
        assert(Aidx < (size_t) M*N*P);
        // New coordinates are offset ...
        size_t Cidx = (aa - m0) +
          (bb - n0)*m +
          (cc - p0)*m*n;
        assert(Cidx < (size_t) m*n*p);
        C[Cidx] = A[Aidx];
      }
    }
  }
  return C;
}

afloat * fim_subregion(afloat * restrict A, const int64_t M, const int64_t N, const int64_t P, const int64_t m, const int64_t n, const int64_t p)
{
    ((void) P);

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
        assert(Aidx < (size_t) M*N*P);
        assert(Sidx < (size_t) m*n*p);
        S[Sidx] = A[Aidx];
      }
    }
  }
  return S;
}

afloat * fim_subregion_ref(afloat * A, int64_t M, int64_t N, int64_t P, int64_t m, int64_t n, int64_t p)
{
    ((void) P);
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
        assert(Aidx < (size_t) M*N*P);
        assert(Sidx < (size_t) m*n*p);
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
  for(size_t pp = 0; pp<(size_t) N; pp++)
  {
    buffer[pp] = V[pp*S];
  }
  for(size_t pp = 0; pp < (size_t) N; pp++)
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
  for(size_t kk = 0; kk < (size_t) M*N*P; kk++)
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
    printf("shift: %" PRId64 " -> ", k);
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
  size_t N = 0;
  float sigma = 1;
  float * K = gaussian_kernel(sigma, &N);
  assert(N>0);
  printf("gaussian_kernel, sigma=%f\n", sigma);
  for(size_t kk = 0; kk<N; kk++)
  {
      printf("%f ", K[kk]);
  }
  printf("\n");
  free(K);
  sigma = 0.5;
  K = gaussian_kernel(sigma, &N);
  printf("gaussian_kernel, sigma=%f\n", sigma);
  for(size_t kk = 0; kk<N; kk++)
  {
      printf("%f ", K[kk]);
  }
  printf("\n");
  free(K);

}

#if 0
static size_t min_size_t(size_t a, size_t b)
{
  if(a < b)
    return a;
  return b;
}
#endif

static size_t max_size_t(size_t a, size_t b)
{
  if(a > b)
    return a;
  return b;
}

void conv1(float * restrict V, int stride, float * restrict W,
    const size_t nV,
    const float * restrict K, const size_t nKu)
{
  const size_t k2 = (nKu-1)/2;
  const size_t N = nV;
  size_t bpos = 0;

  // First part
  for(size_t vv = 0;vv<k2; vv++)
  {
    double acc0 = 0;
    for(size_t kk = k2-vv; kk<nKu; kk++)
    {
        acc0 = acc0 + K[kk]*V[(vv-k2+kk)*stride];
    }
    W[bpos++] = acc0;
  }

  // Central part where K fits completely
  for(size_t vv = k2 ; vv+k2 < N; vv++)
  {
    double acc = 0;
    for(size_t kk = 0; kk<nKu; kk++)
    {
      acc = acc + K[kk]*V[(vv-k2+kk)*stride];
     }
    W[bpos++] = acc;
  }

  // Last part
  for(size_t vv = N-k2;vv<N; vv++)
  {
    double acc0 = 0;
 for(size_t kk = 0; kk<N-vv+k2; kk++)
    {
        acc0 = acc0 + K[kk]*V[(vv-k2+kk)*stride];
    }
    W[bpos++] = acc0;
  }


for(size_t pp = 0; pp<nV; pp++)
{
  V[pp*stride] = W[pp];
}
  return;
}

static float * gaussian_kernel(float sigma, size_t * nK)
{
    /* A Gaussian kernel */

    /* Determine the size so that most of the signal is captured */
    int len = 1; /* The total number of elements will be at least 3 */
    while(erf((len+1.0)/sigma) < 1.0-1e-8)
    {
        len++;
    }
    int N = 2*len + 1;

    float * K = malloc(N*sizeof(float));
    float mid = (N-1)/2;

    float s2 = pow(sigma, 2);
    for(int kk = 0; kk<N; kk++)
    {
        float x = (float) kk-mid;
        K[kk] = exp(-0.5*pow(x,2)/s2);
    }

/* Normalize the sum to 1 */
    float sum = 0;
    for(int kk = 0; kk<N; kk++)
        sum+=K[kk];

    for(int kk = 0; kk<N; kk++)
        K[kk]/=sum;

    nK[0] = N;
    return K;
}

void fim_gsmooth(float * restrict V, size_t M, size_t N, size_t P, float sigma)
{
    /* Convolve V by an isotropic Gaussian
     * implemented as a separated convolution. */

    if(sigma < 0)
    {
        printf("fim_gsmooth sigma=%f does not make sense.", sigma);
    }

  size_t nW = max_size_t(M, max_size_t(N, P));
  printf("fim_gsmooth: M: %zu, N: %zu, P: %zu, nW: %zu\n", M, N, P, nW); fflush(stdout);

  /* Temporary storage/buffer for conv1 */
float * W = malloc(nW*sizeof(float));

/* Create a kernel  */
size_t nK = 0;
float * K = gaussian_kernel(sigma, &nK);
assert(nK > 0);

// X
for(size_t pp = 0; pp < P; pp++)
{
for(size_t nn = 0; nn < N; nn++)
{
  conv1(V+pp*(M*N)+nn*M, 1, W, M, K, nK);
}
}

if(1){
// Y
for(size_t pp = 0; pp<P; pp++)
{
for(size_t mm = 0; mm<M; mm++)
{
  conv1(V + pp*(M*N) + mm, M, W, N, K, nK);
}
}
}

if(1){

// Z
for(size_t mm = 0; mm<M; mm++)
{
for(size_t nn = 0; nn<N; nn++)
{
  conv1(V+mm+M*nn, M*N, W, P, K, nK);
}
}
}

free(W);
return;
}
