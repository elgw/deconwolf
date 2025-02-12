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

#include <assert.h>
#include <inttypes.h>
#include <fftw3.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>


#ifdef __linux__
#include <sys/mman.h>
#endif

#include "fim.h"
#include "quickselect.h"

typedef uint64_t u64;
typedef int64_t i64;

static int fim_verbose = 0;

static float * gaussian_kernel(float sigma, size_t * nK);
static void cumsum_array(float * A, size_t N, size_t stride);
static void fim_show(float * A, size_t M, size_t N, size_t P);
static void fim_show_int(int * A, size_t M, size_t N, size_t P);

static const char * fim_boundary_condition_str(fim_boundary_condition bc)
{
    switch(bc)
    {
    case FIM_BC_VALID:
        return "VALID";
    case FIM_BC_SYMMETRIC_MIRROR:
        return "SYMMETRIC MIRROR";
    case FIM_BC_PERIODIC:
        return "PERIODIC";
    case FIM_BC_WEIGHTED:
        return "WEIGHTED";
    case FIM_BC_ZEROS:
        return "ZEROS";
    }
    assert(0);
    return "UNKNOWN";
}

#ifdef __linux__
void * __attribute__((__aligned__(FIM_ALIGNMENT))) fim_malloc(size_t nbytes)
{
    void * p;
    if(posix_memalign(&p, FIM_ALIGNMENT, nbytes))
    {
        fprintf(stderr, "fim_malloc: unable to allocate %zu bytes\n", nbytes);
        assert(0);
        exit(EXIT_FAILURE);
    }

#if 0
    const size_t HPAGE_SIZE  = (1 << 21); // 2 Mb
    if(FIM_ALIGNMENT == HPAGE_SIZE)
    {
        /* Has to be done before writing the first byte */
        if(madvise( p, nbytes, MADV_HUGEPAGE ))
        {
            fprintf(stderr, "madvise failed\n");
        }
        if(0){ /* Would this decrease the chances of fragmentation? */
            /* Set one byte to allocate */
            /* Set one byte of each page to allocate */
            for (size_t kk = 0; kk < nbytes; kk += HPAGE_SIZE) {
                memset(p + kk, 0, 1);
            }
        }
    }
#endif

    memset(p, 0, nbytes);

    return p;
}
#else
void * fim_malloc(size_t nbytes)
{
#ifdef WINDOWS
    void * p = _aligned_malloc(nbytes, FIM_ALIGNMENT);
    memset(p, 0, nbytes);
    return p;
#else
    void * p;
    if(posix_memalign(&p, FIM_ALIGNMENT, nbytes))
    {
        fprintf(stderr, "Unable to allocate %zu bytes\n", nbytes);
        assert(0);
        exit(EXIT_FAILURE);
    }
    assert(p != NULL);
    return p;
#endif
}
#endif

#ifdef WINDOWS
void * fim_realloc(void * p, size_t nbytes)
{
    void * ret = _aligned_realloc(p, nbytes, FIM_ALIGNMENT);
    return ret;
}
#else
void * __attribute__((__aligned__(FIM_ALIGNMENT))) fim_realloc(void * p, size_t nbytes)
{
    void * out = realloc(p, nbytes);
    /* If address didn't change, we are good */
    if(p == out)
    {
        return out;
    }

    if( (uintptr_t)(const void *)(out) % FIM_ALIGNMENT != 0 )
    {
        void * out2 = fim_malloc(nbytes);
        memcpy(out2, out, nbytes);
        fim_free(out);
        return out2;
    }

    return out;
}
#endif

void fim_free(void * p)
{
#ifdef WINDOWS
    _aligned_free(p);
#else
    free(p);
#endif
}

void fim_set_verbose(int v)
{
    fim_verbose = v;
}

static double clockdiff(struct timespec* end, struct timespec * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}

/* https://en.wikipedia.org/wiki/Anscombe_transform  */
void fim_anscombe(float * x, size_t n)
{
#pragma omp parallel for
    for(size_t kk = 0; kk<n; kk++)
    {
        x[kk] = 2.0*sqrt(x[kk] + 3.0/8.0);
    }
}

void fim_ianscombe(float * x, size_t n)
{
#pragma omp parallel for
    for(size_t kk = 0; kk<n; kk++)
    {
        x[kk] = pow(x[kk] / 2.0, 2) - 3.0/8.0;
    }
}


void fimo_free(fimo * F)
{
    if(F != NULL)
    {
        fim_free(F->V);
        F->V = NULL;
        free(F);
    }
    return;
}


fimo * fim_wrap_array(float * V, size_t M, size_t N, size_t P)
{
    fimo * I = malloc(sizeof(fimo));
    assert(I != NULL);
    I->V = V;
    I->M = M;
    I->N = N;
    I->P = P;
    return I;
}

fimo * fim_image_from_array(const float * restrict V, size_t M, size_t N, size_t P)
{
    fimo * I = malloc(sizeof(fimo));
    assert(I != NULL);
    I->V = fim_malloc(M*N*P*sizeof(float));

    I->V = memcpy(I->V, V, M*N*P*sizeof(float));
    I->M = M;
    I->N = N;
    I->P = P;
    return I;
}

int fim_maxAtOrigo(const float * restrict V, const int64_t M, const int64_t N, const int64_t P)
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

    float maxV = fim_max(V, M*N*P);
    float midValue = V[mM + mN*M + mP*M*N];

    if(maxV > midValue)
    {
        printf("max I = %f > mid I(%" PRId64 ", %" PRId64 ", %" PRId64 ")=%f\n",
               maxV, mM, mN, mP, midValue);
        return 0;
    }

    return 1;
}

void fim_argmax(const float * I,
                const size_t M, const size_t N, const size_t P,
                int64_t * _aM, int64_t *_aN, int64_t *_aP)
{
    float v = 0;
    fim_argmax_max(I, M, N, P, _aM, _aN, _aP, &v);
    return;
}


void fim_argmax_max(const float * I,
                    size_t M, size_t N, size_t P,
                    int64_t * _aM, int64_t *_aN, int64_t *_aP,
                    float * max_value)
{
    float max = I[0];
    size_t amax = 0;

    size_t MNP = M*N*P;

#pragma omp parallel
    {
        size_t n = MNP/omp_get_num_threads();   // step = MAX/number of threads.
        size_t id = omp_get_thread_num();       // id is one of 0, 1, ..., (num_threads-1)
        size_t start = id * n;
        size_t stop = start;
        if ( id + 1 == (size_t) omp_get_num_threads() ) {
            stop = MNP;
        }
        else {
            stop = start + n;
        }

        size_t local_amax = start;
        float local_max = I[local_amax];

        for(size_t kk = start; kk<stop; kk++)
        {
            if(I[kk] > local_max)
            {
                local_max = I[kk];
                local_amax = kk;
            }
        }
#pragma omp critical
        {
            // printf("local_max = %f, local_amax = %zu, th=%zu [%zu, %zu]\n",
            //       local_max, local_amax, id, start, stop-1);
            if(local_max == max)
            {
                if(local_amax < amax)
                {
                    max = local_max;
                    amax = local_amax;
                }
            }
            if(local_max > max)
            {
                max = local_max;
                amax = local_amax;
            }
        }
    }

    int64_t amP = amax / (M*N);
    int64_t amN = (amax - (amP*M*N)) / M;
    int64_t amM = amax - amP*M*N - amN*M;

    _aM[0] = amM;
    _aN[0] = amN;
    _aP[0] = amP;
    max_value[0] = max;
    return;
}


void fim_argmax0_max(const float * I,
                     size_t M, size_t N, size_t P,
                     int64_t * _aM, int64_t *_aN, int64_t *_aP,
                     float * max_value)
{
    float max = I[0];
    int64_t aM = 0;
    int64_t aN = 0;
    int64_t aP = 0;
    size_t pos = 0;
    for(size_t pp = 0; pp<P; pp++)
    {
        for(size_t nn = 0; nn<N; nn++)
        {
            for(size_t mm = 0; mm<M; mm++)
            {
                if(I[pos] > max)
                {
                    max = I[pos];
                    aM = mm;
                    aN = nn;
                    aP = pp;
                }
                pos++;
            }
        }
    }
    _aM[0] = aM;
    _aN[0] = aN;
    _aP[0] = aP;
    max_value[0] = max;
}

float fim_sum(const float * restrict A, size_t N)
{
    double sum = 0;
#pragma omp parallel for shared(A) reduction(+:sum)
    for(size_t kk = 0; kk<N; kk++)
    {
        sum+=(double) A[kk];
    }

    return (float) sum;
}

float fim_sum_double(const double * restrict A, size_t N)
{
    double sum = 0;
#pragma omp parallel for shared(A) reduction(+:sum)
    for(size_t kk = 0; kk<N; kk++)
    {
        sum+=(double) A[kk];
    }

    return sum;
}


float fim_mean(const float * restrict A, size_t N)
{
    return fim_sum(A, N)/(float) N;
}

float fim_min(const float * restrict A, size_t N)
{
    float amin = A[0];
#pragma omp parallel for reduction(min:amin)
    for(size_t kk = 0; kk<N; kk++)
    {
        if(A[kk] < amin)
            amin = A[kk];
    }
    return amin;
}

void fim_div(float * restrict  A,
             const float * restrict B,
             const float * restrict C,
             const size_t N)
/* A = B/C */
{
    size_t kk = 0;

#pragma omp parallel for shared(A, B, C)
    for(kk = 0; kk<N; kk++)
    {
        A[kk] = B[kk]/C[kk];
    }
    return;
}



void fim_minus(float * restrict  A,
               const float * restrict B,
               const float * restrict C,
               const size_t N)
/* A = B - C */
{
    size_t kk = 0;

#pragma omp parallel for shared(A,B,C)
    for(kk = 0; kk<N; kk++)
    {
        A[kk] = B[kk] - C[kk];
    }
    return;
}

void fimo_add(fimo * A, const fimo * B)
{
    assert(A->M == B->M);
    assert(A->N == B->N);
    assert(A->P == B->P);
    fim_add(A->V, B->V, A->M * A->N * A->P);
    return;
}

void fim_add(float * restrict  A,
             const float * restrict B,
             const size_t N)
/* A[kk] += B[kk] */
{
    size_t kk = 0;

#pragma omp parallel for shared(A,B)
    for(kk = 0; kk<N; kk++)
    {
        A[kk] += B[kk];
    }
    return;
}


float fim_max(const float * restrict A, size_t N)
{
    float amax = A[0];
#pragma omp parallel for reduction(max:amax)
    for(size_t kk = 0; kk<N; kk++)
    {
        if(A[kk] > amax)
            amax = A[kk];
    }
    return amax;
}

float fimo_max(const fimo * A)
{
    return fim_max(A->V, fimo_nel(A));
}

void fim_stats(const float * A, const size_t N)
{
    printf("min: %f mean: %f, max: %f\n",
           fim_min(A, N),
           fim_mean(A, N),
           fim_max(A, N));
    return;
}

float fim_percentile(const float * A, size_t N, float prct)
{
    size_t k = (size_t) (prct / 100.0 * ((double) N));
    return qselect_f32(A, N, k);
}

float fimo_percentile(fimo * A, float prct)
{
    return fim_percentile(A->V, fimo_nel(A), prct);
}

float fim_mse(float * A, float * B, size_t N)
/* mean( (A(:)-B(:)).^(1/2) )
 */
{
    double mse = 0;
#pragma omp parallel for reduction(+:mse) shared(A, B)
    for(size_t kk = 0; kk<N; kk++)
    {
        mse += pow(A[kk]-B[kk], 2);
    }
    return mse/N;
}

void fim_flipall(float * restrict T, const float * restrict A, const int64_t a1, const int64_t a2, const int64_t a3)
/* Equivalent to T = flip(flip(flip(A,1),2),3) in matlab */
{
#pragma omp parallel for shared(T, A)
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


void fim_insert(float * restrict T,
                const int64_t t1, const int64_t t2, const int64_t t3,
                const float * restrict F,
                const int64_t f1, const int64_t f2, const int64_t f3)
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

    if(f1 == t1 && f2 == t2 && f3 == t3)
    {
        memcpy(T, F, t1*t2*t3*sizeof(float));
        return;
    }

#pragma omp parallel for shared(T, F)
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

void fim_insert_ref(float * T, int64_t t1, int64_t t2, int64_t t3,
                    float * F, int64_t f1, int64_t f2, int64_t f3)
/* Insert F [f1xf2xf3] into T [t1xt2xt3] in the "upper left" corner
 * Target, T is larger than F*/
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
#pragma omp parallel for shared(F, T)
    for(int64_t pp = 0; pp<f3; pp++)
    {
        for(int64_t nn = 0; nn<f2; nn++)
        {
            for(int64_t mm = 0; mm<f1; mm++)
            {
                float x = F[mm + nn*f1 + pp*f1*f2];
                T[mm + nn*t1 + pp*t1*t2] = x;
            }
        }
    }
    return;
}


float * fim_get_cuboid(float * restrict A, const int64_t M, const int64_t N, const int64_t P,
                       const int64_t m0, const int64_t m1, const int64_t n0, const int64_t n1, const int64_t p0, const int64_t p1)
{


    ((void) P);

    /* Create a new array from V using [m0, m1]x[n0, n1]x[p0, p1] */
    int64_t m = m1-m0+1;
    int64_t n = n1-n0+1;
    int64_t p = p1-p0+1;

    float * C = fim_malloc(m*n*p*sizeof(float));

#pragma omp parallel for shared(C, A)
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

float * fim_subregion(const float * restrict A, const int64_t M, const int64_t N, const int64_t P, const int64_t m, const int64_t n, const int64_t p)
{
    ((void) P);

    /* Extract sub region starting at (0,0,0) */
    float * S = fim_malloc(m*n*p*sizeof(float));

#pragma omp parallel for shared(S, A)
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

float * fim_subregion_ref(float * A, int64_t M, int64_t N, int64_t P, int64_t m, int64_t n, int64_t p)
{
    ((void) P);
    float * S = fim_malloc(m*n*p*sizeof(float));

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

void fim_set_min_to_zero(float * restrict I, const size_t N)
{
    float min = fim_min(I, N);

#pragma omp parallel for
    for(size_t kk = 0; kk<N; kk++)
    {
        I[kk] -= min;
    }
    return;
}

void fim_mult_scalar(float * restrict I, size_t N, float x)
{
#pragma omp parallel for shared(I)
    for(size_t kk = 0; kk < N ; kk++)
    {
        I[kk]*=x;
    }
    return;
}

void fimo_mult_scalar(fimo * A, float x)
{
    fim_mult_scalar(A->V, fimo_nel(A), x);
}

int fimo_mult_image(fimo * A, const fimo * B)
{
    int mode = 0;
    if(A->M != B->M)
    {
        return -1;
    }
    if(A->N != B->N)
    {
        return -1;
    }
    if(A->P > 1)
    {
        if(B->P == 1)
        {
            mode = 1;
        }
        if(B->P == A->P)
        {
            mode = 2;
        }
    }

    // per plane
    if(mode == 1)
    {
#pragma omp parallel for
        for(size_t pp = 0; pp < A->P; pp++)
        {
            float * P = A->V + pp*A->M*A->N;
            for(size_t kk = 0; kk < fimo_nel(B); kk++)
            {
                P[kk] *= B->V[kk];
            }
        }
        return 0;
    }

    // same number of elements
    if(mode == 2)
    {
#pragma omp parallel for
        for(size_t kk = 0; kk < fimo_nel(A); kk++)
        {
            A->V[kk] *= B->V[kk];
        }
        return 0;
    }
    return -1;
}

int fimo_div_image(fimo * A, const fimo * B)
{
    int mode = 0;
    if(A->M != B->M)
    {
        return -1;
    }
    if(A->N != B->N)
    {
        return -1;
    }
    if(A->P == 1)
    {
        if(B->P == 1)
        {
            mode = 2;
        }
    }
    if(A->P > 1)
    {
        if(B->P == 1)
        {
            mode = 1;
        }
        if(B->P == A->P)
        {
            mode = 2;
        }
    }

    // per plane
    if(mode == 1)
    {
#pragma omp parallel for
        for(size_t pp = 0; pp < A->P; pp++)
        {
            float * P = A->V + pp*A->M*A->N;
            for(size_t kk = 0; kk < fimo_nel(B); kk++)
            {
                P[kk] /= B->V[kk];
            }
        }
        return 0;
    }

    // same number of elements
    if(mode == 2)
    {
#pragma omp parallel for
        for(size_t kk = 0; kk < fimo_nel(A); kk++)
        {
            A->V[kk] /= B->V[kk];
        }
        return 0;
    }
    return -1;
}


void fim_project_positive(float * I, size_t N)
{
#pragma omp parallel for shared(I)
    for(size_t kk = 0; kk < N ; kk++)
    {
        I[kk] < 0 ? I[kk] = 0 : 0;
    }
    return;
}

void fim_add_scalar(float * restrict I, size_t N, float x)
{
#pragma omp parallel for shared(I)
    for(size_t kk = 0; kk < N ; kk++)
    {
        I[kk]+=x;
    }
    return;
}


void fim_normalize_sum1(float * restrict psf, int64_t M, int64_t N, int64_t P)
{
    const size_t pMNP = M*N*P;;
    double psf_sum = 0;
#pragma omp parallel for shared(psf) reduction(+:psf_sum)
    for(size_t kk = 0; kk<pMNP; kk++)
    { psf_sum += psf[kk]; }
    //  printf("psf_sum: %f\n", psf_sum);
#pragma omp parallel for shared(psf)
    for(size_t kk = 0; kk<pMNP; kk++)
    { psf[kk]/=psf_sum; }
    return;
}

float * fim_copy(const float * restrict V, const size_t N)
// Return a newly allocated copy of V
{
    float * C = fim_malloc(N*sizeof(float));

    memcpy(C, V, N*sizeof(float));
    return C;
}

fimo * fimo_copy(const fimo * restrict F)
{
    fimo * C = malloc(sizeof(fimo));
    assert(C != NULL);
    C->V = fim_copy(F->V, F->M*F->N*F->P);
    C->M = F->M;
    C->N = F->N;
    C->P = F->P;
    return C;
}

float * fim_zeros(const size_t N)
// Allocate and return an array of N floats
{
    assert(N > 0);
    float * A = fim_malloc(N*sizeof(float));
    assert(A[0] == 0);
    assert(A[N-1] == 0);
    return A;
}

fimo * fimo_zeros(const size_t M, const size_t N, const size_t P)
// Allocate and return an array of N floats
{
    size_t n = M*N*P;
    fimo * F = malloc(sizeof(fimo));
    assert(F != NULL);
    F->V = fim_malloc(n*sizeof(float));

    F->M = M;
    F->N = N;
    F->P = P;
    //memset(F->V, 0, n*sizeof(float));
    return F;
}

size_t fimo_nel(const fimo * F)
{
    return F->M*F->N*F->P;
}

float fimo_sum(const fimo * F)
{
    return fim_sum(F->V, fimo_nel(F));
}

float * fim_constant(const size_t N, const float value)
// Allocate and return an array of N floats sets to a constant value
{
    float * A = fim_malloc(N*sizeof(float));

#pragma omp parallel for shared(A)
    for(size_t kk = 0; kk<N; kk++)
    {
        A[kk] = value;
    }
    return A;
}

void fim_circshift(float * restrict A,
                   const int64_t M, const int64_t N, const int64_t P,
                   const int64_t sm, const int64_t sn, const int64_t sp)
/* Shift the image A [MxNxP] by sm, sn, sp in each dimension */
{



    int nThreads = 0;

    /* Start a parallel region to figure out how many threads that will be used
       and how much memory to allocate for the buffers */
#pragma omp parallel
    {
        nThreads = omp_get_num_threads();
    }

    const size_t bsize = fmax(fmax(M, N), P);
    float * restrict buf = fim_malloc(bsize*sizeof(float)*nThreads);
    assert(buf != NULL);


    /* Dimension 1 */
#pragma omp parallel for shared(A)
    for(int64_t cc = 0; cc<P; cc++)
    {
        float * tbuf = buf + bsize*omp_get_thread_num();
        for(int64_t bb = 0; bb<N; bb++)
        {
            //shift_vector(A + bb*M + cc*M*N, 1, M, sm);
            shift_vector_buf(A + bb*M + cc*M*N, // start
                             1, // stride
                             M, // number of elements
                             sm, // shift
                             tbuf); // buffer
        }
    }

    /* Dimension 2 */
#pragma omp parallel for shared(A)
    for(int64_t cc = 0; cc<P; cc++)
    {
        float * tbuf = buf + bsize*omp_get_thread_num();
        //printf("Thread num: %d\n", omp_get_thread_num());
        for(int64_t aa = 0; aa<M; aa++)
        {
            //shift_vector(A + aa+cc*M*N, M, N, sn);
            shift_vector_buf(A + aa+cc*M*N, M, N, sn, tbuf);
        }
    }

    /* Dimension 3 */
#pragma omp parallel for shared(A)
    for(int64_t bb = 0; bb<N; bb++)
    {
        float * tbuf = buf + bsize*omp_get_thread_num();
        for(int64_t aa = 0; aa<M; aa++)
        {
            //shift_vector(A + aa+bb*M, M*N, P, sp);
            shift_vector_buf(A + aa+bb*M, M*N, P, sp, tbuf);
        }
    }

    fim_free(buf);

#if 0
    char name[100];
    sprintf(name, "shift%d.raw", nThreads);
    printf("Debugdump to %s\n", name);
    FILE *fid = fopen(name, "w");
    fwrite(A, sizeof(float), M*N*P, fid);
    fclose(fid);
#endif

    return;
}

/* A linear interpolation kernel for -1<delta<1 */
static float * kernel_linear_shift(float delta, int * _N)
{
    if(delta == 0)
    {
        *_N = 0;
        return NULL;
    }
    *_N = 3; /* Always three elements */
    float * K = fim_malloc(3*sizeof(float));
    assert(K != NULL);
    K[2] = delta;
    K[2] < 0 ? K[2] = 0 : 0;
    K[1] = fabs(1.0-fabs(delta));
    K[0] = -delta;
    K[0] < 0 ? K[0] = 0 : 0;
    printf("Kernel : [%f, %f, %f]\n", K[0], K[1], K[2]);
    return K;
}

/* Shift the image A [MxNxP] by dm, dn, dp in each dimension,
 * What is outside of the image is interpreted as zero */
void fim_shift(float * restrict A,
               const int64_t M, const int64_t N, const int64_t P,
               const float dm, const float dn, const float dp)
{
    int nThreads = 0;


    int sm = round(dm);
    int sn = round(dn);
    int sp = round(dp);

    float deltam = dm-sm;
    float deltan = dn-sn;
    float deltap = dp-sp;

    int nkernelx = -1, nkernely = -1, nkernelz = -1;
    float * kernelx = kernel_linear_shift(deltam, &nkernelx);
    float * kernely = kernel_linear_shift(deltan, &nkernely);
    float * kernelz = kernel_linear_shift(deltap, &nkernelz);


    //    printf("X Shift: %f = %d+%f\n", dm, sm, deltam);
    //    printf("Y Shift: %f = %d+%f\n", dn, sn, deltan);
    //    printf("Z Shift: %f = %d+%f\n", dp, sp, deltap);

    /* Start a parallel region to figure out how many threads that will be used
       and how much memory to allocate for the buffers */
#pragma omp parallel
    {
        nThreads = omp_get_num_threads();
    }

    const size_t bsize = fmax(fmax(M, N), P);
    float * restrict buf = fim_malloc(bsize*sizeof(float)*nThreads);
    assert(buf != NULL);


    /* Dimension 1 */
#pragma omp parallel for
    for(int64_t cc = 0; cc<P; cc++)
    {
        float * tbuf = buf + bsize*omp_get_thread_num();
        for(int64_t bb = 0; bb<N; bb++)
        {
            //shift_vector(A + bb*M + cc*M*N, 1, M, sm);
            shift_vector_float_buf(A + bb*M + cc*M*N, // start
                                   1, // stride
                                   M, // number of elements
                                   sm, // shift
                                   kernelx,
                                   nkernelx,
                                   tbuf); // buffer
        }
    }

    /* Dimension 2 */
#pragma omp parallel for
    for(int64_t cc = 0; cc<P; cc++)
    {
        float * tbuf = buf + bsize*omp_get_thread_num();
        //printf("Thread num: %d\n", omp_get_thread_num());
        for(int64_t aa = 0; aa<M; aa++)
        {
            //shift_vector(A + aa+cc*M*N, M, N, sn);
            shift_vector_float_buf(A + aa+cc*M*N,
                                   M,
                                   N,
                                   sn, kernely, nkernely, tbuf);
        }
    }

    /* Dimension 3 */
#pragma omp parallel for
    for(int64_t bb = 0; bb<N; bb++)
    {
        float * tbuf = buf + bsize*omp_get_thread_num();
        for(int64_t aa = 0; aa<M; aa++)
        {
            //shift_vector(A + aa+bb*M, M*N, P, sp);
            shift_vector_float_buf(A + aa+bb*M,
                                   M*N,
                                   P,
                                   sp, kernelz, nkernelz, tbuf);
        }
    }

    fim_free(buf);

    fim_free(kernelx);
    fim_free(kernely);
    fim_free(kernelz);

#if 0
    char name[100];
    sprintf(name, "shift%d.raw", nThreads);
    printf("Debugdump to %s\n", name);
    FILE *fid = fopen(name, "w");
    fwrite(A, sizeof(float), M*N*P, fid);
    fclose(fid);
#endif

    return;
}


static int64_t mod_int(const int64_t a, const int64_t b)
{
    int64_t r = a % b;
    return r < 0 ? r + b : r;
}


/* Shift vector by interpolation */
void shift_vector_float_buf(float * restrict V, // data
                            const int64_t S, // stride
                            const int64_t N, // elements
                            int n, // integer shift
                            float * restrict kernel, // centered kernel used for sub pixels shift
                            const int nkernel, // kernel size (odd!)
                            float * restrict buffer)
{
    // 1. Sub pixel shift by convolution (conv1) of signal and kernel
    // 2. Integer part of the shift
    /* TODO: Interpolation here ... by interpolation kernel?*/
    /* First integer part and then sub pixel? */
    if(kernel == NULL)
    {
        for(size_t pp = 0; pp < (size_t) N; pp++)
        {
            buffer[pp] = V[pp*S];
        }
    } else {
        fim_conv1_vector(V, S, buffer, N, kernel, nkernel, 1);
    }

    for(size_t pp = 0; pp<(size_t) N; pp++)
    {
        int q = pp+n;
        if(q>= 0 && q<N)
        {
            buffer[pp] = V[q*S];
        } else {
            buffer[pp] = 0;
        }
    }

    for(size_t pp = 0; pp < (size_t) N; pp++)
    {
        V[pp*S] = buffer[pp];
    }
    return;
}


void shift_vector_buf(float * restrict V, // data
                      const int64_t S, // stride
                      const int64_t N, // elements
                      int64_t k, // shift
                      float * restrict buffer)
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

void shift_vector(float * restrict V,
                  const int64_t S,
                  const int64_t N,
                  const int64_t k)
/* Circular shift of a vector of length N with stride S by step k */
{

    float * buffer = fim_malloc(N*sizeof(float));
    assert(buffer != NULL);
    shift_vector_buf(V, S, N, k, buffer);
    fim_free(buffer);
    return;
}


float * fim_expand(const float * restrict in,
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
    assert(in != NULL);

    float * out = fim_malloc(M*N*P*sizeof(float));

    for(size_t kk = 0; kk < (size_t) M*N*P; kk++)
        out[kk] = 0;
    fim_insert(out, M, N, P, in, pM, pN, pP);
    return out;
}

void fim_flipall_ut()
{

    float * a = fim_malloc(3*3*3*sizeof(float));
    assert(a != NULL);
    float * b = fim_malloc(3*3*3*sizeof(float));
    assert(b != NULL);
    float * c = fim_malloc(3*3*3*sizeof(float));
    assert(c != NULL);

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

    fim_free(a); fim_free(b); fim_free(c);
    return;
}

void shift_vector_ut()
{
    int64_t N = 5;
    int64_t S = 1; // stride
    float * V = fim_malloc(N*sizeof(float));
    assert(V != NULL);

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
    fim_free(V);
}

void fim_invert(float * restrict A, size_t N)
{
    for(size_t kk = 0; kk < N; kk++)
    {
        A[kk] = 1-A[kk];
    }
    return;
}

static void show_vec(float * V, int N)
{
    for(int kk = 0; kk<N; kk++)
    {
        printf("%f ", V[kk]);
    }
    printf("\n");
}

void fim_local_sum_ut()
{
    printf("fim_local_sum_ut\n");
    size_t M = 3;
    size_t N = 5;
    float * A = fim_constant(M*N, 1);
    fim_show(A, M, N, 1);
    int filtM = 2; int filtN = 2;
    float * L = fim_local_sum(A, M, N, filtM, filtN);
    size_t lM = M + filtM-1;
    size_t lN = N + filtN-1;
    fim_show(L, lM, lN, 1);
    fim_free(A);
    fim_free(L);

}

void fim_cumsum_ut()
{
    size_t M = 3;
    size_t N = 5;
    float * A = fim_zeros(M*N);
    fim_show(A, M, N, 1);
    fim_cumsum(A, M, N, 0);
    for(size_t kk = 0; kk<M*N; kk++)
    {
        assert(A[kk] == 0);
    }
    for(size_t kk = 0; kk<M*N; kk++)
    {
        A[kk] = kk;
    }
    fim_show(A, M, N, 1);
    printf("After cumsum along dim 0\n");
    fim_cumsum(A, M, N, 0);
    fim_show(A, M, N, 1);
    printf("After cumsum along dim 1\n");
    fim_cumsum(A, M, N, 1);
    fim_show(A, M, N, 1);
    fim_free(A);
}

void fim_xcorr2_ut()
{
    size_t M = 3;
    size_t N = 7;
    float * T = fim_zeros(M*N);
    float * A = fim_zeros(M*N);
    T[0] = 1;
    T[1] = 1;
    A[M+1] = 1;
    A[M+2] = 0.5;
    for(size_t kk = 0; kk<M*N; kk++)
    {
        A[kk] = (float) rand()/ (float) RAND_MAX;
        T[kk] = A[kk];
    }
    printf("A=\n");
    fim_show(A, M, N, 1);
    float * X = fim_xcorr2(T, A, M, N);
    printf("T=\n");
    fim_show(T, M, N, 1);
    printf("A=\n");
    fim_show(A, M, N, 1);
    printf("xcorr2(T,A)=\n");
    fim_show(X, 2*M-1, 2*N-1, 1);
    fim_free(A);
    fim_free(T);
    fim_free(X);
}

void fim_otsu_ut()
{
    size_t M = 10;
    size_t N = 10;
    size_t P = 1;
    float * Im = calloc(M*N*P, sizeof(float));
    assert(Im != NULL);
    for(size_t kk = 0; kk<M*N; kk++)
    {
        if(kk < M*N/2)
        {
            Im[kk] = 1 + 1*((float) rand() / (float) RAND_MAX - 0.5);
        } else {
            Im[kk] = 2 + 1*((float) rand() / (float) RAND_MAX - 0.5);
        }
    }
    fim_show(Im, M, N, P);
    float * B = fim_otsu(Im, M, N);
    fim_show(B, M, N, P);
    fim_free(Im);
    fim_free(B);
}

void fim_conncomp6_ut()
{
    printf("fim_conncomp6_ut\n");
    size_t M = 10;
    size_t N = 10;
    float * Im = calloc(M*N, sizeof(float));
    assert(Im != NULL);
    for(size_t kk = 0; kk < M*N/4; kk++)
    {
        size_t idx = rand() % (M*N);
        Im[idx] = 1;
    }
    printf("Input image\n");
    fim_show(Im, M, N, 1);
    int * L = fim_conncomp6(Im, M, N);
    printf("Labeled array\n");
    fim_show_int(L, M, N, 1);
    fim_free(L);
    fim_free(Im);
}

void fim_conv1_vector_ut()
{
    printf("->fim_conv1_vector_ut()\n");
    int normalized = 1;
    const size_t nV = 5;
    float * V = malloc(nV*sizeof(float));
    assert(V != NULL);
    for(size_t kk = 0; kk<nV; kk++)
    {
        V[kk] = kk+1;
    }
    int stride = 1;
    float * W = NULL;
    printf("V=\n");
    fim_show(V, 1, nV, 1);
    for(int nK = 3; nK<= (int) nV; nK+=2)
    {
        for(size_t kk = 0; kk<nV; kk++)
        {
            V[kk] = kk+1;
        }
        float * K = malloc(nK*sizeof(float));
        assert(K != NULL);
        for(int kk = 0; kk<nK; kk++)
        {
            K[kk] = 1.0;
        }
        printf("K=");
        fim_show(K, 1, nK, 1);
        fim_conv1_vector(V, stride, W, nV, K, nK, normalized);
        printf("V*K=");
        fim_show(V, 1, nV, 1);
        fim_free(K);
    }
    fim_free(V);
}

void fim_LoG_ut()
{
    printf("-> fim_LoG_ut\n");
    struct timespec tstart, tend;
    size_t M = 256;
    size_t N = 512;
    size_t P = 128;
    float sigma_l = 3.1;
    float sigma_a = 1.1;
    float * V = fim_malloc(M*N*P*sizeof(float));

    for(size_t kk = 0; kk<M*N*P; kk++)
    {
        V[kk] = kk % 100;
    }


    fimo * T = fim_image_from_array(V, M, N, P);
    fimo * S1 = fim_shiftdim(T);
    fimo * S2 = fim_shiftdim(S1);
    fimo * S3 = fim_shiftdim(S2);
    for(size_t kk = 0; kk<M*N*P; kk++)
    {
        if(T->V[kk] != S3->V[kk])
        {
            printf("fim_shiftdim does not work at index %zu %f != %f\n",
                   kk, T->V[kk], S3->V[kk]);
            exit(EXIT_FAILURE);
        }
    }
    fimo_free(T);
    T = NULL;
    assert(S1->V != S2->V);
    fimo_free(S1);
    S1 = NULL;
    fimo_free(S2);
    S2 = NULL;
    fimo_free(S3);


    for(size_t kk = 0; kk<M*N*P; kk++)
    {
        V[kk] = 0;
    }
    V[(M/2) + M*(N/2) + M*N*(P/2)] = 1;


    float * Vpre = fim_copy(V, M*N*P);

    dw_gettime(&tstart);
    float * LoG2 = fim_LoG_S(V, M, N, P, sigma_l, sigma_a);
    dw_gettime(&tend);
    for(size_t kk = 0; kk<M*N*P; kk++)
    {
        if(V[kk] != Vpre[kk])
        {
            fprintf(stderr, "fim_LoG_S alters the input array\n");
            exit(EXIT_FAILURE);

        }
    }
    float tLoG_S = clockdiff(&tend, &tstart);

    dw_gettime(&tstart);
    float * LoG = fim_LoG(V, M, N, P, sigma_l, sigma_a);
    dw_gettime(&tend);
    for(size_t kk = 0; kk<M*N*P; kk++)
    {
        if(V[kk] != Vpre[kk])
        {
            fprintf(stderr, "fim_LoG alters the input array\n");
            exit(EXIT_FAILURE);
        }
    }
    fim_free(Vpre);

    float tLoG = clockdiff(&tend, &tstart);

    printf("LoG : %f s (min %f, max %f)\n", tLoG,
           fim_min(LoG, M*N*P),
           fim_max(LoG, M*N*P));
    printf("LoG_S : %f s\n", tLoG_S);

    float maxabs = fabs(LoG[0]-LoG2[0]);
    float * Delta = fim_malloc(M*N*P*sizeof(float));

    for(size_t kk = 0; kk<M*N*P; kk++)
    {
        float diff = fabs(LoG[kk]-LoG2[kk]);
        diff > maxabs ? maxabs = diff : 0;
        Delta[kk] = diff;
    }
    printf("Max abs difference: %e\n", maxabs);

    if(1){
        printf("Writing to V.tif\n");
        fim_tiff_write_float("V.tif", V, NULL,  N, M, P);
        printf("Writing to Diff.tif\n");
        fim_tiff_write_float("Diff.tif", Delta, NULL,  N, M, P);
        printf("Writing to LoG.tif\n");
        fim_tiff_write_float("LoG.tif", LoG, NULL,  N, M, P);
        printf("Writing to LoG_S.tif\n");
        fim_tiff_write_float("LoG_S.tif", LoG2, NULL, N, M, P);
    }
    fim_free(Delta);

    if(M*N*P < 200)
    {
        printf("LoG = \n");
        fim_show(LoG, M, N, P);
        printf("LoG_S = \n");
        fim_show(LoG2, M, N, P);
    }

    dw_gettime(&tstart);
    float * LoGS2 = fim_LoG_S2(V, M, N, P, sigma_l, sigma_a);
    dw_gettime(&tend);
    free(LoGS2);
    printf("fim_LoG_S2 took %f s\n", clockdiff(&tend, &tstart));
    fim_free(V);
    fim_free(LoG);
    fim_free(LoG2);
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


static float fim_interp_symmetric_mirror(const float * V, int64_t nV,
                                         int64_t stride,
                                         int64_t idx)
{
    idx < 0 ? idx = -idx : 0;
    idx >= nV ? idx = (nV-1)-(idx+1-nV) : 0;
    assert(idx >= 0);
    assert(idx < nV);
    return V[stride*idx];
}

static float fim_interp_periodic(const float * V, int64_t nV,
                                 int64_t stride,
                                 int64_t idx)
{
    return V[stride*(idx % nV)];
}

void fimo_conv1_x(fimo * V, fimo * K, fim_boundary_condition bc)
{

#pragma omp parallel
    {
        float * B = fim_malloc(V->M*sizeof(float));

        for(size_t pp = 0; pp<V->P; pp++)
        {
#pragma omp for schedule(dynamic)
            for(size_t nn = 0; nn<V->N; nn++)
            {
                float * line = V->V+pp*(V->M*V->N) + nn*(V->M);
                fim_conv1(line, V->M,
                          1, // stride
                          K->V,
                          fimo_nel(K),
                          B,
                          bc);
            }
        }
        free(B);
    }
}

/* Boundary conditions like here:
 * https://diplib.org/diplib-docs/boundary.html#dip-BoundaryCondition */
void fim_conv1(float * restrict V, const size_t nV, const int stride,
               const float * restrict K, const size_t nK,
               float * restrict buffer,
               fim_boundary_condition bc)
{
    assert(V != NULL);
    assert(buffer != NULL);
    if(nK % 2 == 0)
    {
        fprintf(stderr,
                "fim_conv1 error: will only work with kernels of odd size\n");
        exit(EXIT_FAILURE);
        return;
    }

    if( (nK+1)/2 > nV)
    {
        fprintf(stderr, "\n"
                "fim_conv1 errror:\n"
                "   Kernel size: %zu, Vector size: %zu\n"
                "   The kernel is too large: unable to perform the convolution\n"
                "   Even if this could be done, it would probably be a bad idea\n"
                "   The program will crash now\n",
                nK, nV);
        exit(EXIT_FAILURE);
    }

    memset(buffer, 0, nV*sizeof(float));

    /* Tight case, when the kernel is larger than the input image */
    const int64_t mid = (nK-1)/2;
    if(nK > nV)
    {
        if(bc == FIM_BC_SYMMETRIC_MIRROR)
        {
            for(int64_t ii = 0 ; ii < (int64_t) nV; ii++)
            {
                for(int64_t kk = 0; kk < (int64_t ) nK; kk++)
                {
                    int64_t idx = ii + kk - mid;
                    buffer[ii] += K[kk]*fim_interp_symmetric_mirror(V, nV, stride, idx);
                }
            }
            return;
        }
        if(bc == FIM_BC_PERIODIC)
        {
            for(int64_t ii = 0 ; ii < (int64_t) nV; ii++)
            {
                for(int64_t kk = 0; kk < (int64_t ) nK; kk++)
                {
                    int64_t idx = ii + kk - mid;
                    buffer[ii] += K[kk]*fim_interp_periodic(V, nV, stride, idx);
                }
            }
            return;
        }
        fprintf(stderr, "\n"
                "fim_conv1 errror:\n"
                "   Kernel size: %zu, Vector size: %zu\n"
                "   The boundary condition is not implemented for this problem size\n"
                "   The program will crash now\n",
                nK, nV);
        exit(EXIT_FAILURE);
        return;
    }

    /* Normal case, the kernel is smaller than the vector */

    if(K == NULL) { return; }
    int buffer_allocation = 0;
    if(buffer == NULL)
    {
        buffer_allocation = 1;
        buffer = calloc(nV, sizeof(float));
    }

    double Wtotal = 0;
    if(bc == FIM_BC_WEIGHTED)
    {
        for(size_t kk = 0 ; kk<nK; kk++)
        {
            Wtotal += K[kk];
        }
    }

    /* Part I: The kernel is overlapping the edge
       if bc == FIM_BC_VALID, this part is skipped */


    if(bc == FIM_BC_SYMMETRIC_MIRROR)
    {
        for(int64_t ii = 0 ; ii < mid; ii++)
        {
            for(int64_t kk = 0; kk < (int64_t) nK; kk++)
            {
                int64_t idx = ii + kk - mid;
                idx < 0 ? idx = -idx : 0;
                buffer[ii] += K[kk]*V[stride*idx];
            }
        }
    }

    if(bc == FIM_BC_ZEROS)
    {
        for(int64_t ii = 0 ; ii < mid; ii++)
        {
            for(int64_t kk = 0; kk < (int64_t) nK; kk++)
            {
                int64_t idx = ii + kk - mid;
                if(idx >= 0)
                {
                    buffer[ii] += K[kk]*V[stride*idx];
                }
            }
        }
    }

    if(bc == FIM_BC_PERIODIC)
    {
        for(int64_t ii = 0 ; ii < mid; ii++)
        {
            for(int64_t kk = 0; kk < (int64_t) nK; kk++)
            {
                int64_t idx = (ii + kk - mid);
                idx < 0 ? idx = nV +idx : 0;
                buffer[ii] += K[kk]*V[stride*idx];
            }
        }
    }

    if(bc == FIM_BC_WEIGHTED)
    {
        for(int64_t ii = 0 ; ii < mid; ii++)
        {
            float W = 0;
            for(int64_t kk = 0; kk < (int64_t) nK; kk++)
            {
                int64_t idx = (ii + kk - mid);
                if(idx >= 0)
                {
                    buffer[ii] += K[kk]*V[stride*idx];
                    W+=K[kk];
                }
            }
            buffer[ii]*=Wtotal/W;
        }
    }

    /* Part II: Central, the kernels is inside the image */
    for(int64_t ii = mid ; ii < (int64_t) nV-mid; ii++)
    {
        for(int64_t kk = 0; kk < (int64_t) nK; kk++)
        {
            int64_t idx = ii + kk - mid;
            //printf("kk: %ld, idx: %ld, ii: %ld\n", kk, idx, ii);
            //fflush(stdout);
            buffer[ii] += K[kk]*V[stride*idx];
        }
    }

    /* Part III: Last, where kernel is beyond the last pixel
       if bc == FIM_BC_VALID, this part is skipped */
    if(bc == FIM_BC_SYMMETRIC_MIRROR)
    {
        int64_t inV = nV;
        for(int64_t ii = inV-mid ; ii < inV; ii++)
        {
            for(int64_t kk = 0; kk < (int64_t) nK; kk++)
            {
                int64_t idx = ii + kk - mid;
                idx >= inV ? idx = (inV-1)-(idx+1-inV) : 0;
                buffer[ii] += K[kk]*V[stride*idx];
            }
        }
    }

    if(bc == FIM_BC_ZEROS)
    {
        int64_t inV = nV;
        for(int64_t ii = inV-mid ; ii < inV; ii++)
        {
            for(int64_t kk = 0; kk < (int64_t) nK; kk++)
            {
                int64_t idx = ii + kk - mid;
                if(idx < inV)
                { buffer[ii] += K[kk]*V[stride*idx]; }
            }
        }
    }

    if(bc == FIM_BC_PERIODIC)
    {
        int64_t inV = nV;
        for(int64_t ii = inV-mid ; ii < inV; ii++)
        {
            for(int64_t kk = 0; kk < (int64_t) nK; kk++)
            {
                int64_t idx = (ii + kk - mid) % nV;
                buffer[ii] += K[kk]*V[stride*idx];
            }
        }
    }

    if(bc == FIM_BC_WEIGHTED)
    {
        int64_t inV = nV;
        for(int64_t ii = inV-mid ; ii < inV; ii++)
        {
            float W = 0;
            assert(buffer[ii] == 0);
            for(int64_t kk = 0; kk < (int64_t) nK; kk++)
            {
                int64_t idx = (ii + kk - mid);
                if(idx < inV)
                {
                    buffer[ii] += K[kk]*V[stride*idx];
                    W+=K[kk];
                }
            }
            buffer[ii]*=Wtotal/W;
        }
    }

    /* Write back the result to the input array */
    for(size_t ii = 0; ii<nV; ii++)
    {
        V[ii*stride] = buffer[ii];
    }

    if(buffer_allocation)
    {
        free(buffer);
    }
    return;
}

void fim_conv1_vector(float * restrict V, const int stride, float * restrict W,
                      const size_t nV,
                      const float * restrict K, const size_t nKu, const int normalized)
{
    if(V == NULL || K == NULL)
    {
        return;
    }

    if(nKu > nV)
    {
        fprintf(stderr,
                "fim_conv1_vector: error - kernel can't be longer than data\n");
        return;
    }

    /* Allocate buffer if not provided */
    int Walloc = 0;
    if(W == NULL)
    {
        W = fim_malloc(nV*sizeof(float));

        Walloc = 1;
    }


    const size_t k2 = (nKu-1)/2;

    size_t bpos = 0;

    if(nKu > nV)
    {
        for(size_t kk = 0; kk<nV; kk++)
        {
            double acc0 = 0;
            double kacc = 0;
            size_t from = kk-k2;
            k2 > kk ? from = 0 : 0;
            size_t to = kk+k2+1;
            to > nV ? to = nV: 0;
            for(size_t ll = from; ll < to ; ll++)
            {
                int kpos = ll-kk+k2;
                int vpos = ll;
                //printf("kk=%d, ll=%d, kpos=%d vpos=%d V=%f, K=%f\n", kk, ll, kpos, vpos, V[vpos*stride], K[kpos]);
                acc0 += V[vpos*stride]*K[kpos];
                kacc += K[kpos];
            }
            //printf("acc0=%f\n", acc0);
            if(normalized)
            {
                W[kk] = acc0/kacc;
            } else {
                W[kk] = acc0;
            }
        }
    }

    /* Kernel is smaller than signal */
    if(nKu <= nV)
    {
        /* Crossing the first edge */
        for(size_t vv = 0; vv<k2; vv++)
        {
            double acc0 = 0;
            double kacc = 0;
            for(size_t kk = k2-vv; kk < nKu; kk++)
            {

                assert((vv-k2+kk) < nV);
                acc0 = acc0 + K[kk]*V[(vv-k2+kk)*stride];

                kacc += K[kk];
            }

            if(normalized)
            {
                W[bpos++] = acc0/kacc;
            } else {
                W[bpos++] = acc0;
            }
        }

        /* Central part where K fits completely */
        for(size_t vv = k2 ; vv+k2 < nV; vv++)
        {
            double acc = 0;
            for(size_t kk = 0; kk < nKu; kk++)
            {
                size_t vpos = ((vv-k2)+kk)*stride;
                assert(vpos/stride < nV);
                acc = acc + K[kk]*V[vpos];
            }
            W[bpos++] = acc;
        }

        /* Last part */
        for(size_t vv = nV-k2; vv< nV; vv++)
        {
            double kacc = 0;
            double acc0 = 0;
            for(size_t kk = 0; kk<nV-vv+k2; kk++)
            {
                acc0 = acc0 + K[kk]*V[(vv-k2+kk)*stride];
                kacc += K[kk];
            }
            if(normalized)
            {
                W[bpos++] = acc0/kacc;
            } else {
                W[bpos++] = acc0;
            }
        }
    }

    /* Write back */
    for(size_t pp = 0; pp<nV; pp++)
    {
        V[pp*stride] = W[pp];
    }

    if(Walloc)
    {
        fim_free(W);
    }
    return;
}

static void flip_sign(float * X, size_t N)
{
    for(size_t kk = 0; kk<N; kk++)
    {
        X[kk]*=-1;
    }
}

/* A Gaussian kernel
 * Note: Always normalized to have sum 1.0 */
static float * gaussian_kernel(float sigma, size_t * nK)
{
    /* Determine the size so that most of the signal is captured */

    int len = ceil(3.0*sigma);
    len < 1 ? len = 1 : 0;

    int N = 2*len + 1;

    float * K = fim_malloc(N*sizeof(float));

    float mid = (N-1)/2;

    float s2 = pow(sigma, 2);
    float k0 = 1.0/sigma/sqrt(2.0*M_PI);
    double sum = 0;
    for(int kk = 0; kk<N; kk++)
    {
        float x = (float) kk - mid;
        K[kk] = k0*exp(-0.5*pow(x,2)/s2);
        sum+=K[kk];
    }
    for(int kk = 0; kk<N; kk++)
    {
        K[kk] /= sum;
    }

    nK[0] = N;
    return K;
}

static float * gaussian_kernel_d1(float sigma, size_t * nK)
{
    /* First derivative of a Gaussian kernel */

    /* Determine the size so that most of the signal is captured */
    int len = ceil(3.0*sigma);
    len < 1 ? len = 1 : 0;

    /* Total number of elements */
    int N = 2*len + 1;

    float * K = fim_malloc(N*sizeof(float));

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

    for(int kk = 0; kk<N; kk++)
    {
        float x = (float) kk-mid;
        K[kk] *= -x/s2;
    }

    nK[0] = N;
    return K;
}


/** @brief 1D Laplacian of Gaussian
 *
 * The filter is centered in the returned array. nK
 * is the number of elements.
 *
 * Note: Besides being truncated, the Laplacian is only sampled at the
 * middle of the pixels, not integrated. Hence it might not integrate
 * to 0.
 */
static float * gaussian_kernel_d2(float sigma, size_t * nK)
{
    /* d2/dx2 Gaussian kernel */
    int m = ceil(4.0*sigma);
    m < 1 ? m = 1 : 0;
    int n = 2*m + 1;

    float * K = fim_malloc(n*sizeof(float));

    nK[0] = n;

    float s2 = pow(sigma, 2.0);
    float s4 = pow(sigma, 4.0);
    float k0 = 1.0/sigma/sqrt(2.0*M_PI);
    for(int kk = 0; kk < n; kk++)
    {
        float x = (float) kk- (float) m;
        float G = k0*exp(-0.5*pow(x,2)/s2);
        K[kk] = (pow(x, 2.0)-s2)/s4*G;
    }

    return K;
}


void fim_gsmooth_old(float * restrict V, size_t M, size_t N, size_t P, float sigma)
{
    /* Convolve V by an isotropic Gaussian
     * implemented as a separated convolution. */

    if(sigma < 0)
    {
        printf("fim_gsmooth sigma=%f does not make sense.", sigma);
    }

    size_t nW = max_size_t(M, max_size_t(N, P));
    //printf("fim_gsmooth: M: %zu, N: %zu, P: %zu, nW: %zu\n", M, N, P, nW); fflush(stdout);

    /* Temporary storage/buffer for conv1 */
    float * W = fim_malloc(nW*sizeof(float));


    /* Create a kernel  */
    size_t nK = 0;
    float * K = gaussian_kernel(sigma, &nK);
    assert(nK > 0);

    // X
    for(size_t pp = 0; pp < P; pp++)
    {
        for(size_t nn = 0; nn < N; nn++)
        {
            fim_conv1_vector(V+pp*(M*N)+nn*M, 1, W, M, K, nK, 1);
        }
    }

    if(1){
        // Y
        for(size_t pp = 0; pp<P; pp++)
        {
            for(size_t mm = 0; mm<M; mm++)
            {
                fim_conv1_vector(V + pp*(M*N) + mm, M, W, N, K, nK, 1);
            }
        }
    }

    if(1){
        // Z
        for(size_t mm = 0; mm<M; mm++)
        {
            for(size_t nn = 0; nn<N; nn++)
            {
                fim_conv1_vector(V+mm+M*nn, M*N, W, P, K, nK, 1);
            }
        }
    }

    fim_free(W);
    return;
}

void fim_gsmooth_aniso(float * restrict V,
                       size_t M, size_t N, size_t P,
                       float lsigma, float asigma)
{
    // TODO: use the shiftdim approach as in fim_LoG_S

    if(lsigma > 0)
    {
        size_t nKl = 0;
        float * Kl = gaussian_kernel(lsigma, &nKl);
        fim_convn1(V, M, N, P, Kl, nKl, 0, 1);
        fim_convn1(V, M, N, P, Kl, nKl, 1, 1);
        fim_free(Kl);
    }
    if(asigma > 0 && P > 1)
    {
        size_t nKa = 0;
        float * Ka = gaussian_kernel(asigma, &nKa);
        fim_convn1(V, M, N, P, Ka, nKa, 2, 1);
        fim_free(Ka);
    }

    return;
}

/* Isotropic smoothing with Gaussian filter */
void fim_gsmooth(float * restrict V,
                 size_t M, size_t N, size_t P,
                 float sigma)
{
    fim_gsmooth_aniso(V, M, N, P, sigma, sigma);
}

void fimo_gsmooth(fimo * I, float sigma)
{
    fim_gsmooth(I->V, I->M, I->N, I->P, sigma);
    return;
}

static void cumsum_array(float * A, size_t N, size_t stride)
{
    for(size_t kk = 1; kk<N; kk++)
    {
        A[kk*stride] += A[(kk-1)*stride];
    }
}



void fim_cumsum(float * A, const size_t M, const size_t N, const int dim)
{
    assert(dim >= 0);
    assert(dim <= 1);
    if(dim == 0)
    {
#pragma omp parallel for
        for(size_t nn = 0; nn<N; nn++)
        {
            cumsum_array(A+M*nn, M, 1);
        }
    }
    if(dim == 1)
    {
#pragma omp parallel for
        for(size_t mm = 0; mm<M; mm++)
        {
            cumsum_array(A+mm, N, M);
        }
    }
}

/* Local sum in pM x pN windows, pads the input by (pM, pN), i.e. the
 * output is M+pM-1 x N+pN-1. Equivalent to convolution by a constant
 * array of size pM x pN */

float * fim_local_sum(const float * A,
                      size_t M, size_t N,
                      size_t pM, size_t pN)
{
    /* Work size */
    size_t wM = M + 2*pM;
    size_t wN = N + 2*pN;
    //printf("pM = %zu, pN = %zu wM = %zu wN = %zu\n", pM, pN, wM, wN);
    float * B = fim_zeros(wM*wN);

    /* Insert A, centered, into B */
    // TODO: fim_padarray
#pragma omp parallel for shared(A, B)
    for(size_t mm = 0; mm<M; mm++)
    {
        for(size_t nn = 0; nn<N; nn++)
        {
            B[mm+pM + (nn+pN)*wM] = A[mm + M*nn];
        }
    }
    //printf("B=\n");
    //fim_show(B, wM, wN, 1);
    fim_cumsum(B, wM, wN, 0);
    //c = s(1+m:end-1,:)-s(1:end-m-1,:);
    size_t cM = wM-pM-1;
    size_t cN = wN;
    //printf("cM = %zu, cN = %zu\n", cM, cN);
    float * C = fim_zeros(cM*cN);
    /* Difference along 0th dimension */
    //printf("B= cumsum dim 0\n");
    //fim_show(B, wM, wN, 1);
#pragma omp parallel for shared(C, B)
    for(size_t nn = 0; nn<cN; nn++)
    {
        float * B0 = B+wM*nn;
        for(size_t mm = 0; mm<cM; mm++)
        {
            C[nn*cM+mm] = B0[mm+pM] - B0[mm];
        }
    }
    //printf("C=\n");
    //fim_show(C, cM, cN, 1);
    fim_free(B);
    /* Second dimension */
    size_t dM = (wM-pM-1);
    size_t dN = (wN-pN-1);
    float * D = fim_zeros(dM*dN);
    fim_cumsum(C, cM, cN, 1);
    //printf("cumsum C, 1=\n");
    //fim_show(C, cM, cN, 1);
#pragma omp parallel for shared(D, C)
    for(size_t mm = 0; mm<dM; mm++)
    {
        float * C0 = C+mm;
        for(size_t nn = 0; nn<dN; nn++)
        {
            D[nn*dM+mm] = C0[(nn+pN)*cM] - C0[nn*cM];
        }
    }
    //printf("D=\n");
    //fim_show(D, dM, dN, 1);
    fim_free(C);
    return D;
}

static void fim_show(float * A, size_t M, size_t N, size_t P)
{
    for(size_t pp = 0; pp<P; pp++)
    {
        if(P > 1)
        {
            printf("z=%zu\n", pp);
        }
        for(size_t mm = 0; mm<M; mm++)
        {
            for(size_t nn = 0; nn<N; nn++)
            {
                printf("% .3f ", A[mm+nn*M+pp*M*N]);
            }
            printf("\n");
        }
    }
    return;
}

static void fim_show_int(int * A, size_t M, size_t N, size_t P)
{
    for(size_t pp = 0; pp<P; pp++)
    {
        if(P > 1)
        {
            printf("z=%zu\n", pp);
        }
        for(size_t mm = 0; mm<M; mm++)
        {
            for(size_t nn = 0; nn<N; nn++)
            {
                printf("% d ", A[mm+nn*M+pp*M*N]);
            }
            printf("\n");
        }
    }
    return;
}


float * fim_xcorr2(const float * T, const float * A,
                   const size_t M, const size_t N)
{


    /* Work size and size of output */
    size_t wM = 2*M-1;
    size_t wN = 2*N-1;

    fft_train(wM, wN, 1,
              1, 0,
              stdout);

    float * xA = fim_zeros(wM*wN);
    fim_insert(xA, wM, wN, 1,
               A, M, N, 1);

    float * temp = fim_zeros(M*N);
    //printf("Allocated %zu x %zu = %zu elements\n", M, N, M*N);
#pragma omp parallel for shared(temp, T)
    for(size_t kk = 0; kk < M*N; kk++)
    {
        temp[kk] = T[M*N-1-kk];
    }
    //printf("temp=\n");
    //fim_show(temp, M, N, 1);
    float * xT = fim_zeros(wM*wN);
    fim_insert(xT, wM, wN, 1,
               temp, M, N, 1);
    fim_free(temp);
    temp = NULL;

    fftwf_complex * fxT = fft(xA, wM, wN, 1);
    fim_free(xA);
    fftwf_complex * fxA = fft(xT, wM, wN, 1);
    fim_free(xT);
    float * C = fft_convolve_cc(fxT, fxA, wM, wN, 1);
    // TODO fft_convolve_cc_conj instead?
    //fim_tiff_write_float("C.tif", C, NULL, wM, wN, 1);
    fim_free(fxT);
    fim_free(fxA);

    //printf("C=\n");
    //fim_show(C, 2*M-1, 2*N-1, 1);

    float * LS = fim_local_sum(A, M, N,
                               M, N);
    //printf("LS=\n");
    //fim_show(LS, 2*M-1, 2*N-1, 1);
    float * numerator = fim_zeros(wM*wN);
    const float MN = M*N;
    const float sumT = fim_sum(T, M*N);
#pragma omp parallel for shared(C, LS, numerator)
    for(size_t kk = 0; kk<wM*wN; kk++)
    {
        numerator[kk] = (C[kk] - LS[kk]*sumT/MN);
    }
    //printf("Numerator = \n");
    //fim_show(numerator, wM, wN, 1);

    float * A2 = fim_zeros(M*N);
#pragma omp parallel for shared(A, A2)
    for(size_t kk = 0; kk<M*N; kk++)
    {
        A2[kk] = powf(A[kk], 2);
    }

    float * LS2 = fim_local_sum(A2, M, N,
                                M, N);
    fim_free(A2);
    float * denom_I = fim_zeros(wM*wN);

#pragma omp parallel for shared(LS2, denom_I)
    for(size_t kk = 0; kk<wM*wN; kk++)
    {
        float dLS = ( LS2[kk] - powf(LS[kk],2)/MN );
        float v = sqrt(dLS);
        v < 0 ? v = 0 : 0;
        denom_I[kk] = v;
    }
    fim_free(LS);
    fim_free(LS2);
    //printf("denom_I=\n");
    //fim_show(denom_I, wM, wN, 1);

    float denom_J = sqrt(M*N-1.0)*fim_std(T, M*N);
    //printf("std = %f\n", fim_std(T, M*N));
    //printf("denom_J = %f\n", denom_J);
#pragma omp parallel for shared(denom_I, C, numerator)
    for(size_t kk = 0; kk<wM*wN; kk++)
    {
        float d = denom_I[kk]*denom_J;
        if(d > 1e-9)
        {
            C[kk] = numerator[kk] / d;
        }
    }
    fim_free(numerator);
    fim_free(denom_I);
    return C;
}


float fim_std(const float * V, size_t N)
{

    float mean = fim_sum(V, N)/N;
    double s = 0;
#pragma omp parallel for reduction(+:s)
    for(size_t kk = 0; kk<N; kk++)
    {
        s += pow(V[kk]-mean, 2);
    }

    return (float) sqrt(s/ (double) (N - 1.0) );
}

/* Don't parallelize, called by fim_maxproj */
float fim_array_max(const float * A, size_t N, size_t stride)
{
    float m = A[0];
    for(size_t kk = 0; kk<N; kk++)
    {
        A[kk*stride] > m ? m = A[kk*stride] : 0;
    }
    return m;
}

float * fim_maxproj_ref(const float * V, size_t M, size_t N, size_t P)
{
    float * Pr = fim_zeros(M*N);
#pragma omp parallel for shared(Pr, V)
    for(size_t mm = 0; mm<M ; mm++)
    {
        for(size_t nn = 0; nn<N; nn++)
        {
            Pr[mm+M*nn] = fim_array_max(V+mm+M*nn, P, M*N);
        }
    }
    return Pr;
}

fimo * fimo_maxproj(const fimo * F)
{
    float * _M = fim_maxproj(F->V, F->M, F->N, F->P);
    fimo * M = malloc(sizeof(fimo));
    assert(M != NULL);
    M->V = _M;
    M->M = F->M;
    M->N = F->N;
    M->P = 1;
    return M;
}

fimo * fimo_sumproj(const fimo * F)
{
    if(F->P == 1)
    {
        return fimo_copy(F);
    }
    fimo * P = fimo_zeros(F->M, F->N, 1);
    for(size_t pp = 0; pp < F->P; pp++)
    {
        float * plane = F->V + F->M*F->N*pp;
        for(size_t kk = 0 ; kk < fimo_nel(P); kk++)
        {
            P->V[kk] += plane[kk];
        }
    }
    return P;
}

float * fim_maxproj(const float * V, size_t M, size_t N, size_t P)
{
    /* Initialize to first plane */
    float * Pr = fim_copy(V, M*N);

    /* Compare to remaining planes */
    for(size_t pp = 1; pp<P ; pp++)
    {
        const float * Vp = V + pp*M*N;
        for(size_t kk = 0; kk<M*N; kk++)
        {
            Pr[kk] < Vp[kk] ? Pr[kk] = Vp[kk] : 0;
        }
    }
    return Pr;
}


float * fim_sumproj(const float * V, size_t M, size_t N, size_t P)
{
    float * Pr = fim_zeros(M*N);
#pragma omp parallel for shared(Pr, V)
    for(size_t mm = 0; mm<M ; mm++)
    {
        for(size_t nn = 0; nn<N; nn++)
        {
            float sum = 0;
            for(size_t pp = 0; pp<P; pp++)
            {
                sum += V[mm+M*nn+pp*M*N];
            }
            Pr[mm+M*nn] = sum;
        }
    }
    return Pr;
}



float strel333_max(const float * I, size_t M, size_t N,
                   __attribute__((unused))  size_t P,
                   const float * strel)
{

    float max = 0; /* ok since we handle image data */
    for(int aa = 0 ; aa < 3; aa++)
    {
        for(int bb = 0 ; bb < 3; bb++)
        {
            for(int cc = 0 ; cc < 3; cc++)
            {
                if(strel[aa+3*bb+9*cc] == 1)
                {
                    float v = I[(aa-1) + (bb-1)*M + (cc-1)*M*N];
                    v > max ? max = v : 0;
                }
            }
        }
    }
    return max;
}

/** @brief Locate the max of a vector using 2nd deg poly2
 *
 * i = arg max y[i]
 * if i==0 return x[i]
 * if i==n-1 return x[i]
 * For anything in between, use a polynomial model
 * y = c_0 + c_1 x + c_2 x^2
 * based on x[i-1], x[i], x[i+1], y[i-1], y[i], y[i+1]
 */
static float locate_max_poly2(const float * x, const float * y, int n,
                              int use_log)
{
    int maxpos = 0;
    int max = y[0];
    for(int kk = 1; kk < n; kk++)
    {
        if(y[kk] > max)
        {
            max = y[kk];
            maxpos = kk;
        }
    }
    if(maxpos == 0)
    {
        return x[0];
    }
    if(maxpos == n-1)
    {
        return x[n-1];
    }

    /* Find coefficients */
    float a = x[maxpos-1];
    float b = x[maxpos];
    float c = x[maxpos+1];
    if(use_log)
    {
        a = log(a);
        b = log(b);
        c = log(c);
    }
    float a2 = pow(a, 2);
    float b2 = pow(b, 2);
    float c2 = pow(c, 2);

    /* c0 not needed */
    float C1 = y[maxpos-1]*((-b - c)/(a2 - a*b - a*c + b*c))
        +y[maxpos]*((a + c)/(a*b - a*c - b2 + b*c))
        +y[maxpos+1]*((-a - b)/(a*b - a*c - b*c + c2));
    float C2 = y[maxpos-1]*(1.0/(a2 - a*b - a*c + b*c))
        + y[maxpos]*(-1.0/(a*b - a*c - b2 + b*c))
        + y[maxpos+1]*(1.0/(a*b - a*c - b*c + c2));
    /* d/dx y(x) = 0 -> 0 = c1 + 2*c2*x, x=-c1/(2*c2) */
    if(use_log)
    {
        return exp(-C1 / (2.0 * C2));
    }
    return -C1 / (2.0 * C2);
}

ftab_t * fim_lmax_multiscale(float ** II, float * scales, size_t nscales,
                             size_t M, size_t N, size_t P)
{
    int ncol = 4 + 1 + nscales;
    ftab_t * T = ftab_new(ncol);
    ftab_set_colname(T, 0, "x");
    ftab_set_colname(T, 1, "y");
    ftab_set_colname(T, 2, "z");
    ftab_set_colname(T, 3, "value");
    ftab_set_colname(T, 4, "LoG_scale");
    for(size_t kk = 0; kk<nscales; kk++)
    {
        char * cname = calloc(128, 1);
        assert(cname != NULL);
        sprintf(cname, "LoG_%f", scales[kk]);
        ftab_set_colname(T, 5+kk, cname);
        free(cname);
    }

#pragma omp parallel
    {
        float * row = calloc(ncol, sizeof(float));
        assert(row != NULL);
        float * log_values = calloc(ncol, sizeof(float));
        assert(log_values != NULL);

        /* Detect if a pixel is a local maxima
         * at any scale, i.e. if any of the nscales pixels
         * is a local maxima in the 4D neighbourhood
         * we skip the border pixels
         **/
        float * strel = calloc(27,sizeof(float));
        assert(strel != NULL);
        for(int kk = 0; kk<27; kk++)
        {
            strel[kk] = 1;
        }
        strel[13] = 0;

        /* For each scale, save the local max */
        float * local_max = calloc(nscales, sizeof(float));
        assert(local_max != NULL);
        /* Indicator if the central pixel is the largest for the given scale */
        int * is_local_max = calloc(nscales, sizeof(int));
        assert(is_local_max != NULL);
        int * is_scale_max = calloc(nscales, sizeof(int));
        assert(is_scale_max != NULL);

#pragma omp for
        for(size_t pp = 1; pp < P-1; pp++)
        {
            for(size_t nn = 1; nn+1 < N; nn++)
            {
                for(size_t mm = 1; mm+1 < M; mm++)
                {
                    size_t idx = pp*M*N + nn*M + mm;

                    /* For each scale, determine if the pixel is a local maxima
                     * also save the maxima of the region */
                    for(size_t ss = 0; ss < nscales; ss++)
                    {
                        local_max[ss] = strel333_max(II[ss]+idx, M, N, P, strel);
                        is_local_max[ss] = II[ss][idx] > local_max[ss];
                        II[ss][idx] > local_max[ss] ? local_max[ss] = II[ss][idx] : 0 ;
                    }

                    /* To be a maxima in the scale space, the pixel
                     * has to be a maxima also compared to the adjacent
                     * scales. */
                    for(int ss = 0; ss < (int) nscales; ss++)
                    {
                        is_scale_max[ss] = is_local_max[ss];

                        if(ss > 0) /* Compare with the previous scale */
                        {
                            if(local_max[ss-1] > local_max[ss])
                            {
                                is_scale_max[ss] = 0;
                            }
                        }
                        if(ss + 1 < (int) nscales)
                        {
                            if(local_max[ss+1] > local_max[ss])
                            {
                                is_scale_max[ss] = 0;
                            }
                        }
                    }

                    /* Save the coordinates of the maximas */
                    for(size_t ss = 0; ss < nscales; ss++)
                    {
                        if(is_scale_max[ss])
                        {
                            /* Pos is s a local maxima */
                            row[0] = mm; row[1] = nn; row[2] = pp;
                            float max_value = II[0][idx];
                            for(size_t kk = 0; kk < nscales; kk++)
                            {
                                float value = II[kk][idx];
                                log_values[kk] = value;
                                row[5+kk] = value;
                                value > max_value ? max_value = value : 0;
                            }
                            row[3] = max_value;
                            int use_log = 1;
                            row[4] = locate_max_poly2(scales, log_values, nscales, use_log);
#pragma omp critical
                            ftab_insert(T, row);
                        }
                    }
                }
            }
        }

        free(is_local_max);
        free(local_max);
        free(strel);
        free(row);
        free(log_values);
        free(is_scale_max);
    }
    return T;
}

static float max_float(float a, float b)
{
    if(a > b) {
        return a;
    }
    return b;
}

static int lmax_2d(const float * I, size_t stride)
{
    float other = I[-1];
    other = max_float(other, I[1]);
    other = max_float(other, I[ stride]);
    other = max_float(other, I[-stride]);
    other = max_float(other, I[1+stride]);
    other = max_float(other, I[1-stride]);
    other = max_float(other, I[-1+stride]);
    other = max_float(other, I[-1-stride]);

    return I[0] > other;
}

ftab_t * fim_lmax_2d(const float * I, size_t M, size_t N, size_t P)
{
    if(P != 1)
    {
        printf("Warning: fim_lmax_2d called with a 3D image\n"
               "Will only use the first plane\n");
    }
    ftab_t * T = ftab_new(4);
    ftab_set_colname(T, 0, "x");
    ftab_set_colname(T, 1, "y");
    ftab_set_colname(T, 2, "z");
    ftab_set_colname(T, 3, "value");

    // TODO
    /* 3x3x3 structuring element with 26-connectivity */
    float * strel = malloc(27*sizeof(float));
    assert(strel != NULL);
    for(int kk = 0; kk<27; kk++)
    {
        strel[kk] = 1;
    }
    strel[13] = 0;

    assert(P > 0);
#pragma omp parallel for
    for(size_t nn = 1; nn < N-1; nn++)
    {
        for(size_t mm = 1; mm+1 < M; mm++)
        {
            size_t pos = mm + nn*M;
            if(lmax_2d(I + pos, M))
            {
                /* Pos is s a local maxima */
                float row[4] = {mm, nn, 0.0, I[pos]};
#pragma omp critical
                ftab_insert(T, row);
            }
        }
    }
    free(strel);
    return T;
}

ftab_t * fim_lmax(const float * I, size_t M, size_t N, size_t P)
{
    if(P == 1)
    {
        return fim_lmax_2d(I, M, N, P);
    }

    ftab_t * T = ftab_new(4);
    ftab_set_colname(T, 0, "x");
    ftab_set_colname(T, 1, "y");
    ftab_set_colname(T, 2, "z");
    ftab_set_colname(T, 3, "value");

    /* 3x3x3 structuring element with 26-connectivity */
    float * strel = malloc(27*sizeof(float));
    assert(strel != NULL);
    for(int kk = 0; kk<27; kk++)
    {
        strel[kk] = 1;
    }
    strel[13] = 0;

    assert(P > 0);
#pragma omp parallel for
    for(size_t pp = 1; pp < P-1; pp++)
    {
        for(size_t nn = 1; nn+1 < N; nn++)
        {
            for(size_t mm = 1; mm+1 < M; mm++)
            {
                size_t pos = mm + nn*M + pp*M*N;
                if(I[pos] > strel333_max(I + pos, M, N, P, strel))
                {
                    /* Pos is s a local maxima */
                    float row[4] = {mm, nn, pp, I[pos]};
#pragma omp critical
                    ftab_insert(T, row);
                }
            }
        }
    }

    fim_free(strel);
    return T;
}


fim_histogram_t * fim_histogram(const float * Im, size_t N)
{
    if(N < 2)
    {
        return NULL;
    }
    float min = fim_min(Im, N);
    float max = fim_max(Im, N);
    size_t nbin = pow(2, 16)+1;
    float delta = (max-min) / ((float) nbin);
    float left = min - 0.5*delta-1e-6;
    float right = max + 0.5*delta+1e-6;
    //delta = (right-left)/((float) nbin);
    //printf("Histogram: #=%zu [%f, %f], delta=%f\n", nbin, left, right, delta);
    fim_histogram_t * H = malloc(sizeof(fim_histogram_t));
    assert(H != NULL);
    H->left = left;
    H->right = right;
    H->nbin = nbin;
    H->C = calloc(nbin, sizeof(double));
    assert(H->C != NULL);


    for(size_t kk = 0; kk<N; kk++)
    {
        float v = Im[kk]; /* Value */
        //        printf("v=%f\n", v);
        float fbin = (v-left)/(right-left)*((float) (nbin));
        assert(fbin < nbin);
        assert(fbin >= 0);
        //printf("v=%f, bin=%f\n", v, fbin);
        size_t bin = round(fbin);
        if(bin >= (size_t) nbin)
        {
            bin = nbin-1;
        }
        H->C[bin]++;
    }
    return H;
}

void fim_histogram_log(fim_histogram_t * H)
{
    for(size_t kk = 0; kk<H->nbin; kk++)
    {
        H->C[kk] = log(1+H->C[kk]);
    }
}

void fim_histogram_free(fim_histogram_t * H)
{
    if(H == NULL)
    {
        return;
    }
    if(H->C != NULL)
    {
        fim_free(H->C);
    }
    free(H);
    return;
}

size_t otsu(double * C, size_t N)
{
    double total = fim_sum_double(C, N);
    double sumB = 0;
    double wB = 0;
    double maximum = 0;
    double sum1 = 0;
    for(float ii = 0; ii<N; ii++)
    {
        sum1 += ii*C[(int) ii];
    }
    float best_level = 0;
    for(size_t ii = 0 ; ii < N; ii++)
    {
        double wF = total - wB;

        if(wB > 0 && wF > 0)
        {
            double mF = (sum1 - sumB) / wF;
            double val = wB * wF * ((sumB / wB) - mF) * ((sumB / wB) - mF);

            if ( val >= maximum )
            {
                maximum = val;
                best_level = ii;
            }
        }
        wB += C[(int) ii];
        sumB += ii*C[(int) ii];
    }

    return best_level;
}

float fim_histogram_percentile(const fim_histogram_t * H, const float p)
{
    if(p <= 0)
    {
        return H->left;
    }
    if(p >= 1)
    {
        return H->right;
    }

    float p0 = 0;
    size_t pos = 0;
    float sum = fim_sum_double(H->C, H->nbin);

    while(p0 < p*sum)
    {
        p0 += H->C[pos++];
    }

    float rpos = (float) pos / H->nbin; /* Relative value in [0,1] */
    float th = H->left + (float) rpos*(H->right-H->left);
    return th;
}

float fim_histogram_otsu(fim_histogram_t * H)
{
    assert(H != NULL);
    double level = otsu(H->C, H->nbin);
    double slevel = H->left + level/(pow(2,16)-1)*(H->right-H->left);
    return (float) slevel;
}

float * fim_otsu(float * Im, size_t M, size_t N)
{
    //float mean = fim_sum(Im, M*N)/( (float) M*N);
    //printf("Mean=%f\n", mean);
    fim_histogram_t * H = fim_histogram(Im, M*N);
    assert(H != NULL);
    float th = fim_histogram_otsu(H);
    //    printf("Threshold = %f\n", th);
    fim_histogram_free(H);
    float * B = malloc(M*N*1*sizeof(float));
    assert(B != NULL);
    for(size_t kk = 0; kk < M*N; kk++)
    {
        if( Im[kk] > th )
        {
            B[kk] = 1;
        } else {
            B[kk] = 0;
        }
    }
    return B;
}

float * fim_remove_small(const float * im,
                         size_t M, size_t N,
                         float min_pixels)
{
    int * L = fim_conncomp6(im, M, N);
    int ncomp = 0;
    for(size_t kk = 0; kk<M*N; kk++)
    {
        L[kk] > ncomp ? ncomp = L[kk] : 0;
    }
    printf("%d objects\n", ncomp);

    /* Histogram on the number of pixels per label */
    uint32_t * H = fim_malloc( (ncomp+1)*sizeof(uint32_t));
    assert(H != NULL);
    for(size_t kk = 0; kk<M*N; kk++)
    {
        H[L[kk]]++;
    }

    /* Set  */
    size_t npix = 0;
    float * out = fim_malloc(M*N*sizeof(float));
    assert(out != NULL);
    for(size_t kk = 0; kk<M*N; kk++)
    {
        out[kk] = im[kk];
        if((L[kk] > 0) & (out[kk] > 0))
        {
            if(H[L[kk]] < min_pixels)
            {
                out[kk] = 0;
                npix++;
            }
        }
    }
    printf("Cleared %zu pixels\n", npix);
    fim_free(H);
    fim_free(L);
    return out;
}

float * fim_fill_holes(const float * im, size_t M, size_t N, float max_size)
{
    /* Any region in !im that is not connected to the boundary
     * is a hole
     */
    float * nim = fim_malloc(sizeof(float)*M*N);
    if(nim == NULL)
    {
        exit(EXIT_FAILURE);
    }
    for(size_t kk = 0 ; kk<M*N; kk++)
    {
        nim[kk] = (im[kk] == 0);
    }
    int * L = fim_conncomp6(nim, M, N);
    fim_free(nim);

    /* Count the number of components
     * or actually the largest label used
     */

    int ncomp = 0;
    for(size_t kk = 0; kk<M*N; kk++)
    {
        L[kk] > ncomp ? ncomp = L[kk] : 0;
    }
    printf("%d potential holes\n", ncomp);

    /* Find regions connected to the boundary */
    int * boundary = calloc(ncomp+1, sizeof(int));
    assert(boundary != NULL);

    for(size_t kk = 0; kk<M; kk++)
    {
        boundary[L[kk]] = 1;
        boundary[L[kk + (N-1)*M]] = 1;
    }
    for(size_t kk = 0; kk<N; kk++)
    {
        boundary[L[kk*M]] = 1;
        boundary[L[kk*M + (M-1)]] = 1;
    }
    /* Clear all labels connected to boundary */
    for(size_t kk = 0; kk<M*N; kk++)
    {
        if(boundary[L[kk]] == 1)
        {
            L[kk] = 0;
        }
    }
    fim_free(boundary);

    /* Histogram on the number of pixels per label */
    uint32_t * H = fim_malloc((ncomp+1)*sizeof(uint32_t));
    assert(H != NULL);
    for(size_t kk = 0; kk<M*N; kk++)
    {
        H[L[kk]]++;
    }

    /* Finally, return a copy of the input where
     * the regions corresponding to labeled component
     * in the inverted image that are not connected to the
     * boundary are set to zero. */
    size_t nfilled = 0;
    float * out = fim_malloc(M*N*sizeof(float));
    assert(out != NULL);
    for(size_t kk = 0; kk<M*N; kk++)
    {
        out[kk] = im[kk];
        if(im[kk] == 0)
        {
            if( (L[kk] > 0) & (H[L[kk]] < max_size) )
            {
                out[kk] = 1;
                nfilled++;
            }
        }
    }
    printf("Filled %zu pixels\n", nfilled);
    fim_free(H);
    fim_free(L);
    return out;
}

int * fim_conncomp6(const float * im, size_t M, size_t N)
{
    int * lab = fim_malloc(M*N*sizeof(float));
    assert(lab != NULL);
    int label = 0;

    //    printf("First pass\n");
    /* First label along the first dimension */
    for(size_t nn = 0; nn<N; nn++)
    {
        int currlabel = 0;
        const float * Aim = im + nn*M;
        int * Alab = lab + nn*M;
        for(size_t mm = 0; mm<M; mm++)
        {
            if(currlabel == 0)
            {
                if(Aim[mm] > 0)
                {
                    Alab[mm] = ++label;
                    currlabel = label;
                }
            } else {
                if(Aim[mm] > 0)
                {
                    Alab[mm] = label;
                } else {
                    currlabel = 0;
                }
            }
        }
    }

    //    fim_show_int(lab, M, N, 1);
    label++;
    //    printf("Second pass, set up equivalences\n");
    int * E = fim_malloc(label*sizeof(int));
    assert(E != NULL);
    for(size_t kk = 0; kk< (size_t) label; kk++)
    {
        E[kk] = kk;
    }
    int changed = 1;
    while(changed > 0)
    {
        //   printf("."); fflush(stdout);
        changed = 0;
        for(size_t mm = 0; mm<M; mm++)
        {
            int * Alab = lab + mm;
            for(size_t nn = 0; nn+1<N; nn++)
            {
                int a = Alab[nn*M];
                if(a > 0)
                {
                    int b = Alab[nn*M+M];
                    if(b > 0)
                    {
                        if(E[a] < E[b])
                        {
                            E[b] = E[a];
                            changed++;
                        }
                        if(E[a] > E[b])
                        {
                            E[a] = E[b];
                            changed++;
                        }
                    }
                }
            }
        }
        //printf("%d\n", changed);
    }
    //printf("\n");


    //for(size_t kk = 0; kk<M*N; kk++)
    //{
    //lab[kk] = E[lab[kk]];
    //}


    //    printf("Third pass, set pixel values\n");
    int * E2 = fim_malloc(label*sizeof(int));
    assert(E2 != NULL);
    int newlabel = 1;
    for(size_t kk = 0; kk<M*N; kk++)
    {
        if(lab[kk] > 0)
        {
            int e = E[lab[kk]];
            int e2 = E2[e];
            if(e2 == 0)
            {
                e2 = newlabel;
                E2[e] = newlabel++;
            }
            lab[kk] = e2;
        }
    }
    fim_free(E2); E2 = NULL;

    fim_free(E); E = NULL;
    return lab;
}

/* 1D convolution along 1 dimension in a 3D image */
int fim_convn1(float * restrict V,
               size_t M, size_t N, size_t P, // image size
               const float * K, size_t nK, // Kernel
               int dim, const int normalized)
{
    if(dim < 0 || dim > 2)
    {
        return EXIT_FAILURE;
    }

    /* Temporary storage/buffer for conv1_vector */
    size_t nBuff = max_size_t(M, max_size_t(N, P));

#pragma omp parallel
    {

        float * buff = fim_malloc(nBuff*sizeof(float));
        assert(buff != NULL);

        if(dim == 0)
        {
#pragma omp for
            for(size_t pp = 0; pp < P; pp++)
            {
                for(size_t nn = 0; nn < N; nn++)
                {
                    fim_conv1_vector(V+pp*(M*N)+nn*M, 1, buff, M, K, nK, normalized);
                }
            }
        }

        if(dim == 1)
        {
#pragma omp for
            for(size_t pp = 0; pp<P; pp++)
            {
                for(size_t mm = 0; mm<M; mm++)
                {
                    fim_conv1_vector(V + pp*(M*N) + mm, M, buff, N, K, nK, normalized);
                }
            }
        }

        if(dim == 2)
        {
#pragma omp for
            for(size_t nn = 0; nn<N; nn++)
            {
                for(size_t mm = 0; mm<M; mm++)
                {
                    fim_conv1_vector(V+mm+M*nn, M*N, buff, P, K, nK, normalized);
                }
            }
        }

        fim_free(buff);
    }
    return EXIT_SUCCESS;
}

/** Separable convolution.
 *
 * Convolve V by K1 in the 1st dimension, K2
 * in the 2nd dimension and K3 in the third dimension.  This version
 * use fim_shiftdim to reduce the computational load.
 */
float * conv1_3(const float * restrict V, size_t M, size_t N, size_t P,
                const float * K1, size_t nK1,
                const float * K2, size_t nK2,
                const float * K3, size_t nK3)
{
    const int dim = 0;
    const int norm = 0;

    fimo * F = fim_image_from_array(V, M, N, P);

    fim_convn1(F->V, F->M, F->N, F->P, K1, nK1, dim, norm);
    fimo * F2 = fim_shiftdim(F);
    fim_free(F->V);
    fim_free(F);
    fim_convn1(F2->V, F2->M, F2->N, F2->P, K2, nK2, dim, norm);
    fimo * F3 = fim_shiftdim(F2);
    fim_free(F2->V);
    fim_free(F2);
    fim_convn1(F3->V, F3->M, F3->N, F3->P, K3, nK3, dim, norm);
    fimo * F4 = fim_shiftdim(F3);
    fim_free(F3->V);
    fim_free(F3);
    float * out = F4->V;
    fim_free(F4);

    return out;
}

float * fim_LoG_S(const float * V0, const size_t M, const size_t N, const size_t P0,
                  const float sigmaxy, const float sigmaz)
{

    /* Set up filters */
    /* Lateral filters */
    size_t nlG = 0;
    float * lG = gaussian_kernel(sigmaxy, &nlG);
    size_t nl2;
    float * l2 = gaussian_kernel_d2(sigmaxy, &nl2);
    flip_sign(l2, nl2);
    /* Axial filters */
    size_t naG = 0;
    float * aG = gaussian_kernel(sigmaz,  &naG);
    size_t na2;
    float * a2 = gaussian_kernel_d2(sigmaz,  &na2);
    flip_sign(a2, na2);

    /* Padding */
    int apad = (naG-1)/2;
    int lpad = (nlG-1)/2;
    if(fim_verbose > 1)
    {
        printf("fim_LoG_S, apad: %d, lpad: %d\n", apad, lpad);
    }

    /* Pad in Z-direction */
    size_t P = P0 + 2*apad;
    float * V = fim_malloc(M*N*P*sizeof(float));

    for(size_t kk = 0; kk< (size_t) apad; kk++)
    {
        memcpy(V+kk*M*N, V0, M*N*sizeof(float));
    }
    memcpy(V+apad*M*N, V0, M*N*P0*sizeof(float));
    for(size_t kk = apad+P0; kk<P; kk++)
    {
        memcpy(V+kk*M*N, V0+M*N*(P0-1), M*N*sizeof(float));
    }

    float * GGL = conv1_3(V, M, N, P,
                          lG, nlG, lG, nlG, a2, na2);
    float * GLG = conv1_3(V, M, N, P,
                          lG, nlG, l2, nl2, aG, naG);
    float * LGG = conv1_3(V, M, N, P,
                          l2, nl2, lG, nlG, aG, naG);
    fim_free(V);
    float * LoG = GGL;
    for(size_t kk = 0; kk<M*N*P; kk++)
    {
        LoG[kk] += GLG[kk] + LGG[kk];
    }
    fim_free(GLG);
    fim_free(LGG);

    /* Set border to 0 */
    size_t pos = 0;
    for(int pp = 0; pp< (int) P; pp++)
    {
        int uz = 1;
        if(pp<apad || pp+apad >= (int) P)
        {
            uz = 0;
        }
        for(int nn = 0; nn< (int) N; nn++)
        {
            int uy = 1;
            if(nn<lpad || nn+lpad >= (int) N)
            {
                uy = 0;
            }
            for(int mm = 0; mm< (int) M; mm++)
            {
                int ux = 1;
                if(mm<lpad || mm+lpad >= (int) M)
                {
                    ux = 0;
                }
                if(ux*uy*uz != 1)
                {
                    LoG[pos] = 0;
                }
                pos++;
            }
        }
    }
    fim_free(aG);
    fim_free(lG);
    fim_free(a2);
    fim_free(l2);

    /* Unpad in the axial direction */
    float * uLoG = fim_malloc(M*N*P0*sizeof(float));

    memcpy(uLoG, LoG+M*N*apad, M*N*P0*sizeof(float));
    fim_free(LoG);
    return uLoG;
}

float *
fim_LoG_S2(const float * V0,
           const size_t M, const size_t N, const size_t P,
           const float sigmaxy, const float sigmaz)
{

    fim_boundary_condition bc = FIM_BC_SYMMETRIC_MIRROR;

    /* Set up filters */
    /* Lateral filters */
    size_t nlG = 0;
    float * _lG = gaussian_kernel(sigmaxy, &nlG);
    fimo * lG = fim_wrap_array(_lG, nlG, 1, 1);
    size_t nl2;
    float * _l2 = gaussian_kernel_d2(sigmaxy, &nl2);
    flip_sign(_l2, nl2);
    fimo * l2 = fim_wrap_array(_l2, nl2, 1, 1);

    /* 2D */
    if(P == 1)
    {
        assert(fimo_nel(lG) > 1);
        assert(fimo_nel(l2) > 1);

        printf("2D path\n");
        /** Gaussian, Laplacian **/
        fimo * GI = fim_image_from_array(V0, M, N, P);
        fimo_conv1_x(GI, lG, bc);
        fimo * GL = fim_shiftdim2(GI);
        assert(GL->P == 1);
        fim_free(GI);
        fimo_conv1_x(GL, l2, bc);

        /** Laplacian, Gaussian **/
        fimo * LI = fim_image_from_array(V0, M, N, P);
        fimo_conv1_x(LI, l2, bc);
        assert(LI->P == 1);
        fimo * LG = fim_shiftdim2(LI);
        fim_free(LI);
        fimo_conv1_x(LG, lG, bc);

        /* Free filters */
        fimo_free(lG);
        fimo_free(l2);

        /* Add together filter responses */
        fimo * LoG = GL;
        GL = NULL;
        fimo_add(LoG, LG);
        fimo_free(LG);

        /* Shift back */
        fimo * result = fim_shiftdim2(LoG);

        fimo_free(LoG);
        /* And we are done */

        float * pLoG = result->V;

        result->V = NULL;
        fimo_free(result);
        return pLoG;
    }

    /* 3D */
    if(P > 1)
    {
        /* Axial filters */
        fimo * aG = NULL;
        fimo * a2 = NULL;

        {
            size_t naG = 0;
            size_t na2 = 0;
            float * _aG = gaussian_kernel(sigmaz,  &naG);
            aG = fim_wrap_array(_aG, naG, 1, 1);
            float * _a2 = gaussian_kernel_d2(sigmaz,  &na2);
            flip_sign(_a2, na2);
            a2 = fim_wrap_array(_a2, na2, 1, 1);
        }


        /** First dimension -> GII, LII */
        fimo * GII = fim_image_from_array(V0, M, N, P);
        fimo_conv1_x(GII, lG, bc);

        fimo * LII = fim_image_from_array(V0, M, N, P);
        fimo_conv1_x(LII, l2, bc);

        /** 2nd dimension -> GGI, GLI, LGI */
        /* Prepare buffers */
        fimo * GGI = fim_shiftdim(GII);
        fimo_free(GII);
        fimo * GLI = fimo_copy(GGI);
        fimo * LGI = fim_shiftdim(LII);
        fimo_free(LII);
        /* Apply filters */
        fimo_conv1_x(GGI, lG, bc);
        fimo_conv1_x(GLI, l2, bc);
        fimo_conv1_x(LGI, lG, bc);

        /** 3rd dimension -> GGL, GLG, LGG */
        fimo * GGL = fim_shiftdim(GGI);
        fimo_free(GGI);
        fimo * GLG = fim_shiftdim(GLI);
        fimo_free(GLI);
        fimo * LGG = fim_shiftdim(LGI);
        fimo_free(LGI);

        fimo_conv1_x(GGL, a2, bc);
        fimo_conv1_x(GLG, aG, bc);
        fimo_conv1_x(LGG, aG, bc);

        /* Free the filters */
        fimo_free(lG);
        fimo_free(l2);
        fimo_free(aG);
        fimo_free(a2);

        /** Merge results */
        fimo * LoG = GGL;
        GGL = NULL;
        fimo_add(LoG, GLG);
        fimo_free(GLG);
        fimo_add(LoG, LGG);
        fimo_free(LGG);

        /** Shift back to original shape */
        fimo * result = fim_shiftdim(LoG);
        fimo_free(LoG);

        float * pLoG = result->V;
        result->V = NULL;
        fimo_free(result);
        return pLoG;
    }

    return NULL; // We should not reach this
}

/* Determinant of Hessian for 2D or 3D images */
float *
fim_DoH(const float * V0,
        const size_t M, const size_t N, const size_t P,
        const float sigmaxy, const float sigmaz)
{
    fim_boundary_condition bc = FIM_BC_SYMMETRIC_MIRROR;

    /* Set up filters for the lateral dimensions */
    fimo * lG = NULL; /* Gaussian */
    fimo * l1 = NULL; /* First derivative */
    fimo * l2 = NULL; /* 2nd derivative */
    {
        size_t nlG = 0;
        float * _lG = gaussian_kernel(sigmaxy, &nlG);
        assert(nlG > 0);
        lG = fim_wrap_array(_lG, nlG, 1, 1);

        size_t nl1;
        float * _l1 = gaussian_kernel_d1(sigmaxy, &nl1);
        l1 = fim_wrap_array(_l1, nl1, 1, 1);

        size_t nl2;
        float * _l2 = gaussian_kernel_d2(sigmaxy, &nl2);
        assert(nl2 > 0);
        l2 = fim_wrap_array(_l2, nl2, 1, 1);
    }

    if(P == 1)
    {
        fimo * GI = fim_image_from_array(V0, M, N, P);
        fimo_conv1_x(GI, lG, bc);

        fimo * DI = fim_image_from_array(V0, M, N, P);
        fimo_conv1_x(DI, l1, bc);

        fimo * LI = fim_image_from_array(V0, M, N, P);
        fimo_conv1_x(LI, l2, bc);

        /* Shifted one dimension */
        fimo * sGL = fim_shiftdim(GI);
        fimo_free(GI);
        fimo_conv1_x(sGL, l2, bc);

        fimo * sDD = fim_shiftdim(DI);
        fimo_free(DI);
        fimo_conv1_x(sDD, l1, bc);

        fimo * sLG = fim_shiftdim(LI);
        fimo_free(LI);
        fimo_conv1_x(sLG, lG, bc);

        /* Free kernels */
        fimo_free(lG);
        fimo_free(l1);
        fimo_free(l2);

        /* Compute DoH (still shifted) */
#pragma omp parallel for
        for(size_t kk = 0; kk < M*N; kk++)
        {
            sGL->V[kk] = sGL->V[kk]*sLG->V[kk] - pow(sDD->V[kk], 2.0);
        }
        fimo_free(sLG);
        fimo_free(sDD);

        /* Shift back */
        fimo * _DoH = fim_shiftdim(sGL);
        fimo_free(sGL);
        float * DoH = _DoH->V;
        _DoH->V = NULL;
        fimo_free(_DoH);
        return DoH;
    }

    if(P > 1)
    {
        printf("3D DOH not implemented. %f\n", sigmaz);
        assert(0);
    }

    return NULL;
}

/* Parital derivative in dimension dim */
fimo * fimo_partial(const fimo * F, const int dim, const float sigma)
{

    if(F->P != 1)
    {
        fprintf(stderr, "fimo_partial only supports 2D images (please fix me)\n");
        exit(EXIT_FAILURE);
    }

    size_t nG = 0;
    float * G = gaussian_kernel(sigma, &nG);
    size_t nD1;
    float * D1 = gaussian_kernel_d1(sigma, &nD1);


    fimo * D = fimo_copy(F);
    if(dim == 0)
    {
        fim_convn1(D->V, D->M, D->N, D->P, D1, nD1, 0, 0);
        fim_convn1(D->V, D->M, D->N, D->P, G, nG, 1, 0);
        fim_free(D1);
        fim_free(G);
        return D;
    }
    if(dim == 1)
    {
        fim_convn1(D->V, D->M, D->N, D->P, G, nG, 0, 0);
        fim_convn1(D->V, D->M, D->N, D->P, D1, nD1, 1, 0);
        fim_free(D1);
        fim_free(G);
        return D;
    }

    fprintf(stderr, "Something went wrong in fimo_partial dim=%d\n", dim);
    exit(EXIT_FAILURE);


}

float * fim_LoG(const float * V, const size_t M, const size_t N, const size_t P,
                const float sigmaxy, const float sigmaz)
{

    /* Set up filters */
    /* Lateral filters */
    size_t nlG = 0;
    float * lG = gaussian_kernel(sigmaxy, &nlG);
    size_t nl2;
    float * l2 = gaussian_kernel_d2(sigmaxy, &nl2);
    flip_sign(l2, nl2);
    /* Axial filters */
    size_t naG = 0;
    float * aG = gaussian_kernel(sigmaz,  &naG);
    size_t na2;
    float * a2 = gaussian_kernel_d2(sigmaz,  &na2);
    flip_sign(a2, na2);

    /* 1st dimension */
    float * LoG = NULL;
    {
        float * GGL = fim_copy(V, M*N*P);
        fim_convn1(GGL, M, N, P, lG, nlG, 0, 0);
        fim_convn1(GGL, M, N, P, lG, nlG, 1, 0);
        fim_convn1(GGL, M, N, P, a2, na2, 2, 0);
        LoG = GGL;
    }

    /* 2nd dimension */
    {
        float * GLG = fim_copy(V, M*N*P);
        fim_convn1(GLG, M, N, P, lG, nlG, 0, 0);
        fim_convn1(GLG, M, N, P, l2, nl2, 1, 0);
        fim_convn1(GLG, M, N, P, aG, naG, 2, 0);
        fim_add(LoG, GLG, M*N*P);
        fim_free(GLG);
    }

    /* 3rd dimension */
    {
        float * LGG = fim_copy(V, M*N*P);
        fim_convn1(LGG, M, N, P, l2, nl2, 0, 0);
        fim_convn1(LGG, M, N, P, lG, nlG, 1, 0);
        fim_convn1(LGG, M, N, P, aG, naG, 2, 0);
        fim_add(LoG, LGG, M*N*P);
        fim_free(LGG);
    }

    /* Free the kernels */
    fim_free(lG);
    fim_free(l2);
    fim_free(aG);
    fim_free(a2);


    /* Set border to 0 */
    int apad = (naG-1)/2;
    int lpad = (nlG-1)/2;
    size_t pos = 0;
    for(int pp = 0; pp< (int) P; pp++)
    {
        int uz = 1;
        if(pp<apad || pp+apad >= (int) P)
        {
            uz = 0;
        }
        for(int nn = 0; nn< (int) N; nn++)
        {
            int uy = 1;
            if(nn<lpad || nn+lpad >= (int) N)
            {
                uy = 0;
            }
            for(int mm = 0; mm< (int) M; mm++)
            {
                int ux = 1;
                if(mm<lpad || mm+lpad >= (int) M)
                {
                    ux = 0;
                }
                if(ux*uy*uz != 1)
                {
                    LoG[pos] = 0;
                }
                pos++;
            }
        }
    }

    return LoG;
}

double * fim_get_line_double(fimo * I,
                             int x, int y, int z,
                             int dim, int nPix)
{
    const float * V = I->V;
    int M = I->M;
    int N = I->N;
    int P = I->P;
    int pos = z*M*N + y*M + x;
    double * L = fim_malloc(nPix*sizeof(double));
    assert(L != NULL);
    if(dim == 0)
    {
        for(int kk = 0; kk<nPix; kk++)
        {
            int d = kk-(nPix-1)/2;
            if(x + d < 0 || x + d >= M)
                continue;
            L[kk] = V[pos + d];
        }
    }
    if(dim == 1)
    {
        for(int kk = 0; kk<nPix; kk++)
        {
            int d = kk-(nPix-1)/2;
            if(y + d < 0 || y + d >= N)
                continue;
            L[kk] = V[pos + d*M];
        }
    }
    if(dim == 2)
    {
        for(int kk = 0; kk<nPix; kk++)
        {
            int d = kk-(nPix-1)/2;
            if(z + d < 0 || z + d >= P)
                continue;
            L[kk] = V[pos + d*M*N];
        }
    }
    return L;
}

fimo * fim_shiftdim2(const fimo * restrict I)
{
    const float * V = I->V;
    const size_t M = I->M;
    const size_t N = I->N;
    const size_t P = I->P;
    assert(P == 1);
    /* Output image */
    fimo * O = malloc(sizeof(fimo));
    assert(O != NULL);
    O->V = fim_malloc(M*N*P*sizeof(float));
    O->M = N;
    O->N = M;
    O->P = P;
    float * S = O->V;
    const size_t blocksize = 2*64;

    /* Blocks are of size blocksize x blocksize in M and N */
    for(size_t bm = 0; bm<M; bm = bm+blocksize)
    {
        size_t bm_end = bm+blocksize;
        bm_end > M ? bm_end = M : 0;

        const size_t bm_end_c = bm_end;

        for(size_t bn = 0; bn<N; bn = bn+blocksize)
        {
            size_t bn_end = bn+blocksize;
            bn_end > N ? bn_end = N : 0;

            const size_t bn_end_c = bn_end;

            for(size_t nn = bn; nn< bn_end_c; nn++)
            {
                for(size_t mm = bm; mm< bm_end_c; mm++)
                {
                    assert(nn + mm*N < M*N);
                    S[nn + mm*N] = V[mm + nn*M];
                }
            }
        }
    }

    return O;
}

fimo * fim_shiftdim(const fimo * restrict I)
{
    const float * V = I->V;
    const size_t M = I->M;
    const size_t N = I->N;
    const size_t P = I->P;

    /* Output image */
    fimo * O = calloc(1, sizeof(fimo));
    assert(O != NULL);
    O->V = fim_malloc(M*N*P*sizeof(float));

    O->M = N;
    O->N = P;
    O->P = M;
    //printf("%zu, %zu, %zu -> %zu, %zu, %zu\n", I->M, I->N, I->P, O->M, O->N, O->P);
    /* The shifted volume */
    float * S = O->V;

    const size_t blocksize = 2*64;

#pragma omp parallel for shared(V, S) schedule(dynamic)
    for(size_t pp = 0; pp< P; pp++)
    {
        /* Blocks are of size blocksize x blocksize in M and N */
        for(size_t bm = 0; bm<M; bm = bm+blocksize)
        {
            size_t bm_end = bm+blocksize;
            bm_end > M ? bm_end = M : 0;

            const size_t bm_end_c = bm_end;

            for(size_t bn = 0; bn<N; bn = bn+blocksize)
            {
                size_t bn_end = bn+blocksize;
                bn_end > N ? bn_end = N : 0;

                const size_t bn_end_c = bn_end;

                for(size_t nn = bn; nn< bn_end_c; nn++)
                {
                    for(size_t mm = bm; mm< bm_end_c; mm++)
                    {
                        S[nn + pp*N + mm*N*P] = V[mm + nn*M + pp*M*N];
                    }
                }
            }
        }
    }

    return O;
}

/* 1st Eigenvalue of the matrix [a, c ; c, b]*/
static float eig_sym_22_1st(float a, float b, float c)
{
    return 0.5*(a+b + sqrt(pow(a-b,2) + 4*pow(c,2)) );
}

/* 2nd Eigenvalue of the matrix [a, c ; c, b]*/
static float eig_sym_22_2nd(float a, float b, float c)
{
    return 0.5*(a+b - sqrt(pow(a-b,2) + 4*pow(c,2)) );
}

static float total_gm(const float * I0, size_t M, size_t N, float sigma)
{
    fimo * I = fim_image_from_array(I0, M, N, 1);
    fimo * dx = fimo_partial(I, 0, sigma);
    fimo * dy = fimo_partial(I, 1, sigma);
    fimo_free(I);

    double gm = 0;
    for(size_t kk = 0; kk<M*N; kk++)
    {
        gm += sqrt( pow(dx->V[kk], 2) + pow(dy->V[kk], 2));
    }

    fimo_free(dx);
    fimo_free(dy);
    return (float) gm;
}

float * fim_focus_gm(const fimo * I, float sigma)
{
    float * gm = malloc(I->P*sizeof(float));
    assert(gm != NULL);

    /* Since the gpartial functions are parallel per plane we can
     * loop over the planes in parallel here instead
     */
#pragma omp parallel for
    for(size_t kk = 0; kk<I->P; kk++)
    {
        gm[kk] = total_gm(I->V + kk*I->M*I->N, I->M, I->N, sigma);
    }
    return gm;
}

float * fim_auto_zcrop(const float * V,
                       const size_t M, const size_t N, const size_t P,
                       const size_t newP)
{
    fimo * fV = calloc(1, sizeof(fimo));
    assert(fV != NULL);
    fV->V = (float*) V;
    fV->M = M;
    fV->N = N;
    fV->P = P;
    float * focus = fim_focus_gm(fV, 3);
    int slice = float_arg_max(focus, P);
    free(focus);
    int first = slice - (newP-1)/2;
    first < 0 ? first = 0 : 0;
    int last = first + newP - 1;
    if(last <= (int) P)
    {
        last = last - (last - P) -1;
        first = last - newP;
    }
    free(fV);
    float * IZ = fim_malloc(M*N*newP*sizeof(float));
    if(IZ == NULL)
    {
        fprintf(stderr, "Failed to allocate memory in fim_auto_zcrop\n");
        return NULL;
    }

    for(size_t kk = 0; kk < newP; kk++)
    {
        memcpy(IZ + kk*M*N,
               V + (first+kk)*M*N,
               M*N*sizeof(float));
    }
    return IZ;
}

float * fim_zcrop(const float * V,
                  const size_t M, const size_t N, const size_t P,
                  const size_t zcrop)
{
    if(2*zcrop >= P)
    {
        fprintf(stderr, "Impossible zcrop value passed to fim_zcrop\n");
        return NULL;
    }
    size_t newP = P - 2*zcrop;
    float * IZ = fim_malloc(M*N*newP*sizeof(float));
    if(IZ == NULL)
    {
        fprintf(stderr, "Failed to allocate memory in fim_zcrop\n");
        return NULL;
    }

    for(size_t kk = 0; kk < newP; kk++)
    {
        memcpy(IZ + kk*M*N,
               V + (kk+zcrop)*M*N,
               M*N*sizeof(float));
    }
    return IZ;
}

fimo * fimo_get_plane(const fimo * A, int plane)
{
    if(plane < 0)
    {
        printf("fimo_get_plane: Can't extract plane %d\n", plane);
        return NULL;
    }
    if(plane >= (int) A->P)
    {
        printf("fimo_get_plane: Can't extract plane %d from an image with %d planes\n",
               plane, (int) A->P);
        return NULL;
    }
    float * P = fim_malloc(A->M*A->N*sizeof(float));
    if(P == NULL)
    {
        return NULL;
    }
    memcpy(P,
           A->V + plane*A->M*A->N,
           A->M*A->N*sizeof(float));
    fimo * Z = calloc(1, sizeof(fimo));
    if(Z == NULL)
    {
        free(P);
        return NULL;
    }
    Z->V = P;
    Z->M = A->M;
    Z->N = A->N;
    Z->P = 1;
    return Z;
}

ftab_t * fim_features_2d(const fimo * fI,
                         const float * sigmas, int nsigma)
{
    int debug = 0; /* Write out the features  */
    const char debug_image_name[] = "fim_features_2d_debug_image.tif";

    if(fI->P != 1)
    {
        fprintf(stderr, "fim_features_2d can only work with 2D images, supplied image was %zu x %zu x %zu\n",
                fI->M, fI->N, fI->P);
        exit(EXIT_FAILURE);
    }
    /* For varying sigmas */
    //float sigmas[] = {0.3, 0.7, 1, 1.6, 3.5, 5, 10};
    //int nsigma = 7;
    //float sigmas[] = {0.3, 1, 3.5, 10};
    //int nsigma = 4;
    int f_per_s = 7; /* Features per sigma */
    int nfeatures = nsigma*f_per_s;
    //printf("Will produce %d features\n", nfeatures);
    ftab_t * T = malloc(sizeof(ftab_t));
    assert(T != NULL);
    T->nrow = fI->M*fI->N;
    T->ncol = nfeatures;
    T->T = malloc(T->ncol*T->nrow*sizeof(float));
    assert(T->T != NULL);
    T->nrow_alloc = T->nrow;
    T->colnames = NULL;

    const size_t M = fI->M;
    const size_t N = fI->N;
    const size_t P = 1;

    float * debug_image = NULL;
    if(debug)
    {
        debug_image = malloc(M*N*nfeatures*sizeof(float));
        assert(debug_image != NULL);
    }
    printf("Extracting features ... \n");
    int col = 0;
    char * sbuff = malloc(1024);
    assert(sbuff != NULL);
    for(int ss = 0; ss<nsigma; ss++)
    {
        float sigma = sigmas[ss];
        //printf("%f ", sigma); fflush(stdout);

        fimo * G = fimo_copy(fI);

        fim_gsmooth(G->V, M, N, P, sigma);

        fimo * dx = fimo_partial(fI, 0, sigma);
        fimo * dy = fimo_partial(fI, 1, sigma);
        fimo * ddx = fimo_partial(dx, 0, sigma);
        fimo * ddy = fimo_partial(dy, 1, sigma);
        fimo * dxdy = fimo_partial(dx, 1, sigma);

        float * value = G->V;
        /* Gaussian, 1f */
        ftab_set_coldata(T, col, value);
        debug == 1 ? memcpy(debug_image + M*N*col, value, M*N*sizeof(float)) : 0;
        sprintf(sbuff, "s%.1f_Gaussian", sigma);
        ftab_set_colname(T, col++, sbuff);
        debug == 1 ? memcpy(debug_image + M*N*col, value, M*N*sizeof(float)) : 0;


        /* LoG, 1f, ddx+ddy*/
        for(size_t kk = 0; kk<M*N; kk++)
        {
            value[kk] = ddx->V[kk] + ddy->V[kk];
        }
        ftab_set_coldata(T, col, value);
        debug == 1 ? memcpy(debug_image + M*N*col, value, M*N*sizeof(float)) : 0;
        sprintf(sbuff, "s%.1f_LoG", sigma);
        ftab_set_colname(T, col++, sbuff);

        /* Gradient Magnitude, 1f (dx.^2 + dy.^2).^(1/2) */
        for(size_t kk = 0; kk<M*N; kk++)
        {
            value[kk] = sqrt( pow(dx->V[kk], 2) + pow(dy->V[kk], 2));
        }
        ftab_set_coldata(T, col, value);
        debug == 1 ? memcpy(debug_image + M*N*col, value, M*N*sizeof(float)) : 0;
        sprintf(sbuff, "s%.1f_GM", sigma);
        ftab_set_colname(T, col++, sbuff);

        /* Eigenvalues of the Structure Tensor [dx*dx, dx*dy; dx*dy, dy*dy], 2f*/
        for(size_t kk = 0; kk<M*N; kk++)
        {
            float a = dx->V[kk]*dx->V[kk];
            float b = dy->V[kk]*dy->V[kk];
            float c = dx->V[kk]*dy->V[kk];
            value[kk] = eig_sym_22_1st(a, b, c);
        }
        ftab_set_coldata(T, col, value);
        debug == 1 ? memcpy(debug_image + M*N*col, value, M*N*sizeof(float)) : 0;
        sprintf(sbuff, "s%.1f_ST_EV_1", sigma);
        ftab_set_colname(T, col++, sbuff);

        for(size_t kk = 0; kk<M*N; kk++)
        {
            float a = dx->V[kk]*dx->V[kk];
            float b = dy->V[kk]*dy->V[kk];
            float c = dx->V[kk]*dy->V[kk];
            value[kk] = eig_sym_22_2nd(a, b, c);
        }
        ftab_set_coldata(T, col, value);
        debug == 1 ? memcpy(debug_image + M*N*col, value, M*N*sizeof(float)) : 0;
        sprintf(sbuff, "s%.1f_ST_EV_2", sigma);
        ftab_set_colname(T, col++, sbuff);


        /* Eigenvalues of the Hessian of Gaussians [ddx, dxdy; dxdy ddy], 2f */
        for(size_t kk = 0; kk<M*N; kk++)
        {
            float a = ddx->V[kk];
            float b = ddy->V[kk];
            float c = dxdy->V[kk];
            value[kk] = eig_sym_22_1st(a, b, c);
        }
        ftab_set_coldata(T, col, value);
        debug == 1 ? memcpy(debug_image + M*N*col, value, M*N*sizeof(float)) : 0;
        sprintf(sbuff, "s%.1f_HE_EV_1", sigma);
        ftab_set_colname(T, col++, sbuff);

        for(size_t kk = 0; kk<M*N; kk++)
        {
            float a = ddx->V[kk];
            float b = ddy->V[kk];
            float c = dxdy->V[kk];
            value[kk] = eig_sym_22_2nd(a, b, c);
        }
        ftab_set_coldata(T, col, value);
        debug == 1 ? memcpy(debug_image + M*N*col, value, M*N*sizeof(float)) : 0;
        sprintf(sbuff, "s%.1f_HE_EV_2", sigma);
        ftab_set_colname(T, col++, sbuff);
        fimo_free(dx);
        fimo_free(dy);
        fimo_free(ddx);
        fimo_free(ddy);
        fimo_free(dxdy);
        fimo_free(G);
    }
    free(sbuff); /* Free string buffer */
    printf("\n");
    if(debug)
    {
        fim_tiff_write_float(debug_image_name,
                             debug_image, NULL,
                             M, N, nfeatures);
        fim_free(debug_image);
    }
    return T;
}

fimo * fimo_tiff_read(const char * file)
{
    i64 M, N, P;

    float * V = fim_tiff_read(file, NULL, &M, &N, &P, 0);
    if(V == NULL)
    {
        return NULL;
    }
    fimo * I = calloc(1, sizeof(fimo));
    assert(I != NULL);
    I->V = V;
    I->M = M;
    I->N = N;
    I->P = P;
    return I;
}


void fim_features_2d_ut()
{
    printf("-> fim_features_2d_ut\n");
    char file[] = "features_test.tif";
    FILE * f = fopen(file, "r");
    if(f == NULL)
    {
        printf("can't open %s aborting fim_features_2d_ut\n", file);
        return;
    }
    fclose(f);

    int64_t M = 0;
    int64_t N = 0;
    int64_t P = 0;

    float * V = fim_tiff_read(file, NULL, &M, &N, &P, 0);
    printf("%s is %" PRId64 " %" PRId64 " %" PRId64 " image\n", file, M, N, P);
    fimo * I = malloc(sizeof(fimo));
    assert(I != NULL);
    I->V = V;
    I->M = M;
    I->N = N;
    I->P = P;

    float * sigma = calloc(2, sizeof(float));
    assert(sigma != NULL);
    sigma[0] = 3.5;
    sigma[1] = sqrt(2)*sigma[0];
    ftab_t * T = fim_features_2d(I, sigma, 2);
    free(sigma);
    T->nrow = 10; // just to print the first 10 ...
    ftab_print(stdout, T, ",");
    ftab_free(T);
    fimo_free(I);
}


void fim_argmax_max_ut()
{
    struct timespec tstart, tend;
    size_t M = 2048;
    size_t N = 2048;
    size_t P = 60;
    size_t MNP = M*N*P;
    float * A = fim_malloc(MNP*sizeof(float));
    assert(A != NULL);

    for(size_t kk = 0 ; kk < MNP; kk ++)
    {
        A[kk] = (float) rand()/ (float) RAND_MAX;
    }
    dw_gettime(&tstart);

    float ref_max = A[0];
    size_t ref_argmax = 0;
    for(size_t kk = 0 ; kk < MNP; kk++)
    {
        if(A[kk] > ref_max)
        {
            ref_max = A[kk];
            ref_argmax = kk;
        }
    }
    int64_t ref_amP = ref_argmax / (M*N);
    int64_t ref_amN = (ref_argmax - (ref_amP*M*N)) / M;
    int64_t ref_amM = ref_argmax - ref_amP*M*N - ref_amN*M;
    dw_gettime(&tend);

    double t_ref_amax = clockdiff(&tend, &tstart);


    dw_gettime(&tstart);
    int64_t amM, amN, amP;
    float lib_max;
    fim_argmax_max(A, M, N, P, &amM, &amN, &amP, &lib_max);

    dw_gettime(&tend);
    double t_lib_amax = clockdiff(&tend, &tstart);

    if( (lib_max == ref_max) &&
        (amM + amN*M + amP*M*N == ref_argmax))
    {
        printf("lib: A(%" PRId64 ", %" PRId64 ", %" PRId64 ") = %f\n", amM, amN, amP, lib_max);
        printf("ref: A(%zu) = A(%" PRId64 ", %" PRId64 ", %" PRId64 ") = %f\n", ref_argmax,
               ref_amM, ref_amN, ref_amP, ref_max);
        printf("fim_argmax_max_ut passed (%e s, ref-implementation: %e s)\n",
               t_lib_amax, t_ref_amax);
    } else {
        fprintf(stderr, "fim_argmax_max_ut failed\n");
        printf("lib: A(%" PRId64 ", %" PRId64 ", %" PRId64 ") = %f\n", amM, amN, amP, lib_max);
        printf("ref: A(%zu) = A(%" PRId64 ", %" PRId64 ", %" PRId64 ") = %f\n", ref_argmax,
               ref_amM, ref_amN, ref_amP, ref_max);
        exit(EXIT_FAILURE);
    }
    fim_free(A);
    return;
}


void fim_max_ut()
{
    struct timespec tstart, tend;
    size_t N = 2024*1024*60;
    float * A = fim_malloc(N*sizeof(float));
    assert(A != NULL);

    for(size_t kk = 0 ; kk < N; kk ++)
    {
        A[kk] = (float) rand()/ (float) RAND_MAX;
    }

    dw_gettime(&tstart);
    float ref_max = A[0];
    for(size_t kk = 0 ; kk < N; kk++)
    {
        if(A[kk] > ref_max)
        {
            ref_max = A[kk];
        }
    }

    dw_gettime(&tend);
    double t_ref_max = clockdiff(&tend, &tstart);

    dw_gettime(&tstart);
    float lib_max = fim_max(A, N);
    dw_gettime(&tend);
    double t_lib_max = clockdiff(&tend, &tstart);

    if(lib_max == ref_max)
    {
        printf("fim_max_ut passed (%e s, ref-implementation: %e s)\n", t_lib_max, t_ref_max);
    } else {
        fprintf(stderr, "fim_max_ut failed\n");
        exit(EXIT_FAILURE);
    }
    fim_free(A);
    return;
}

void fim_min_ut()
{
    struct timespec tstart, tend;
    size_t N = 2024*1024*60;
    float * A = fim_malloc(N*sizeof(float));
    assert(A != NULL);

    for(size_t kk = 0 ; kk < N; kk ++)
    {
        A[kk] = (float) rand()/ (float) RAND_MAX;
    }
    dw_gettime(&tstart);
    float ref_min = A[0];
    for(size_t kk = 0 ; kk < N; kk++)
    {
        if(A[kk] < ref_min)
        {
            ref_min = A[kk];
        }
    }
    dw_gettime(&tend);
    double t_ref_min = clockdiff(&tend, &tstart);

    dw_gettime(&tstart);
    float lib_min = fim_min(A, N);
    dw_gettime(&tend);
    double t_lib_min = clockdiff(&tend, &tstart);

    if(lib_min == ref_min)
    {
        printf("fim_min_ut passed (%e s, ref-implementation: %e s)\n",
               t_lib_min, t_ref_min);
    } else {
        fprintf(stderr, "fim_min_ut failed\n");
        exit(EXIT_FAILURE);
    }
    fim_free(A);
    return;
}

void fim_conv1_ut(fim_boundary_condition bc)
{
    printf("fim_conv1_ut(), bc: %s\n", fim_boundary_condition_str(bc));
    size_t nV = 7;
    size_t nK = 3;
    size_t stride = 1;
    float * V = calloc(nV, sizeof(float));
    assert(V != NULL);
    for(size_t kk = 0; kk < nV; kk++)
    { V[kk] = kk+1; }
    printf("V=");fim_show(V, 1, nV, 1);
    float * K = calloc(nK, sizeof(float));
    assert(K != NULL);
    K[0] = 1;
    K[1] = 1;
    K[2] = 1;
    printf("K=");fim_show(K, 1, nK, 1);
    fim_conv1(V, nV, stride,
              K, nK,
              NULL, bc);
    printf("V*K=");fim_show(V, 1, nV, 1);
    if(bc == FIM_BC_ZEROS)
    {
        assert(V[0] == 1+2);
        assert(V[6] == 6+7);
    }
    if(bc == FIM_BC_SYMMETRIC_MIRROR)
    {
        assert(V[0] == 2+1+2);
        assert(V[6] == 6+7+6);
    }
    if(bc == FIM_BC_VALID)
    {
        assert(V[0] == 0);
        assert(V[6] == 0);
    }
    if(bc == FIM_BC_PERIODIC)
    {
        assert(V[0] == 1 + 2 + 7);
        assert(V[6] == 6 + 7 + 1);
    }
    if(bc == FIM_BC_WEIGHTED)
    {
        assert(fabs(V[0] - (1.0+2.0)*3.0/2.0) < 1e-5);
        assert(fabs(V[6] - (6.0+7.0)*3.0/2.0) < 1e-5);
    }
    free(V);
    free(K);
    return;
}


fimo * fimo_transpose(const fimo * restrict A)
{
    fimo * B = fimo_copy(A);
    size_t M = A->M;
    size_t N = A->N;
    size_t P = A->P;

    B->N = M;
    B->M = N;
    B->P = P;

    for(size_t pp = 0; pp < P; pp++)
    {
        for(size_t mm = 0; mm < M; mm++)
        {
            for(size_t nn = 0; nn < N; nn++)
            {
                B->V[nn + mm*N + pp*M*N] =
                    A->V[mm + nn*M + pp*M*N];
            }
        }
    }
    return B;
}

int fimo_tiff_write(const fimo * I, const char * fName)
{
    return fim_tiff_write_float(fName, I->V, NULL, I->M, I->N, I->P);
}

void fimo_blit_2D(fimo * A, const fimo * B, size_t x0, size_t y0)
{
    assert(A != NULL);
    assert(B != NULL);
    assert(A->P == 1);
    assert(B->P == 1);

    assert(A->M >= B->M);
    assert(A->N >= B->N);

    assert(A->M >= x0 + B->M);
    assert(A->N >= y0 + B->N);

    for(size_t xx = 0; xx < B->M; xx++)
    {
        for( size_t yy = 0; yy < B->N ; yy++)
        {
            A->V[xx+x0 + A->M*(yy+y0)] = B->V[xx + B->M*yy];
        }
    }

    return;
}


float
fim_interp3_trilinear(const float * restrict A,
                      const size_t M, const size_t N, const size_t P,
                      const float x, const float y, const float z)
{

    /* We would like to calculate
     * A[x + M*y + N*z] by linear interpolation.
     * For that we need 8 points from A
     */


    if(x < 0 || x > (M - 1)
       || y < 0 || y > (N - 1)
       || z < 0 || z > (P - 1))
    {
        return 0;
    }

    int x0 = x;
    float wx = 1.0 - ( x-(float) x0 );
    int y0 = y;
    float wy = 1.0 - ( y-(float) y0 );
    int z0 = z;
    float wz = 1.0 - ( z- (float) z0 );

    return wz*(
        wy*(
            wx*A[x0 + y0*M + z0*M*N]
            + (1.0-wx)*A[x0+1 + y0*M + z0*M*N])
        +(1.0-wy)*(wx*A[x0 + (y0+1)*M + z0*M*N] +
                   (1.0-wx)*A[x0+1 + (y0+1)*M + z0*M*N]))
        + (1.0-wz)*(
            wy*(
                wx*A[x0 + y0*M + (z0+1)*M*N]
                + (1.0-wx)*A[x0+1 + y0*M + (z0+1)*M*N])
            + (1.0-wy)*(wx*A[x0 + (y0+1)*M + (z0+1)*M*N]
                        + (1.0-wx)*A[x0+1 + (y0+1)*M + (z0+1)*M*N]));
}


/* Low precision covariance calculation */
float
fim_covariance_lp(const float * X, const float * Y, size_t n)
{
    assert(X != NULL);
    assert(Y != NULL);
    double mx = fim_mean(X, n);
    double my = mx;
    if(X != Y)
    {
        my = fim_mean(Y, n);
    }
    double covar = 0;
    for(size_t kk = 0; kk < n; kk++)
    {
        covar += (X[kk]-mx)*(Y[kk]-my);
    }
    covar /= ((double) n - 1.0);
    return covar;
}

float
fim_dot_lateral_circularity(const float * I,
                            size_t M, size_t N, size_t P,
                            double x, double y, double z,
                            double sigma)
{
    assert(sigma > 0);
    assert(I != NULL);
    // Check z-coordinate
    const int zi = round(z);
    if(zi < 0)
    {
        return -1;
    }
    if((size_t) zi >= P)
    {
        return -1;
    }

    /* Radius of the windowing function */
    sigma < 1.0 ? sigma = 1.0 : 0;

    const double radius = 1.75*sigma;

    const int s = (int) ( (2.0*sigma) + 1.0 );
    const int n = 2*s + 1;
    assert(n > 0);

    float * X = calloc(n*n, sizeof(float));
    assert(X != NULL);

    float * Y = calloc(n*n, sizeof(float));
    assert(Y != NULL);

    const int ix = round(x);
    const int iy = round(y);

    int xlow = ix - s;
    ix - s < 0 ? xlow = 0 : 0;
    int xhigh = ix + s;
    ix + s >= (int) M ? xhigh = M-1 : 0;

    int ylow = iy - s;
    iy - s < 0 ? ylow = 0 : 0;
    int yhigh = iy + s;
    iy + s >= (int) N ? yhigh = N-1 : 0;

    //printf("[%d %d] x [%d %d]\n", xlow, xhigh, ylow, yhigh);

    /* Get the max and min values of the pixels */
    double minI = 1e99;
    double maxI = -1e99;
    for(int yy = ylow; yy <= yhigh; yy++)
    {
        for(int xx = xlow; xx <= xhigh; xx++)
        {
            double p = I[xx + yy*M + zi*M*N];
            p > maxI ? maxI = p : 0;
            p < minI ? minI = p : 0;
        }
    }
    //printf("data in range [%f, %f]\n", minI, maxI);

    if(maxI <= minI)
    {
        /* If the pixel data is constant we can't estimate anything */
        free(X);
        free(Y);
        return -1;
    }
    assert(maxI > minI);
    assert( (yhigh-ylow +1)*(xhigh - xlow + 1) <= n*n);

    size_t idx = 0;
    for(int yy = ylow; yy <= yhigh; yy++)
    {
        for(int xx = xlow; xx <= xhigh; xx++)
        {

            double r = sqrtf(
                powf( (double) yy - y, 2) + powf( (double) xx - x, 2)
                );
            double rw = (radius - r); /* A soft threshold */
            rw  > 1.0 ? rw = 1 : 0;
            rw < 0.0 ? rw = 0 : 0;

            // normalized pixel value to [0, 1]
            double w = rw*(I[xx + yy*M + zi*M*N] - minI)/(maxI-minI);
            X[idx] = w* ((double) xx - x);
            Y[idx] = w* ((double) yy - y);
            if(0)
            {
                printf("I[%d, %d]=%f, rw= %f  x=(%f, %f)  %f, %f\n",
                       xx, yy,
                       I[xx + yy*M + zi*M*N], rw,
                       (double) xx - x, (double) yy - y,
                       X[idx], Y[idx]);
            }
            idx++;
        }
    }

    double a = fim_covariance_lp(X, X, n*n);
    double b = fim_covariance_lp(Y, Y, n*n);
    double c = fim_covariance_lp(X, Y, n*n);

    // printf("cov = [%f %f ; %f %f]\n", a, b, b, c);
    free(X);
    free(Y);

    double l1 = (a+b)/2 + sqrt(pow(c, 2)  + pow(a+b,2)/4 - a*b);
    double l2 = (a+b)/2 - sqrt(pow(c, 2)  + pow(a+b,2)/4 - a*b);
    //printf("l1=%f, l2 = %f\n", l1, l2);
    return l2 / l1;
}

static void fim_covariance_lp_ut()
{
    size_t n = 11;
    float * X = calloc(n, sizeof(float));
    assert(X != NULL);
    float * Y = calloc(n, sizeof(float));
    assert(Y != NULL);
    for(size_t kk = 0; kk < n; kk++)
    {
        X[kk] = (double) rand() / (double) RAND_MAX;
        Y[kk] = (double) rand() / (double) RAND_MAX;
    }


    for(size_t kk = 0; kk < n; kk++)
    {
        if(kk == 0)
        {
            printf("X=[");
        }
        printf("%f, %f", X[kk], Y[kk]);
        if( kk+1 == n)
        {
            printf("];");
        }
        printf("\n");
    }

    printf("cov = [%f, %f; %f, %f]\n",
           fim_covariance_lp(X, X, n),
           fim_covariance_lp(X, Y, n),
           fim_covariance_lp(Y, X, n),
           fim_covariance_lp(Y, Y, n));

    free(X);
    free(Y);
}

void fim_dot_lateral_circularity_ut()
{
    int s = 7;
    int n = 2*s + 1;
    float * G = calloc(n*n, sizeof(float));
    assert(G != NULL);
    for(int xx = 0; xx < n; xx++)
    {
        for(int yy = 0; yy < n; yy++)
        {
            float r2 = n-sqrt((xx-s)*(xx-s) + (yy-s)*(yy-s));
            G[xx + n*yy] = r2;
        }
    }
    double circ = fim_dot_lateral_circularity(G, n, n, 1,
                                              s, s, 0,
                                              1.5);
    printf("circ = %f\n", circ);

    free(G);

    FILE * fid = fopen("G.f32", "rb");
    if(fid == NULL)
    {
        return;
    }
    printf("Testing on data in G.f32 which is assumed to be a square shaped 2D image\n");
    fseek(fid, 0, SEEK_END);
    size_t nb = ftell(fid);
    rewind(fid);
    float * G2 = calloc(nb/sizeof(float), sizeof(float));
    size_t nread = fread(G2, nb/sizeof(float), sizeof(float), fid);

    if(nread != nb/sizeof(float))
    {
        printf("Error reading G.f32\n");
        exit(EXIT_FAILURE);
    }
    size_t side = sqrt(nb/sizeof(float));
    circ = fim_dot_lateral_circularity(G2, side, side, 1,
                                       ((double) side - 1.0) /2,
                                       ((double) side - 1.0) /2,
                                       0,
                                       2);
    free(G2);
    printf("circ = %f\n", circ);
    return;
}

static void fim_DoH_ut(void)
{
    size_t M = 128;
    size_t N = 64;
    size_t P = 1;
    float * X = calloc(M*N*P, sizeof(float));
    assert(X != NULL);
    float * DoH = fim_DoH(X, M, N, P, 1, 1);
    printf("DoH[0] = %f\n", DoH[0]);
    free(X);
    free(DoH);
    return;
}

void fim_ut()
{
    fim_DoH_ut();
    return;
    fim_covariance_lp_ut();
    fim_dot_lateral_circularity_ut();

    fim_conv1_ut(FIM_BC_SYMMETRIC_MIRROR);
    fim_conv1_ut(FIM_BC_ZEROS);
    fim_conv1_ut(FIM_BC_VALID);
    fim_conv1_ut(FIM_BC_PERIODIC);
    fim_conv1_ut(FIM_BC_WEIGHTED);
    fim_argmax_max_ut();
    fim_min_ut();
    fim_max_ut();
    fim_LoG_ut();
    fim_features_2d_ut();
    fim_conv1_vector_ut();
    fim_conncomp6_ut();
    exit(EXIT_SUCCESS);
    fim_otsu_ut();
    exit(EXIT_SUCCESS);
    fim_flipall_ut();
    shift_vector_ut();
    size_t N = 0;
    float sigma = 1;
    float * K = gaussian_kernel(sigma, &N);
    assert(N>0);
    printf("gaussian_kernel, sigma=%f\n", sigma);
    show_vec(K, N);
    fim_free(K);
    sigma = 0.5;
    K = gaussian_kernel(sigma, &N);
    printf("gaussian_kernel, sigma=%f\n", sigma);
    show_vec(K, N);

    int nV = 10;
    float * V = fim_malloc(nV*sizeof(float));
    assert(V != NULL);
    for(int kk = 0; kk<nV; kk++)
    {
        V[kk] = kk;
    }
    printf("V=");
    show_vec(V, nV);
    fim_conv1_vector(V, 1, NULL, nV, K, N, 1);
    printf("V*K = ");
    show_vec(V, nV);

    fim_free(V);
    fim_free(K);

    fim_cumsum_ut();
    fim_local_sum_ut();
    //exit(EXIT_FAILURE);
    myfftw_start(1, 1, stdout);
    fim_xcorr2_ut();
    myfftw_stop();
}

float *
fim_read_npy(const char * filename,
             int64_t * M, int64_t * N, int64_t * P,
             int verbose)
{
    npio_t * npy = npio_load(filename);
    if(npy == NULL)
    {
        fprintf(stderr, "Error reading %s as a npy file\n", filename);
        return NULL;
    }
    if(npy->ndim != 3)
    {
        fprintf(stderr, "Error reading %s, not 3D\n", filename);
        goto fail;
    }
    *M = npy->shape[2];
    *N = npy->shape[1];
    *P = npy->shape[0];

    size_t nel = npy->nel;
    if( (int64_t) nel != M[0]*N[0]*P[0])
    {
        fprintf(stderr, "Internal error in %s %d\n", __FILE__, __LINE__);
        goto fail;
    }

    float * V = fim_malloc(nel*sizeof(float));

    if(npy->dtype == NPIO_F32)
    {
        memcpy(V, npy->data, nel*sizeof(float));
        goto success;
    }

    if(npy->dtype == NPIO_U8)
    {
        uint8_t * IN = (uint8_t * ) npy->data;
        for(size_t kk = 0; kk < nel; kk++)
        {
            V[kk] = (float) IN[kk];
        }
        goto success;
    }

    if(npy->dtype == NPIO_U16)
    {
        uint16_t * IN = (uint16_t * ) npy->data;
        for(size_t kk = 0; kk < nel; kk++)
        {
            V[kk] = (float) IN[kk];
        }
        goto success;
    }

    fprintf(stderr, "Unsupported data type in %s\n", filename);
    goto fail;

fail:
    if(verbose > 0)
    {
        npio_print(stderr, npy);
    }
    npio_free(npy);
    return NULL;
success:
    npio_free(npy);
    return V;
}

float *
fim_imread(const char * filename,
           ttags * T,
           int64_t * M, int64_t * N, int64_t * P,
           int verbose)
{
    const size_t n = strlen(filename);
    assert(n > 0);

    if(filename[n-1] == 'y' || filename[n-1] == 'Y')
    {
        return fim_read_npy(filename, M, N, P, verbose);
    }
    return fim_tiff_read(filename, T, M, N, P, verbose);
}

int
fim_imwrite_f32(const char * outname,
                const float * V,
                const ttags * T,
                int64_t M, int64_t N, int64_t P)
{
    const size_t n = strlen(outname);
    assert(n > 0);

    if(outname[n-1] == 'y' || outname[n-1] == 'Y')
    {
        int shape[3] = {P, N, M};
        return npio_write(outname, 3, shape, (void *) V,
                          NPIO_F32, NPIO_F32);
    }
    return fim_tiff_write_float(outname, V, T, M, N, P);
}

/* Write as uint16 data. If scaling <= 0 automatic scaling is used,
 * i.e. using the full range. Else the provided value is used */
int
fim_imwrite_u16(const char * outname,
                const float * V,
                const ttags * T,
                int64_t M, int64_t N, int64_t P,
                float scaling)
{
    const size_t n = strlen(outname);
    assert(n > 0);
    assert(M > 0);
    assert(N > 0);
    assert(P > 0);

    if(outname[n-1] == 'y' || outname[n-1] == 'Y')
    {
        uint16_t * I = calloc(M*N*P, sizeof(uint16_t));

        for(int64_t kk = 0; kk < M*N*P; kk++)
        {
            I[kk] = V[kk]*scaling;
        }
        int shape[3] = {P, N, M};
        return npio_write(outname, 3, shape, (void *) I, NPIO_U16, NPIO_U16);
        free(I);
    }
    return fim_tiff_write_opt(outname, V, T, M, N, P, scaling);
}

int
fim_imread_size(const char * filename,
                int64_t * M, int64_t * N, int64_t * P)
{
    const size_t n = strlen(filename);
    assert(n > 0);

    if(filename[n-1] == 'y' || filename[n-1] == 'Y')
    {
        npio_t * npy = npio_load_metadata(filename);
        if(npy->ndim != 3)
        {
            npio_free(npy);
            return -1;
        }
        *P = npy->shape[0];
        *N = npy->shape[1];
        *M = npy->shape[2];
        return 0;
    }
    return fim_tiff_get_size(filename, M, N, P);
}
