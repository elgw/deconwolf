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

static int fim_verbose = 0;

static float * gaussian_kernel(float sigma, size_t * nK);
static void cumsum_array(float * A, size_t N, size_t stride);
static void fim_show(float * A, size_t M, size_t N, size_t P);
static void fim_show_int(int * A, size_t M, size_t N, size_t P);

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
    return _aligned_malloc(nbytes, FIM_ALIGNMENT);
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


void fim_delete(fim_t * F)
{
    if(F != NULL)
    {
        fim_free(F->V);
        F->V = NULL;
        free(F);
    }
    return;
}

fim_t * fim_image_from_array(const float * restrict V, size_t M, size_t N, size_t P)
{
    fim_t * I = malloc(sizeof(fim_t));
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


void fim_stats(const float * A, const size_t N)
{
    printf("min: %f mean: %f, max: %f\n",
           fim_min(A, N),
           fim_mean(A, N),
           fim_max(A, N));
    return;
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

fim_t * fimt_copy(const fim_t * restrict F)
{
    fim_t * C = malloc(sizeof(fim_t));
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
    float * A = fim_malloc(N*sizeof(float));
    //memset(A, 0, N*sizeof(float));
    return A;
}

fim_t * fimt_zeros(const size_t M, const size_t N, const size_t P)
// Allocate and return an array of N floats
{
    size_t n = M*N*P;
    fim_t * F = malloc(sizeof(fim_t));
    assert(F != NULL);
    F->V = fim_malloc(n*sizeof(float));

    F->M = M;
    F->N = N;
    F->P = P;
    //memset(F->V, 0, n*sizeof(float));
    return F;
}

size_t fimt_nel(fim_t * F)
{
    return F->M*F->N*F->P;
}

float fimt_sum(fim_t * F)
{
    return fim_sum(F->V, fimt_nel(F));
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
    printf("Labelled array\n");
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
    for(int nK = 3; nK<8; nK+=2)
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
    printf("-> fim_LoG_ut");
    struct timespec tstart, tend;
    size_t M = 12;
    size_t N = 12;
    size_t P = 12;
    float sigma_l = 3.1;
    float sigma_a = 1.1;
    float * V = fim_malloc(M*N*P*sizeof(float));

    for(size_t kk = 0; kk<M*N*P; kk++)
    {
        V[kk] = kk % 100;
    }


    fim_t * T = fim_image_from_array(V, M, N, P);
    fim_t * S1 = fim_shiftdim(T);
    fim_t * S2 = fim_shiftdim(S1);
    fim_t * S3 = fim_shiftdim(S2);
    for(size_t kk = 0; kk<M*N*P; kk++)
    {
        if(T->V[kk] != S3->V[kk])
        {
            printf("fim_shiftdim does not work at index %zu %f != %f\n",
                   kk, T->V[kk], S3->V[kk]);
            exit(EXIT_FAILURE);
        }
    }
    fim_delete(T);
    T = NULL;
    assert(S1->V != S2->V);
    fim_delete(S1);
    S1 = NULL;
    fim_delete(S2);
    S2 = NULL;
    fim_delete(S3);


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

void fim_conv1_vector(float * restrict V, const int stride, float * restrict W,
                      const size_t nV,
                      const float * restrict K, const size_t nKu, const int normalized)
{
    if(V == NULL || K == NULL)
    {
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
        for(size_t vv = 0;vv<k2; vv++)
        {
            double acc0 = 0;
            double kacc = 0;
            for(size_t kk = k2-vv; kk<nKu; kk++)
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

        /* Central part where K fits completely */
        for(size_t vv = k2 ; vv+k2 < nV; vv++)
        {
            double acc = 0;
            for(size_t kk = 0; kk < nKu; kk++)
            {
                size_t vpos = ((vv-k2)+kk)*stride;
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

    nK[0] = N;
    return K;
}

static float * gaussian_kernel_d1(float sigma, size_t * nK)
{
    /* First derivative of a Gaussian kernel */

    /* Determine the size so that most of the signal is captured */
    int len = 1; /* The total number of elements will be at least 3 */
    while(erf((len+1.0)/sigma) < 1.0-1e-8)
    {
        len++;
    }
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


static float * gaussian_kernel_d2(float sigma, size_t * nK)
{
    /* d2/dx2 Gaussian kernel */
    float * K = gaussian_kernel(sigma, nK);

    int n = (int) nK[0];
    int m = (n-1)/2;
    float s2 = pow(sigma, 2);
    float b = 1.0/(2*s2);

    for(int kk = 0; kk < n; kk++)
    {
        float x = kk-m;
        float x2 = pow(x, 2);
        K[kk] *= -2*b*(2*b*x2-1.0);
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

void fim_gsmooth(float * restrict V, size_t M, size_t N, size_t P, float sigma)
{
    fim_gsmooth_aniso(V, M, N, P, sigma, sigma);
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
        printf("\n");
    }
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
        printf("\n");
    }
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

fim_t * fimt_maxproj(const fim_t * F)
{
    float * _M = fim_maxproj(F->V, F->M, F->N, F->P);
    fim_t * M = malloc(sizeof(fim_t));
    assert(M != NULL);
    M->V = _M;
    M->M = F->M;
    M->N = F->N;
    M->P = 1;
    return M;
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

ftab_t * fim_lmax(const float * I, size_t M, size_t N, size_t P)
{
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

    for(size_t mm = 1; mm+1 < M; mm++)
    {
        for(size_t nn = 1; nn+1 < N; nn++)
        {
            for(size_t pp = 1; pp+1 < P; pp++)
            {
                size_t pos = mm + nn*M + pp*M*N;
                if(I[pos] > strel333_max(I + pos, M, N, P, strel))
                {
                    /* Pos is s a local maxima */
                    float row[4] = {mm, nn, pp, I[pos]};
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
    float min = fim_min(Im, N);
    float max = fim_max(Im, N);
    size_t nbin = pow(2, 16)+1;
    float delta = (max-min) / ((float) nbin);
    float left = min - 0.5*delta-1e-6;
    float right = max + 0.5*delta+1e-6;
    delta = (right-left)/((float) nbin);
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
    double level = otsu(H->C, H->nbin);
    double slevel = H->left + level/(pow(2,16)-1)*(H->right-H->left);
    return (float) slevel;
}

float * fim_otsu(float * Im, size_t M, size_t N)
{
    //float mean = fim_sum(Im, M*N)/( (float) M*N);
    //printf("Mean=%f\n", mean);
    fim_histogram_t * H = fim_histogram(Im, M*N);
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
    fim_free(E2);

    fim_free(E);
    return lab;
}

int fim_convn1(float * restrict V, size_t M, size_t N, size_t P,
               float * K, size_t nK,
               int dim, const int normalized)
{
    if(dim < 0 || dim > 2)
    {
        return EXIT_FAILURE;
    }

    /* Temporary storage/buffer for conv1_vector */
    // TODO: one buffer per thread
    size_t nW = max_size_t(M, max_size_t(N, P));

    int nThreads = 1;
#pragma omp parallel
    {
        nThreads = omp_get_num_threads();
    }


    float * W = fim_malloc(nThreads*nW*sizeof(float));
    assert(W != NULL);

    if(dim == 0)
    {
#pragma omp parallel for
        for(size_t pp = 0; pp < P; pp++)
        {
            float * buff = W+omp_get_thread_num()*nW;
            for(size_t nn = 0; nn < N; nn++)
            {
                fim_conv1_vector(V+pp*(M*N)+nn*M, 1, buff, M, K, nK, normalized);
            }
        }
    }

    if(dim == 1)
    {
#pragma omp parallel for
        for(size_t pp = 0; pp<P; pp++)
        {
            float * buff = W+omp_get_thread_num()*nW;
            for(size_t mm = 0; mm<M; mm++)
            {
                fim_conv1_vector(V + pp*(M*N) + mm, M, buff, N, K, nK, normalized);
            }
        }
    }

    if(dim == 2)
    {
#pragma omp parallel for
        for(size_t mm = 0; mm<M; mm++)
        {
            float * buff = W+omp_get_thread_num()*nW;
            for(size_t nn = 0; nn<N; nn++)
            {
                fim_conv1_vector(V+mm+M*nn, M*N, buff, P, K, nK, normalized);
            }
        }
    }

    fim_free(W);
    return EXIT_SUCCESS;
}


float * conv1_3(const float * V, size_t M, size_t N, size_t P,
                float * K1, size_t nK1,
                float * K2, size_t nK2,
                float * K3, size_t nK3)
{
    const int dim = 0;
    const int norm = 0;

    fim_t * F = fim_image_from_array(V, M, N, P);

    fim_convn1(F->V, F->M, F->N, F->P, K1, nK1, dim, norm);
    fim_t * F2 = fim_shiftdim(F);
    fim_free(F->V);
    fim_free(F);
    fim_convn1(F2->V, F2->M, F2->N, F2->P, K2, nK2, dim, norm);
    fim_t * F3 = fim_shiftdim(F2);
    fim_free(F2->V);
    fim_free(F2);
    fim_convn1(F3->V, F3->M, F3->N, F3->P, K3, nK3, dim, norm);
    fim_t * F4 = fim_shiftdim(F3);
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
    /* Axial filters */
    size_t naG = 0;
    float * aG = gaussian_kernel(sigmaz,  &naG);
    size_t na2;
    float * a2 = gaussian_kernel_d2(sigmaz,  &na2);

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


/* Parital derivative in dimension dim */
fim_t * fimt_partial(const fim_t * F, const int dim, const float sigma)
{

    if(F->P != 1)
    {
        fprintf(stderr, "fimt_partial only supports 2D images (please fix me)\n");
        exit(EXIT_FAILURE);
    }

    size_t nG = 0;
    float * G = gaussian_kernel(sigma, &nG);
    size_t nD1;
    float * D1 = gaussian_kernel_d1(sigma, &nD1);


    fim_t * D = fimt_copy(F);
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

    fprintf(stderr, "Something went wrong in fimt_partial dim=%d\n", dim);
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
    /* Axial filters */
    size_t naG = 0;
    float * aG = gaussian_kernel(sigmaz,  &naG);
    size_t na2;
    float * a2 = gaussian_kernel_d2(sigmaz,  &na2);

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

double * fim_get_line_double(fim_t * I,
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

fim_t * fim_shiftdim(const fim_t * restrict I)
{
    const float * V = I->V;
    const size_t M = I->M;
    const size_t N = I->N;
    const size_t P = I->P;

    /* Output image */
    fim_t * O = malloc(sizeof(fim_t));
    assert(O != NULL);
    O->V = fim_malloc(M*N*P*sizeof(float));

    O->M = N;
    O->N = P;
    O->P = M;
    //printf("%zu, %zu, %zu -> %zu, %zu, %zu\n", I->M, I->N, I->P, O->M, O->N, O->P);
    /* The shifted volume */
    float * S = O->V;

    const size_t blocksize = 2*64;

#pragma omp parallel for shared(V, S)
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
    fim_t * I = fim_image_from_array(I0, M, N, 1);
    fim_t * dx = fimt_partial(I, 0, sigma);
    fim_t * dy = fimt_partial(I, 1, sigma);
    fim_delete(I);

    double gm = 0;
    for(size_t kk = 0; kk<M*N; kk++)
    {
        gm += sqrt( pow(dx->V[kk], 2) + pow(dy->V[kk], 2));
    }

    fim_delete(dx);
    fim_delete(dy);
    return (float) gm;
}

float * fim_focus_gm(const fim_t * I, float sigma)
{
    float * gm = malloc(I->P*sizeof(float));
    assert(gm != NULL);
    for(size_t kk = 0; kk<I->P; kk++)
    {
        gm[kk] = total_gm(I->V + kk*I->M*I->N, I->M, I->N, sigma);
    }
    return gm;
}

ftab_t * fim_features_2d(const fim_t * fI)
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
    float sigmas[] = {3.5};
    int nsigma = 1;
    int f_per_s = 7; /* Features per sigma */
    int nfeatures = nsigma*f_per_s;
    printf("Will produce %d features\n", nfeatures);
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
        printf("%f ", sigma); fflush(stdout);

        fim_t * G = fimt_copy(fI);

        fim_gsmooth(G->V, M, N, P, sigma);

        fim_t * dx = fimt_partial(fI, 0, sigma);
        fim_t * dy = fimt_partial(fI, 1, sigma);
        fim_t * ddx = fimt_partial(dx, 0, sigma);
        fim_t * ddy = fimt_partial(dy, 1, sigma);
        fim_t * dxdy = fimt_partial(dx, 1, sigma);

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
        fim_delete(dx);
        fim_delete(dy);
        fim_delete(ddx);
        fim_delete(ddy);
        fim_delete(dxdy);
        fim_delete(G);
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
    fim_t * I = malloc(sizeof(fim_t));
    assert(I != NULL);
    I->V = V;
    I->M = M;
    I->N = N;
    I->P = P;
    ftab_t * T = fim_features_2d(I);
    T->nrow = 10;
    ftab_print(stdout, T);
    ftab_free(T);
    fim_delete(I);
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

void fim_ut()
{

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

fim_t * fimt_transpose(const fim_t * restrict A)
{
    fim_t * B = fimt_copy(A);
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

int fimt_tiff_write(const fim_t * I, const char * fName)
{
    return fim_tiff_write_float(fName, I->V, NULL, I->M, I->N, I->P);
}

void fimt_blit_2D(fim_t * A, const fim_t * B, size_t x0, size_t y0)
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
