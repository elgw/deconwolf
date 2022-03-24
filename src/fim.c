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


static float * gaussian_kernel(float sigma, size_t * nK);
static void conv1(float * restrict V, int stride, float * restrict W,
                  const size_t nV,
                  const float * restrict K, const size_t nKu);

static void cumsum_array(float * A, size_t N, size_t stride);
static void fim_show(float * A, size_t M, size_t N, size_t P);

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
    // TODO reduction(max:maxV)
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

void fim_argmax(const float * I,
                size_t M, size_t N, size_t P,
                int64_t * _aM, int64_t *_aN, int64_t *_aP)
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
}

float fim_sum(const afloat * restrict A, size_t N)
{
    double sum = 0;
#pragma omp parallel for shared(A) reduction(+:sum)
    for(size_t kk = 0; kk<N; kk++)
    {
        sum+=(double) A[kk];
    }

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

void fim_div(afloat * restrict  A,
             const afloat * restrict B,
             const afloat * restrict C,
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



void fim_minus(afloat * restrict  A,
               const afloat * restrict B,
               const afloat * restrict C,
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


float fim_max(const afloat * A, size_t N)
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
#pragma omp parallel for reduction(+:mse) shared(A, B)
    for(size_t kk = 0; kk<N; kk++)
    {
        mse += pow(A[kk]-B[kk], 2);
    }
    return mse/N;
}

void fim_flipall(afloat * restrict T, const afloat * restrict A, const int64_t a1, const int64_t a2, const int64_t a3)
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

void fim_insert_ref(afloat * T, int64_t t1, int64_t t2, int64_t t3,
                    afloat * F, int64_t f1, int64_t f2, int64_t f3)
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

afloat * fim_subregion(const afloat * restrict A, const int64_t M, const int64_t N, const int64_t P, const int64_t m, const int64_t n, const int64_t p)
{
    ((void) P);

    /* Extract sub region starting at (0,0,0) */
    afloat * S = fftwf_malloc(m*n*p*sizeof(float));
    assert(S != NULL);
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
#pragma omp parallel for shared(I)
    for(size_t kk = 0; kk<N; kk++)
    {
        I[kk] -= min;
    }
}

void fim_mult_scalar(afloat * I, size_t N, float x)
{
#pragma omp parallel for shared(I)
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
#pragma omp parallel for shared(psf) reduction(+:psf_sum)
    for(size_t kk = 0; kk<pMNP; kk++)
    { psf_sum += psf[kk]; }
    //  printf("psf_sum: %f\n", psf_sum);
#pragma omp parallel for shared(psf)
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
    assert(A != NULL);
    memset(A, 0, N*sizeof(float));
    return A;
}

afloat * fim_constant(const size_t N, const float value)
// Allocate and return an array of N floats sets to a constant value
{
    afloat * A = fftwf_malloc(N*sizeof(float));
#pragma omp parallel for shared(A)
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



    int nThreads = 0;

    /* Start a parallel region to figure out how many threads that will be used
       and how much memory to allocate for the buffers */
#pragma omp parallel
    {
        nThreads = omp_get_num_threads();
    }

    const size_t bsize = fmax(fmax(M, N), P);
    afloat * restrict buf = malloc(bsize*sizeof(float)*nThreads);


    /* Dimension 1 */
#pragma omp parallel for
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
#pragma omp parallel for
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
#pragma omp parallel for
    for(int64_t bb = 0; bb<N; bb++)
    {
        float * tbuf = buf + bsize*omp_get_thread_num();
        for(int64_t aa = 0; aa<M; aa++)
        {
            //shift_vector(A + aa+bb*M, M*N, P, sp);
            shift_vector_buf(A + aa+bb*M, M*N, P, sp, tbuf);
        }
    }

    free(buf);

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
    float * K = malloc(3*sizeof(float));
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
void fim_shift(afloat * restrict A,
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
    afloat * restrict buf = malloc(bsize*sizeof(float)*nThreads);


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

    free(buf);

    free(kernelx);
    free(kernely);
    free(kernelz);

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


INLINED static int64_t mod_int(const int64_t a, const int64_t b)
{
    int64_t r = a % b;
    return r < 0 ? r + b : r;
}


/* Shift vector by interpolation */
void shift_vector_float_buf(afloat * restrict V, // data
                            const int64_t S, // stride
                            const int64_t N, // elements
                            int n, // integer shift
                            afloat * restrict kernel, // centered kernel used for sub pixels shift
                            const int nkernel, // kernel size (odd!)
                            afloat * restrict buffer)
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
        conv1(V, S, buffer, N, kernel, nkernel);
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


void shift_vector_buf(afloat * restrict V, // data
                      const int64_t S, // stride
                      const int64_t N, // elements
                      int64_t k, // shift
                      afloat * restrict buffer)
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
    fftwf_free(V);
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
    fftwf_free(A);
    fftwf_free(L);

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
    fftwf_free(A);
}

void fim_xcorr2_ut()
{
    size_t M = 5;
    size_t N = 6;
    float * T = fim_zeros(M*N);
    float * A = fim_zeros(M*N);
    T[0] = 1;
    T[1] = 1;
    A[M+1] = 1;
    float * X = fim_xcorr2(T, A, M, N);
    printf("T=\n");
    fim_show(T, M, N, 1);
    printf("A=\n");
    fim_show(A, M, N, 1);
    printf("xcorr2(T,A)=\n");
    fim_show(X, 2*M-1, 2*N-1, 1);
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
    show_vec(K, N);
    free(K);
    sigma = 0.5;
    K = gaussian_kernel(sigma, &N);
    printf("gaussian_kernel, sigma=%f\n", sigma);
    show_vec(K, N);

    int nV = 10;
    float * V = malloc(nV*sizeof(float));
    for(int kk = 0; kk<nV; kk++)
    {
        V[kk] = kk;
    }
    printf("V=");
    show_vec(V, nV);
    conv1(V, 1, NULL, nV, K, N);
    printf("V*K = ");
    show_vec(V, nV);

    free(V);
    free(K);

    fim_cumsum_ut();
    fim_local_sum_ut();
    //exit(EXIT_FAILURE);
    myfftw_start(1, 1, stdout);
    fim_xcorr2_ut();
    myfftw_stop();
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

static void conv1(float * restrict V, int stride, float * restrict W,
                  const size_t nV,
                  const float * restrict K, const size_t nKu)
{
    /* Normalized convolution
     * In MATLAB that would be
     * Y = convn(V, K, 'same') / convn(ones(size(V)), K, 'same')
     */
    int Walloc = 0;
    if(W == NULL)
    {
        W = malloc(nV*sizeof(float));
        Walloc = 1;
    }
    const size_t k2 = (nKu-1)/2;
    const size_t N = nV;
    size_t bpos = 0;

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
        W[bpos++] = acc0/kacc;
    }

    /* Central part where K fits completely */
    for(size_t vv = k2 ; vv+k2 < N; vv++)
    {
        double acc = 0;
        for(size_t kk = 0; kk<nKu; kk++)
        {
            acc = acc + K[kk]*V[(vv-k2+kk)*stride];
        }
        W[bpos++] = acc;
    }

    /* Last part */
    for(size_t vv = N-k2;vv<N; vv++)
    {
        double kacc = 0;
        double acc0 = 0;
        for(size_t kk = 0; kk<N-vv+k2; kk++)
        {
            acc0 = acc0 + K[kk]*V[(vv-k2+kk)*stride];
            kacc += K[kk];
        }
        W[bpos++] = acc0/kacc;
    }

    /* Write back */
    for(size_t pp = 0; pp<nV; pp++)
    {
        V[pp*stride] = W[pp];
    }

    if(Walloc)
    {
        free(W);
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
    //printf("fim_gsmooth: M: %zu, N: %zu, P: %zu, nW: %zu\n", M, N, P, nW); fflush(stdout);

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
    fftwf_free(B);
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
    fftwf_free(C);
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
    fftwf_free(temp);
    temp = NULL;

    fftwf_complex * fxT = fft(xA, wM, wN, 1);
    fftwf_free(xA);
    fftwf_complex * fxA = fft(xT, wM, wN, 1);
    fftwf_free(xT);
    float * C = fft_convolve_cc(fxT, fxA, wM, wN, 1);
    // TODO fft_convolve_cc_conj instead?
    fftwf_free(fxT);
    fftwf_free(fxA);

    //printf("C=\n");
    //fim_show(C, 2*M-1, 2*N-1, 1);

    float * LS = fim_local_sum(A, M, N, M, N);
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

    float * LS2 = fim_local_sum(A2, M, N, M, N);

    float * denom_I = fim_zeros(wM*wN);
#pragma omp parallel for shared(LS2, denom_I)
    for(size_t kk = 0; kk<wM*wN; kk++)
    {
        float dLS = ( LS2[kk] - powf(LS[kk],2)/MN );
        float v = sqrt(dLS);
        v < 0 ? v = 0 : 0;
        denom_I[kk] = v;
    }
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

float * fim_maxproj(const float * V, size_t M, size_t N, size_t P)
{
    float * Pr = fim_zeros(M*N);
#pragma parallel for shared(Pr, A)
    for(size_t mm = 0; mm<M ; mm++)
    {
        for(size_t nn = 0; nn<N; nn++)
        {
            Pr[mm+M*nn] = fim_array_max(V+mm+M*nn, P, M*N);
        }
    }
    return Pr;
}
