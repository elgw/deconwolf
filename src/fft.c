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
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>
#include <string.h>
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "fim.h"

/* This provides some utility functions for using fftw3.
 *
 * Please see the fftw3 alignment requirements before making changes.
 */


/*
 * Settings
 */

/* Should only be changed with fft_set_plan() */
static unsigned int FFTW3_PLANNING = FFTW_MEASURE;
static int fft_nthreads = 1;
/* Enable with fft_set_inplace() */
static int use_inplace = 0;


fftwf_plan plan_r2c = NULL;
fftwf_plan plan_c2r = NULL;
fftwf_plan plan_r2c_inplace = NULL;
fftwf_plan plan_c2r_inplace = NULL;


/*
 * Forward declarations
 */
/* Pad data for inplace FFT transforms  */
static void fft_inplace_pad(float ** pX,
                     const size_t M,
                     const size_t N,
                     const size_t P);

/* Reverse the effect of fft_inplace_pad  */
static void fft_inplace_unpad(float ** pX,
                              const size_t M,
                              const size_t N,
                              const size_t P);

/* Number of complex floats required for the
 * Hermitian representation  */
static size_t nch(size_t M, size_t N, size_t P);

/** @brief Pad for inplace FFT
 *
 * Pad the input array along the first dimension to make room for
 * the inplace FFT. Example: If the input is:
 * [ ----- | ----- | ---- | ... ,where there is a '|' every M elements
 * then they array is converted to
 * [ -----P | -----P | ----P | ... ,where P is the padding
 *
 */
static void fft_inplace_pad(float ** pX,
                     const size_t M,
                     const size_t N,
                     const size_t P)
{

    const size_t nchunk = N*P;
    const size_t chunk_size = M*sizeof(float);

    *pX = fim_realloc(*pX, 2*nch(M, N, P) * sizeof(float));
    assert(*pX != NULL);

    // Think twice before trying to parallelize!

    if(M%2 == 0)
    {
        for(size_t c = nchunk-1; c != (size_t) -1; c--)
        {
            memmove(*pX+c*(M+2), *pX + c*M, chunk_size);
        }
    } else {
        for(size_t c = nchunk-1; c != (size_t) -1; c--)
        {
            memmove(*pX+c*(M+1), *pX + c*M, chunk_size);
        }
    }

    return;
}

/** @brief Reverse the effect of fft_inplace_pad
 *
*/
static void fft_inplace_unpad(float ** pX,
                              const size_t M,
                              const size_t N,
                              const size_t P)
{

    const size_t nchunk = N*P;
    const size_t chunk_size = M*sizeof(float);
    // Think twice before trying to parallelize!

    if(M%2 == 0)
    {
        for(size_t c = 0; c < nchunk; c++)
        {
            memmove(*pX + c*M, *pX+c*(M+2), chunk_size);
        }
    } else {
        for(size_t c = 0; c < nchunk; c++)
        {
            memmove(*pX + c*M, *pX+c*(M+1), chunk_size);
        }
    }

    *pX = fim_realloc(*pX, M*N*P * sizeof(float));
    assert(*pX != NULL);

    return;
}




static size_t nch(size_t M, size_t N, size_t P)
{
    //return ((3+M)/2)*N*P;
    return (1+M/2)*N*P;
}


#ifdef WINDOWS
#include <direct.h>
#include <sysinfoapi.h>
/* Not tested but it compiles */
int clock_gettime(int a , struct timespec *spec)      //C-file part
{  __int64 wintime; GetSystemTimeAsFileTime((FILETIME*)&wintime);
    wintime      -=116444736000000000;  //1jan1601 to 1jan1970
    spec->tv_sec  =wintime / 10000000;           //seconds
    spec->tv_nsec =wintime % 10000000 *100;      //nano-seconds
    return 0;
}
#endif

static char * get_swf_file_name(int nThreads)
{
    char * dir_home = getenv("HOME");
    char * dir_config = malloc(1024*sizeof(char));
    assert(dir_config != NULL);
    char * swf = malloc(1024*sizeof(char));
    assert(swf != NULL);

    if(use_inplace == 0)
    {
        sprintf(swf, "fftw_wisdom_float_threads_%d.dat", nThreads);
    } else {
        sprintf(swf, "fftw_wisdom_float_inplace_threads_%d.dat", nThreads);
    }

    sprintf(dir_config, "%s/.config/", dir_home);
    // printf("dir_config = %s\n", dir_config);
    if( !isdir(dir_config) )
    {
        free(dir_config);
        return swf;
    }

    sprintf(dir_config, "%s/.config/deconwolf/", dir_home);
    //  printf("dir_config = %s\n", dir_config);
    if( ensuredir(dir_config) == 0 )
    {
        char * prefered = malloc(1024*sizeof(char));
        assert(prefered != NULL);
        sprintf(prefered, "%s%s", dir_config, swf);
        free(dir_config);
        free(swf);
        return prefered;
    } else {
        free(dir_config);
        return swf;
    }
}

#ifdef CUDA
void myfftw_start(__attribute__((unused)) const int nThreads,
                  __attribute__((unused)) int verbose,
                  __attribute__((unused)) FILE * log)
{
    return;
}
#endif

void fft_set_plan(unsigned int plan)
{

    int ok = 0;
    switch(plan)
    {
    case FFTW_ESTIMATE:
        ok = 1;
        break;
    case FFTW_MEASURE:
        ok = 1;
        break;
    case FFTW_PATIENT:
        ok = 1;
        break;
    case FFTW_EXHAUSTIVE:
        ok = 1;
        break;
    default:
        ;
    }
    if(ok)
    {
        FFTW3_PLANNING = plan |  FFTW_UNALIGNED;
    } else {
        fprintf(stderr, "ERROR: Unknown FFTW plan, please use FFTW_ESTIMATE, "
                "FFTW_MEASURE, FFTW_PATIENT or FFTW_EXHAUSTIVE. %s %s %d\n",
                __FILE__, __FUNCTION__, __LINE__);
    }
    return;
}

#ifndef CUDA
void myfftw_start(const int nThreads, int verbose, FILE * log)
{
    fft_nthreads = nThreads;

    if(verbose > 1)
    {
#ifndef CUDA
        printf("\t using %s with %d threads\n", fftwf_version, nThreads);
#else
        printf("\t using cuFFT\n");
#endif
    }
    if(log != NULL)
    {
#ifndef CUDA
        fprintf(log, "Using %s with %d threads\n", fftwf_version, nThreads);
#else
        printf(log, "\t using cuFFT\n");
#endif
    }

    fftwf_init_threads();
    fftwf_plan_with_nthreads(nThreads);


    char * swf = get_swf_file_name(nThreads);
    assert(swf != NULL);
    if(log != NULL)
    {
        fprintf(log, "FFTW wisdom file: %s\n", swf);
    }
    if(swf == NULL)
    {
        assert(0);
    }
    else
    {
        if(log != NULL)
        {
            fprintf(log, "Importing FFTW wisdom\n");
        }
        fftwf_import_wisdom_from_filename(swf);
        free(swf);
    }
}
#endif

void myfftw_stop(void)
{
#ifndef CUDA
    fftwf_destroy_plan(plan_r2c);
    fftwf_destroy_plan(plan_c2r);
    fftwf_destroy_plan(plan_r2c_inplace);
    fftwf_destroy_plan(plan_c2r_inplace);
    fftwf_cleanup_threads();
#endif
    fftwf_cleanup();
    // Note: wisdom is only exported by fft_train
}

float * ifft(const fftwf_complex * fX, size_t M, size_t N, size_t P)
{

    float * X = fim_malloc(M*N*P*sizeof(float));
    assert(X != NULL);

    fftwf_execute_dft_c2r(plan_c2r, (fftwf_complex*) fX, X);

#pragma omp parallel for shared(X)
    for(size_t kk = 0 ; kk < M*N*P; kk++)
    {
        X[kk] /= (float) (M*N*P);
    }

    return X;
}


fftwf_complex * fft(const float * restrict in,
                    const int n1, const int n2, const int n3)
{
    assert(in != NULL);
    assert(plan_r2c != NULL);
    size_t N = nch(n1, n2, n3);
    fftwf_complex * out = fim_malloc(N*sizeof(fftwf_complex));
    assert(out != NULL);
    memset(out, 0, N*sizeof(fftwf_complex));

    fftwf_execute_dft_r2c(plan_r2c, (float*) in, out);

    return out;
}

void fft_mul(fftwf_complex * restrict C,
             fftwf_complex * restrict A,
             fftwf_complex * restrict B,
             const size_t n1, const size_t n2, const size_t n3)
{
    size_t N = nch(n1, n2, n3);
    /* C = A*B */
#pragma omp parallel for shared(A,B,C)
    for(size_t kk = 0; kk<N; kk++)
    {
        float a = A[kk][0]; float ac = A[kk][1];
        float b = B[kk][0]; float bc = B[kk][1];
        C[kk][0] = a*b - ac*bc;
        C[kk][1] = a*bc + b*ac;
    }
    return;
}

void fft_mul_inplace(fftwf_complex * restrict A,
                     fftwf_complex * restrict B,
                     const size_t n1, const size_t n2, const size_t n3)
{
    size_t N = nch(n1, n2, n3);
    /* C = A*B */
#pragma omp parallel for shared(A,B)
    for(size_t kk = 0; kk<N; kk++)
    {
        float a = A[kk][0]; float ac = A[kk][1];
        float b = B[kk][0]; float bc = B[kk][1];
        B[kk][0] = a*b - ac*bc;
        B[kk][1] = a*bc + b*ac;
    }
    return;
}



void fft_mul_conj(fftwf_complex * restrict C,
                  fftwf_complex * restrict A,
                  fftwf_complex * restrict B,
                  const size_t n1, const size_t n2, const size_t n3)
/* Multiply and conjugate the elements in the array A
 * i.e. C = conj(A)*B
 * All inputs should have the same size [n1 x n2 x n3]
 * */
{
    size_t N = nch(n1, n2, n3);
    // C = A*B
    size_t kk = 0;
#pragma omp parallel for shared(A, B, C)
    for(kk = 0; kk<N; kk++)
    {
        float a = A[kk][0]; float ac = -A[kk][1];
        float b = B[kk][0]; float bc = B[kk][1];
        C[kk][0] = a*b - ac*bc;
        C[kk][1] = a*bc + b*ac;
    }
    return;
}

float * fft_convolve_cc_f2(fftwf_complex * A, fftwf_complex * B,
                           const int M, const int N, const int P)
{
    fft_mul_inplace(A, B, M, N, P);
    float * out = ifft_and_free(B, M, N, P);
    return out;
}

void fft_mul_conj_inplace(fftwf_complex * restrict A,
                          fftwf_complex * restrict B,
                          const size_t n1, const size_t n2, const size_t n3)
/* Multiply and conjugate the elements in the array A
 * i.e. C = conj(A)*B
 * All inputs should have the same size [n1 x n2 x n3]
 * */
{
    //size_t N = nch(n1, n2, n3);
    size_t N = nch(n1, n2, n3);
    // C = A*B
    size_t kk = 0;
#pragma omp parallel for shared(A, B)
    for(kk = 0; kk<N; kk++)
    {
        float a = A[kk][0]; float ac = -A[kk][1];
        float b = B[kk][0]; float bc = B[kk][1];
        B[kk][0] = a*b - ac*bc;
        B[kk][1] = a*bc + b*ac;
    }
    return;
}


float * fft_convolve_cc_conj_f2(fftwf_complex * A, fftwf_complex * B,
                                const int M, const int N, const int P)
{
    fft_mul_conj_inplace(A, B, M, N, P);
    float * out = ifft_and_free(B, M, N, P);
    return out;
}


float * fft_convolve_cc(fftwf_complex * A, fftwf_complex * B,
                        const int M, const int N, const int P)
{
    size_t n = nch(M, N, P);
    fftwf_complex * C = fim_malloc(n*sizeof(fftwf_complex));
    assert(C != NULL);
    fft_mul(C, A, B, M, N, P);

    float * out = fim_malloc(M*N*P*sizeof(float));
    assert(out != NULL);

    fftwf_execute_dft_c2r(plan_r2c, C, out);
    free(C);

    const size_t MNP = M*N*P;
#pragma omp parallel for shared(out)
    for(size_t kk = 0; kk<MNP; kk++)
    {
        out[kk]/=(MNP);
    }
    return out;
}

float * fft_convolve_cc_conj(fftwf_complex * A, fftwf_complex * B,
                             const int M, const int N, const int P)
{
    size_t n = nch(M, N, P);
    fftwf_complex * C = fim_malloc(n*sizeof(fftwf_complex));
    assert(C != NULL);
    fft_mul_conj(C, A, B, M, N, P);

    float * out = fim_malloc(M*N*P*sizeof(float));
    assert(out != NULL);
    assert(plan_c2r != NULL);

    fftwf_execute_dft_c2r(plan_c2r, C, out);
    free(C);

    const size_t MNP = M*N*P;
#pragma omp parallel for shared(out)
    for(size_t kk = 0; kk<MNP; kk++)
    {
        out[kk]/=(MNP);
    }
    return out;
}

#ifdef CUDA
void fft_train( __attribute__((unused)) const size_t M,
                __attribute__((unused)) const size_t N,
                __attribute__((unused)) const size_t P,
                __attribute__((unused)) const int verbosity,
                __attribute__((unused)) const int nThreads,
                __attribute__((unused)) FILE * log)
{
    return;
}
#endif

#ifndef CUDA
void fft_train(const size_t M, const size_t N, const size_t P,
               const int verbosity, int nThreads,
               FILE * log)
{
    int updatedWisdom = 0;

    if(nThreads < 1)
    {
        nThreads = fft_nthreads;
    }
    nThreads < 1 ? nThreads = 1 : 0;

    if(verbosity > 0){
        printf("creating fftw3 plans ... \n"); fflush(stdout);
    }
    if(log != stdout)
    {
        fprintf(log, "--- fftw3 training ---\n");
    }

    /* Free old plans if exist */
    fftwf_destroy_plan(plan_r2c);
    fftwf_destroy_plan(plan_c2r);
    fftwf_destroy_plan(plan_r2c_inplace);
    fftwf_destroy_plan(plan_c2r_inplace);

    fftwf_complex * C = fim_malloc(nch(M, N, P)*sizeof(fftwf_complex));
    assert(C != NULL);
    float * R = fim_malloc(nch(M,N,P)*2*sizeof(float));
    assert(R != NULL);

    /* Hermitian to Real */

    plan_c2r = fftwf_plan_dft_c2r_3d(P, N, M,
                                     C, R, FFTW3_PLANNING  | FFTW_WISDOM_ONLY);
    if(plan_c2r == NULL)
    {
        plan_c2r = fftwf_plan_dft_c2r_3d(P, N, M,
                                         C, R, FFTW3_PLANNING);
        printf("   c2r plan ... \n");
        fftwf_execute(plan_c2r);
        updatedWisdom = 1;
    }

    plan_c2r_inplace = fftwf_plan_dft_c2r_3d(P, N, M,
                                             C, (float *) C, FFTW3_PLANNING  | FFTW_WISDOM_ONLY);
    if(plan_c2r_inplace == NULL)
    {
        plan_c2r_inplace = fftwf_plan_dft_c2r_3d(P, N, M,
                                                 C, (float *) C, FFTW3_PLANNING);
        printf("   c2r inplace plan ... \n");
        fftwf_execute(plan_c2r_inplace);
        updatedWisdom = 1;
    }

    /* Real to Hermitian */

    plan_r2c = fftwf_plan_dft_r2c_3d(P, N, M,
                                     R, C,
                                     FFTW3_PLANNING | FFTW_WISDOM_ONLY);
    if(plan_r2c == NULL)
    {
        plan_r2c = fftwf_plan_dft_r2c_3d(P, N, M,
                                         R, C,
                                         FFTW3_PLANNING);
        printf("   r2c plan ... \n");
        fftwf_execute(plan_r2c);
        updatedWisdom = 1;
    }
    plan_r2c_inplace = fftwf_plan_dft_r2c_3d(P, N, M,
                                             (float*) C, C,
                                             FFTW3_PLANNING | FFTW_WISDOM_ONLY);
    if(plan_r2c_inplace == NULL)
    {
        plan_r2c_inplace = fftwf_plan_dft_r2c_3d(P, N, M,
                                                 (float*) C, C,
                                                 FFTW3_PLANNING);
        printf("   r2c inplace plan ... \n");
        fftwf_execute(plan_r2c_inplace);
        updatedWisdom = 1;
    }


    free(C);
    free(R);

    if(updatedWisdom)
    {
        char * swf = get_swf_file_name(nThreads);
        assert(swf != NULL);

        fprintf(log, "Exporting fftw wisdom to %s\n", swf);
        int ret = fftwf_export_wisdom_to_filename(swf);

        if(ret != 0)
        {
            if(verbosity > 1)
            {
                printf("Exported fftw wisdom to %s\n", swf);
            }
        } else {
            printf("ERROR; Failed to write fftw wisdom to %s\n", swf);
        }
        free(swf);
    }

    return;
}
#endif

void fft_ut_wisdom_name(void){
    /* Wisdom file names
     * Try this when $HOME/.config/deconwolf/ does not exist
     * and when it does ... jonas.paulsen1
     * could also test it when that dir isn't writeable.
     * */

    int nThreads = 4;
    char * swf = get_swf_file_name(nThreads);
    printf("swf = '%s'\n", swf);
    free(swf);
}

void fft_ut_flipall_conj()
{
    myfftw_start(2, 0, NULL);
    /* Test the identity
     * flip(X) = ifft(conj(fft(X)))
     */
    int M = 12, N = 13, P = 15;
    float * A = fim_malloc(M*N*P*sizeof(float));
    assert(A != NULL);
    for(int kk = 0; kk<M*N*P; kk++)
    { A[kk] = (float) rand() / (float) RAND_MAX; }
    fim_stats(A, M*N*P);
    float * B = fim_malloc(M*N*P*sizeof(float));
    assert(B != NULL);
    memcpy(B, A, M*N*P*sizeof(float));
    float * B_flipall = fim_malloc(M*N*P*sizeof(float));
    assert(B_flipall != NULL);
    fim_flipall(B_flipall, B, M, N, P);


    fftwf_complex * FA = fft(A, M, N, P);
    fftwf_complex * FB = fft(B, M, N, P);
    free(A);
    free(B);
    fftwf_complex * FB_flipall = fft(B_flipall, M, N, P);

    float * Y1 = fft_convolve_cc(FA, FB_flipall, M, N, P);
    float * Y2 = fft_convolve_cc_conj(FA, FB, M, N, P);
    free(FA);
    free(FB);

    float mse = fim_mse(Y1, Y2, M*N*P);
    printf("mse=%f ", mse);
    if(mse < 1e-5)
    { printf("ok!\n"); } else
    { printf("BAD :(\n"); }


    free(B_flipall);
    free(FB_flipall);

    free(Y1);
    free(Y2);

    myfftw_stop();
}


double * fft_bench_1d(int64_t from, int64_t to, int niter)
{
    assert(to >= from);
    assert(from > 0);
    assert(niter > 0);

    struct timespec tictoc_start, tictoc_end;

    double * t = fim_malloc(sizeof(double) * (to + 1 - from));
    assert( t != NULL);
    memset(t, 0, sizeof(double)*(to+1-from));
    for(int kk = 0; kk<niter; kk++)
    {
        for(int64_t s = from; s <= to; s++)
        {
            float * fft_in = fim_malloc(s*sizeof(float));
            assert(fft_in != NULL);
            for(int64_t kk = 0; kk<s; kk++)
            {
                fft_in[kk] = (float)  ((double) rand() / (double) RAND_MAX);
            }
            fftwf_complex * fft_out = fim_malloc(s * sizeof(fftwf_complex));
            assert(fft_out != NULL);
            fftwf_plan plan = fftwf_plan_dft_r2c_1d(s,
                                                    fft_in, fft_out,
                                                    FFTW_ESTIMATE);

            clock_gettime(CLOCK_REALTIME, &tictoc_start);
            fftwf_execute(plan);
            clock_gettime(CLOCK_REALTIME, &tictoc_end);

            double took = timespec_diff(&tictoc_end, &tictoc_start);

            if(kk > 0)
            {
                t[s-from] += took/(double) niter;
            }


            fftwf_destroy_plan(plan);
            free(fft_out);
            free(fft_in);
        }
    }
    return t;
}

/* Only used for fft_ut. Should be renamed to fim_max_rel_error */
static float fim_compare(const float * X, const float * Y,
                  size_t M, size_t N, size_t P)
{
    float rel_err_max = -1;
    float tol = 1e-5;
    for(size_t kk = 0 ; kk<M*N*P; kk++)
    {
        if( fabs(X[kk])> tol)
        {
            float rel_err = fabs( (X[kk] - Y[kk])/X[kk] );
            rel_err > rel_err_max ? rel_err_max = rel_err : 0;
        }
    }
    return rel_err_max;
}

float * test_data_rand(size_t M, size_t N, size_t P)
{
    srand(time(NULL));
    size_t MNP = M*N*P;
    float * X = fim_malloc(M*N*P*sizeof(float));
    assert(X != NULL);
    for(size_t kk = 0; kk<MNP; kk++)
    {
        X[kk] = 5 + (float) rand() / (float) RAND_MAX;
    }
    return X;
}

void test_inplace(void)
{
    printf(" -> Testing out-of-place vs in-place\n");
    for(size_t kk = 0; kk<1000; kk++)
    {
        /* Normal non-inplace version */
        size_t M = 51+(rand() % 10);
        size_t N = 25+(rand() % 10);
        size_t P = 20+(rand() % 10);
        fft_train(M, N, P, 0, 8, stdout);

        //M = 2228; N = 2228; P = 208; // caused strange problems
        printf("Test image size: %zu %zu %zu\n", M, N, P);
        size_t MNP = M*N*P;
        float * X = test_data_rand(M, N, P);
        printf("fft\n");
        fftwf_complex * fX = fft(X, M, N, P);
        printf("ifft\n");
        float * ffX = ifft(fX, M, N, P);
        float relerr_oo = fim_compare(X, ffX, M, N, P);
        free(ffX);
        printf("fft:out-of-place ifft:out-of-place max_rel_err: %e\n", relerr_oo);

        /* Only forward inplace */
        float * X2 = fim_malloc(M*N*P*sizeof(float));
        assert(X2 != NULL);
        memcpy(X2, X, M*N*P*sizeof(float));
        float * dummy1 = malloc(MNP);
        assert(dummy1 != NULL);
        printf("fft_inplace\n");
        fftwf_complex * fX2 = fft_inplace(X2, M, N, P);
        float * dummy2 = malloc(MNP);
        assert(dummy2 != NULL);
        printf("ifft\n");
        float * ffX2 = ifft(fX2, M , N, P);
        float relerr_io = fim_compare(X, ffX2, M, N, P);
        printf("fft:in-place ifft:out-of-place max_rel_err: %e\n", relerr_io);

        /* Forward inplace and backward inplace */
        float * X3 = fim_malloc(M*N*P*sizeof(float));
        assert(X3 != NULL);
        float * dummy3 = malloc(MNP);
        assert(dummy3 != NULL);
        memcpy(X3, X, M*N*P*sizeof(float));
        printf("fft_inplace\n");
        fftwf_complex * fX3 = fft_inplace(X3, M, N, P);
        float * dummy4 = malloc(MNP);
        assert(dummy4 != NULL);
        printf("ifft_inplace\n");
        float * ffX3 = ifft_inplace(fX3, M , N, P);
        float relerr_ii = fim_compare(X, ffX3, M, N, P);
        printf("fft:in-place ifft:out-of-place max_rel_err: %e\n", relerr_ii);

        assert(relerr_oo < 1e-5);
        assert(relerr_io < 1e-5);
        assert(relerr_ii < 1e-5);

        free(dummy1);
        free(dummy2);
        free(dummy3);
        free(dummy4);
        free(X);
        free(ffX3);
        free(fX2);
        free(ffX2);
        free(fX);
    }
    exit(EXIT_SUCCESS);
    return;
}

void fft_ut(void)
{
    fft_set_inplace(0);
    fft_set_plan(FFTW_ESTIMATE);
    myfftw_start(8, 1, stdout);
    fftwf_forget_wisdom();

    test_inplace();
    tictoc;
    tic;
    fft_ut_wisdom_name();
    fft_ut_flipall_conj();
    toc(fft_ut took);

    int64_t from = 1024;
    int64_t to = 1030;
    double * t = fft_bench_1d(from, to, 1000);
    for(int64_t kk = 0; kk<= (to-from); kk++)
    {
        printf("Size: %" PRId64 ", time %e\n", kk+from, t[kk]);
    }
    free(t);

    // Typical z-sizes:
    from = 80;
    to = 90;
    t = fft_bench_1d(from, to, 10000);
    for(int64_t kk = 0; kk<= (to-from); kk++)
    {
        printf("Size: %" PRId64 ", time %e\n", kk+from, t[kk]);
    }
    free(t);



    // Free plans etc
    myfftw_stop();
    return;
}

fftwf_complex * fft_and_free(float * restrict in,
                             const int n1, const int n2, const int n3)
{
    if(use_inplace == 1)
    {
        fftwf_complex * F = fft_inplace(in, n1, n2, n3);
        return F;
    } else {
        fftwf_complex * F = fft(in, n1, n2, n3);
        free(in);
        return F;
    }
}

float * ifft_and_free(fftwf_complex * F,
                      const size_t n1, const size_t n2, const size_t n3)
{
    if(use_inplace == 1)
    {
        float * f = ifft_inplace(F, n1, n2, n3);
        return f;
    } else {
        float * f = ifft(F, n1, n2, n3);
        free(F);
        return f;
    }
}

void fft_set_inplace(int ip)
{
    if(ip == 1)
    {
        use_inplace = 1;
        return;
    }
    if(ip == 0)
    {
        use_inplace = 0;
        return;
    }

    fprintf(stderr,
            "WARNING, %s %s %d, specified value(%d) is "
            "not valid, use either 0 or 1\n",
            __FILE__, __FUNCTION__, __LINE__, ip);
    return;
}

fftwf_complex * fft_inplace(float * X, const size_t M, const size_t N, const size_t P)
{

    fft_inplace_pad(&X, M, N, P);
    assert(plan_r2c_inplace != NULL);
    fftwf_execute_dft_r2c(plan_r2c_inplace, X, (fftwf_complex *) X);
    return (fftwf_complex*) X;
}

float * ifft_inplace(fftwf_complex * fX, const size_t M, const size_t N, const size_t P)
{
    float * X = (float *) fX;

    assert(plan_c2r_inplace != NULL);
    fftwf_execute_dft_c2r(plan_c2r_inplace, fX, (float *) X);

    fft_inplace_unpad(&X, M, N, P);
#pragma omp parallel for shared(X)
    for(size_t kk = 0 ; kk < M*N*P; kk++)
    {
        X[kk] /= (float) (M*N*P);
    }
    return X;
}
