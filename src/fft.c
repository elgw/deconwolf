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

#include <inttypes.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "fim.h"
#include "dw_util.h"

typedef int64_t i64;

/*
 * Forward declarations
 */
/* Pad data for inplace FFT transforms  */
static void fft_inplace_pad(dw_fft *,
                            float ** pX,
                            const size_t M,
                            const size_t N,
                            const size_t P);

/* Reverse the effect of fft_inplace_pad  */
static void fft_inplace_unpad(dw_fft * ff,
                              float ** pX,
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
static void fft_inplace_pad(dw_fft * ff,
                            float ** pX,
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
static void fft_inplace_unpad(dw_fft * ff,
                              float ** pX,
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

static char * get_swf_file_name(dw_fft * fft)
{
    char * dir_home = getenv("HOME");
    char * dir_config = malloc(1024*sizeof(char));
    assert(dir_config != NULL);
    char * swf = malloc(1024*sizeof(char));
    assert(swf != NULL);

    if(fft->use_inplace == 0)
    {
        sprintf(swf, "fftw_wisdom_float_threads_%d.dat", fft->n_thread);
    } else {
        sprintf(swf, "fftw_wisdom_float_inplace_threads_%d.dat", fft->n_thread);
    }

    sprintf(dir_config, "%s/.config/", dir_home);
    // printf("dir_config = %s\n", dir_config);
    if( !dw_isdir(dir_config) )
    {
        free(dir_config);
        return swf;
    }

    sprintf(dir_config, "%s/.config/deconwolf/", dir_home);
    //  printf("dir_config = %s\n", dir_config);
    if( dw_ensuredir(dir_config) == 0 )
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

void fft_set_plan(dw_fft * fft, unsigned int plan)
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
        fft->FFTW3_PLANNING = plan |  FFTW_UNALIGNED;
    } else {
        fprintf(stderr, "ERROR: Unknown FFTW plan, please use FFTW_ESTIMATE, "
                "FFTW_MEASURE, FFTW_PATIENT or FFTW_EXHAUSTIVE. %s %s %d\n",
                __FILE__, __FUNCTION__, __LINE__);
    }
    return;
}

static void fft_train(dw_fft * ff)
{
    assert(ff != NULL);
    if(ff->plan_r2c != NULL)
    {
        printf("fft_train error: expected plan_r2c to be NULL\n");
        exit(EXIT_FAILURE);
    }
    int updatedWisdom = 0;

    if(ff->verbose > 0){
        printf("creating fftw3 plans ... \n"); fflush(stdout);
    }
    if(ff->log != NULL)
    {
        fprintf(ff->log, "--- fftw3 training ---\n");
    }

    i64 M = ff->M;
    i64 N = ff->N;
    i64 P = ff->P;

    fftwf_complex * C = fim_malloc(nch(M, N, P)*sizeof(fftwf_complex));
    assert(C != NULL);
    float * R = fim_malloc(nch(M,N,P)*2*sizeof(float));
    assert(R != NULL);

    /* Hermitian to Real */

    // Se if we have one stored?
    ff->plan_c2r = fftwf_plan_dft_c2r_3d(P, N, M,
                                     C, R, ff->FFTW3_PLANNING  | FFTW_WISDOM_ONLY);
    if(ff->plan_c2r == NULL)
    {
        //
        ff->plan_c2r = fftwf_plan_dft_c2r_3d(P, N, M,
                                         C, R, ff->FFTW3_PLANNING);

        printf("   c2r plan ..."); fflush(stdout);
        fftwf_execute(ff->plan_c2r);
        printf("\n");

        updatedWisdom = 1;
    }

    ff->plan_c2r_inplace = fftwf_plan_dft_c2r_3d(P, N, M,
                                             C, (float *) C, ff->FFTW3_PLANNING  | FFTW_WISDOM_ONLY);
    if(ff->plan_c2r_inplace == NULL)
    {
        ff->plan_c2r_inplace = fftwf_plan_dft_c2r_3d(P, N, M,
                                                 C, (float *) C, ff->FFTW3_PLANNING);
        printf("   c2r inplace plan ..."); fflush(stdout);
        fftwf_execute(ff->plan_c2r_inplace);
        printf("\n");
        updatedWisdom = 1;
    }

    /* Real to Hermitian */

    ff->plan_r2c = fftwf_plan_dft_r2c_3d(P, N, M,
                                     R, C,
                                     ff->FFTW3_PLANNING | FFTW_WISDOM_ONLY);
    if(ff->plan_r2c == NULL)
    {
        ff->plan_r2c = fftwf_plan_dft_r2c_3d(P, N, M,
                                         R, C,
                                         ff->FFTW3_PLANNING);
        printf("   r2c plan ..."); fflush(stdout);
        fftwf_execute(ff->plan_r2c);
        printf("\n");
        updatedWisdom = 1;
    }
    ff->plan_r2c_inplace = fftwf_plan_dft_r2c_3d(P, N, M,
                                             (float*) C, C,
                                             ff->FFTW3_PLANNING | FFTW_WISDOM_ONLY);
    if(ff->plan_r2c_inplace == NULL)
    {
        ff->plan_r2c_inplace = fftwf_plan_dft_r2c_3d(P, N, M,
                                                 (float*) C, C,
                                                 ff->FFTW3_PLANNING);
        printf("   r2c inplace plan ..."); fflush(stdout);
        fftwf_execute(ff->plan_r2c_inplace);
        printf("\n");
        updatedWisdom = 1;
    }

    fim_free(C);
    fim_free(R);

    if(updatedWisdom)
    {
        if(ff->verbose > 1)
        {
            printf("Exporting fftw wisdom to %s\n", ff->wisdom_file);
        }
        if(ff->log != NULL)
        {
            fprintf(ff->log, "Exporting fftw wisdom to %s\n", ff->wisdom_file);
        }
        int ret = fftwf_export_wisdom_to_filename(ff->wisdom_file);

        if(ret != 0)
        {
            if(ff->verbose > 1)
            {
                printf("Exported fftw wisdom to %s\n", ff->wisdom_file);
            }
        } else {
            printf("ERROR; Failed to write fftw wisdom to %s\n", ff->wisdom_file);
        }
    }

    return;
}


dw_fft * dw_fft_new(int n_thread, int verbose, FILE * log,
                    i64 M, i64 N, i64 P,
    int planner_flags)
{
    dw_fft * ff = calloc(1, sizeof(dw_fft));
    ff->use_inplace = 1;
    ff->n_thread = n_thread;
    ff->verbose = verbose;
    ff->FFTW3_PLANNING = planner_flags;
    ff->M = M;
    ff->N = N;
    ff->P = P;

    fftwf_init_threads();
    fftwf_plan_with_nthreads(n_thread);

    /* Wisdom loading */
    ff->wisdom_file = get_swf_file_name(ff);
    assert(ff->wisdom_file != NULL);
    if(log != NULL)
    {
        fprintf(log, "FFTW wisdom file: %s\n", ff->wisdom_file);
    }
    if(ff->verbose > 1)
    {
        printf("FFTW wisdom file: %s\n", ff->wisdom_file);
    }


    if(log != NULL)
    {
        fprintf(log, "Importing FFTW wisdom\n");
    }
    fftwf_import_wisdom_from_filename(ff->wisdom_file);

    fft_train(ff);
    return ff;
}


void dw_fft_destroy(dw_fft * fft)
{
    if(fft == NULL)
    {
        printf("Warning calling dw_fft_destroy(NULL)\n");
        return;
    }
    fftwf_destroy_plan(fft->plan_r2c);
    fft->plan_r2c = NULL;
    fftwf_destroy_plan(fft->plan_c2r);
    fft->plan_c2r = NULL;
    fftwf_destroy_plan(fft->plan_r2c_inplace);
    fft->plan_r2c_inplace = NULL;
    fftwf_destroy_plan(fft->plan_c2r_inplace);
    fft->plan_c2r_inplace = NULL;
    fftwf_cleanup_threads();
    fftwf_cleanup();
    free(fft->wisdom_file);
    free(fft);
}

float * ifft(dw_fft * fft, const fftwf_complex * fX, size_t M, size_t N, size_t P)
{

    float * X = fim_malloc(M*N*P*sizeof(float));
    assert(X != NULL);

    fftwf_execute_dft_c2r(fft->plan_c2r, (fftwf_complex*) fX, X);

#pragma omp parallel for shared(X)
    for(size_t kk = 0 ; kk < M*N*P; kk++)
    {
        X[kk] /= (float) (M*N*P);
    }

    return X;
}


fftwf_complex * fft(dw_fft * fft, const float * restrict in,
                    const int n1, const int n2, const int n3)
{
    assert(in != NULL);
    assert(fft != NULL);
    assert(fft->plan_r2c != NULL);
    size_t N = nch(n1, n2, n3);
    fftwf_complex * out = fim_malloc(N*sizeof(fftwf_complex));
    assert(out != NULL);
    memset(out, 0, N*sizeof(fftwf_complex));

    fftwf_execute_dft_r2c(fft->plan_r2c, (float*) in, out);

    return out;
}

void fft_mul(dw_fft * fft,
             fftwf_complex * restrict C,
             fftwf_complex * restrict A,
             fftwf_complex * restrict B,
             const size_t n1, const size_t n2, const size_t n3)
{
    assert(fft != NULL);
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

void fft_mul_inplace(dw_fft * fft, fftwf_complex * restrict A,
                     fftwf_complex * restrict B,
                     const size_t n1, const size_t n2, const size_t n3)
{
    assert(fft != NULL);
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



void fft_mul_conj(dw_fft * fft,
                  fftwf_complex * restrict C,
                  fftwf_complex * restrict A,
                  fftwf_complex * restrict B,
                  const size_t n1, const size_t n2, const size_t n3)
/* Multiply and conjugate the elements in the array A
 * i.e. C = conj(A)*B
 * All inputs should have the same size [n1 x n2 x n3]
 * */
{
    assert(fft != NULL);
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

float * fft_convolve_cc_f2(dw_fft * fft, fftwf_complex * A, fftwf_complex * B,
                           const int M, const int N, const int P)
{
    fft_mul_inplace(fft, A, B, M, N, P);
    float * out = ifft_and_free(fft, B, M, N, P);
    return out;
}

void fft_mul_conj_inplace(dw_fft * fft,
                          fftwf_complex * restrict A,
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


float * fft_convolve_cc_conj_f2(dw_fft * fft, fftwf_complex * A, fftwf_complex * B,
                                const int M, const int N, const int P)
{
    fft_mul_conj_inplace(fft, A, B, M, N, P);
    float * out = ifft_and_free(fft, B, M, N, P);
    return out;
}


float * fft_convolve_cc(dw_fft * fft, fftwf_complex * A, fftwf_complex * B,
                        const int M, const int N, const int P)
{
    size_t n = nch(M, N, P);
    fftwf_complex * C = fim_malloc(n*sizeof(fftwf_complex));
    assert(C != NULL);
    fft_mul(fft, C, A, B, M, N, P);

    float * out = fim_malloc(M*N*P*sizeof(float));
    assert(out != NULL);

    fftwf_execute_dft_c2r(fft->plan_r2c, C, out);
    fim_free(C);

    const i64 MNP = (i64) M*N*P;
    const float fMNP = (float) M*N*P;
#pragma omp parallel for shared(out)
    for(size_t kk = 0; kk<MNP; kk++)
    {
        out[kk] /= fMNP;
    }
    return out;
}

float * fft_convolve_cc_conj(dw_fft * fft,
                             fftwf_complex * A, fftwf_complex * B,
                             const int M, const int N, const int P)
{
    size_t n = nch(M, N, P);
    fftwf_complex * C = fim_malloc(n*sizeof(fftwf_complex));
    assert(C != NULL);
    fft_mul_conj(fft, C, A, B, M, N, P);

    float * out = fim_malloc(M*N*P*sizeof(float));
    assert(out != NULL);
    assert(fft->plan_c2r != NULL);

    fftwf_execute_dft_c2r(fft->plan_c2r, C, out);
    fim_free(C);

    const size_t MNP = M*N*P;
#pragma omp parallel for shared(out)
    for(size_t kk = 0; kk<MNP; kk++)
    {
        out[kk]/=(MNP);
    }
    return out;
}




void fft_ut_flipall_conj()
{

    /* Test the identity
     * flip(X) = ifft(conj(fft(X)))
     */
    int M = 12, N = 13, P = 15;
    dw_fft * ff = dw_fft_new(2, 1, NULL, M, N, P, FFTW_ESTIMATE);
    float * A = fim_malloc(M*N*P*sizeof(float));

    for(int kk = 0; kk<M*N*P; kk++)
    { A[kk] = (float) rand() / (float) RAND_MAX; }
    fim_stats(A, M*N*P);
    float * B = fim_malloc(M*N*P*sizeof(float));

    memcpy(B, A, M*N*P*sizeof(float));
    float * B_flipall = fim_malloc(M*N*P*sizeof(float));

    fim_flipall(B_flipall, B, M, N, P);


    fftwf_complex * FA = fft(ff, A, M, N, P);
    fim_free(A);
    fftwf_complex * FB = fft(ff, B, M, N, P);
    fim_free(B);

    fftwf_complex * FB_flipall = fft(ff, B_flipall, M, N, P);

    float * Y1 = fft_convolve_cc(ff, FA, FB_flipall, M, N, P);

    float * Y2 = fft_convolve_cc_conj(ff, FA, FB, M, N, P);
    assert(FA != FB);
    fim_free(FA);
    FA = NULL;

    fim_free(FB);
    FB = NULL;

    float mse = fim_mse(Y1, Y2, M*N*P);
    printf("mse=%f ", mse);
    if(mse < 1e-5)
    { printf("ok!\n"); } else
    { printf("BAD :(\n"); }


    fim_free(B_flipall);
    fim_free(FB_flipall);

    fim_free(Y1);
    fim_free(Y2);

    dw_fft_destroy(ff);
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

            dw_gettime(&tictoc_start);
            fftwf_execute(plan);
            dw_gettime(&tictoc_end);

            double took = timespec_diff(&tictoc_end, &tictoc_start);

            if(kk > 0)
            {
                t[s-from] += took/(double) niter;
            }


            fftwf_destroy_plan(plan);
            fim_free(fft_out);
            fim_free(fft_in);
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
        dw_fft * ff = dw_fft_new(2, 0, NULL, M, N, P, FFTW_ESTIMATE);

        fft_train(ff);

        //M = 2228; N = 2228; P = 208; // caused strange problems
        printf("Test image size: %zu %zu %zu\n", M, N, P);
        size_t MNP = M*N*P;
        float * X = test_data_rand(M, N, P);
        printf("fft\n");
        fftwf_complex * fX = fft(ff, X, M, N, P);
        printf("ifft\n");
        float * ffX = ifft(ff, fX, M, N, P);
        float relerr_oo = fim_compare(X, ffX, M, N, P);
        fim_free(ffX);
        printf("fft:out-of-place ifft:out-of-place max_rel_err: %e\n", relerr_oo);

        /* Only forward inplace */
        float * X2 = fim_malloc(M*N*P*sizeof(float));
        assert(X2 != NULL);
        memcpy(X2, X, M*N*P*sizeof(float));
        float * dummy1 = fim_malloc(MNP);
        assert(dummy1 != NULL);
        printf("fft_inplace\n");
        fftwf_complex * fX2 = fft_inplace(ff, X2, M, N, P);
        float * dummy2 = fim_malloc(MNP);
        assert(dummy2 != NULL);
        printf("ifft\n");
        float * ffX2 = ifft(ff, fX2, M , N, P);
        float relerr_io = fim_compare(X, ffX2, M, N, P);
        printf("fft:in-place ifft:out-of-place max_rel_err: %e\n", relerr_io);

        /* Forward inplace and backward inplace */
        float * X3 = fim_malloc(M*N*P*sizeof(float));
        assert(X3 != NULL);
        float * dummy3 = fim_malloc(MNP);
        assert(dummy3 != NULL);
        memcpy(X3, X, M*N*P*sizeof(float));
        printf("fft_inplace\n");
        fftwf_complex * fX3 = fft_inplace(ff, X3, M, N, P);
        float * dummy4 = fim_malloc(MNP);
        assert(dummy4 != NULL);
        printf("ifft_inplace\n");
        float * ffX3 = ifft_inplace(ff, fX3, M , N, P);
        float relerr_ii = fim_compare(X, ffX3, M, N, P);
        printf("fft:in-place ifft:out-of-place max_rel_err: %e\n", relerr_ii);

        assert(relerr_oo < 1e-5);
        assert(relerr_io < 1e-5);
        assert(relerr_ii < 1e-5);

        fim_free(dummy1);
        fim_free(dummy2);
        fim_free(dummy3);
        fim_free(dummy4);
        fim_free(X);
        fim_free(ffX3);
        fim_free(fX2);
        fim_free(ffX2);
        fim_free(fX);
        dw_fft_destroy(ff);
    }
    exit(EXIT_SUCCESS);
    return;
}

void fft_ut(void)
{
    test_inplace();
    tictoc;
    tic;
    fft_ut_flipall_conj();
    toc(fft_ut took);

    int64_t from = 1024;
    int64_t to = 1030;
    double * t = fft_bench_1d(from, to, 1000);
    for(int64_t kk = 0; kk<= (to-from); kk++)
    {
        printf("Size: %" PRId64 ", time %e\n", kk+from, t[kk]);
    }
    fim_free(t);

    // Typical z-sizes:
    from = 80;
    to = 90;
    t = fft_bench_1d(from, to, 10000);
    for(int64_t kk = 0; kk<= (to-from); kk++)
    {
        printf("Size: %" PRId64 ", time %e\n", kk+from, t[kk]);
    }
    fim_free(t);

    return;
}

fftwf_complex * fft_and_free(dw_fft * ff, float * restrict in,
                             const int n1, const int n2, const int n3)
{
    if(ff->use_inplace == 1)
    {
        fftwf_complex * F = fft_inplace(ff, in, n1, n2, n3);
        return F;
    } else {
        fftwf_complex * F = fft(ff, in, n1, n2, n3);
        fim_free(in);
        return F;
    }
}

float * ifft_and_free(dw_fft * ff, fftwf_complex * F,
                      const size_t n1, const size_t n2, const size_t n3)
{
    if(ff->use_inplace == 1)
    {
        float * f = ifft_inplace(ff, F, n1, n2, n3);
        return f;
    } else {
        float * f = ifft(ff, F, n1, n2, n3);
        fim_free(F);
        return f;
    }
}

fftwf_complex * fft_inplace(dw_fft * ff, float * X, const size_t M, const size_t N, const size_t P)
{
    fft_inplace_pad(ff, &X, M, N, P);
    assert(ff->plan_r2c_inplace != NULL);
    fftwf_execute_dft_r2c(ff->plan_r2c_inplace, X, (fftwf_complex *) X);
    return (fftwf_complex*) X;
}

float * ifft_inplace(dw_fft * ff, fftwf_complex * fX, const size_t M, const size_t N, const size_t P)
{
    float * X = (float *) fX;

    assert(ff->plan_c2r_inplace != NULL);
    fftwf_execute_dft_c2r(ff->plan_c2r_inplace, fX, (float *) X);

    fft_inplace_unpad(ff, &X, M, N, P);
    float MNP = (float) M*N*P;
#pragma omp parallel for shared(X)
    for(size_t kk = 0 ; kk < M*N*P; kk++)
    {
        X[kk] /= MNP;
    }
    return X;
}
