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

#define tictoc struct timespec tictoc_start, tictoc_end;
#define tic clock_gettime(CLOCK_REALTIME, &tictoc_start);
#define toc(X) clock_gettime(CLOCK_REALTIME, &tictoc_end); printf(#X); printf(" %f s\n", timespec_diff(&tictoc_end, &tictoc_start));

static double timespec_diff(struct timespec* end, struct timespec * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}


/* Should only be changed with fft_set_plan */
unsigned int FFTW3_PLANNING = FFTW_MEASURE;

static int fft_nthreads;

static int isdir(char * dir)
{
    /* Check if directory exist, do not create if missing
     * returns 1 if it exist
     * */

    DIR* odir = opendir(dir);
    if (odir) {
        /* Directory exists. */
        closedir(odir);
        return 1;
    } else if (ENOENT == errno) {
        /* Directory does not exist. */
        return 0;
    } else {
        /* opendir() failed for some other reason. */
        return 0;
    }
}

static int ensuredir(char * dir)
/* Create dir if it does not exist.
 * Returns 0 if the dir already existed or could be created
 * returns non-zeros if the dir can't be created
 */
{
    if(isdir(dir) == 1)
    {
        return 0;
    }

#ifdef WINDOWS
    if(_mkdir(dir) == ENOENT)
    {
        return 0;
    }
#else
    if(mkdir(dir, 0700) == 0)
    {
        return 0;
    }
#endif

    return 1;
}

static char * get_swf_file_name(int nThreads)
{
    char * dir_home = getenv("HOME");
    char * dir_config = malloc(1024*sizeof(char));
    char * swf = malloc(1024*sizeof(char));
    sprintf(swf, "fftw_wisdom_float_threads_%d.dat", nThreads);

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
        FFTW3_PLANNING = plan;
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
    fftwf_cleanup_threads();
    #endif
    fftwf_cleanup();
    // Note: wisdom is only exported by fft_train
}


fftwf_complex * fft(float * restrict in, const int n1, const int n2, const int n3)
{
    size_t N = (n1+3)/2*n2*n3;
    fftwf_complex * out = fftwf_malloc(N*sizeof(fftwf_complex));
    memset(out, 0, N*sizeof(fftwf_complex));

    fftwf_plan p = fftwf_plan_dft_r2c_3d(n3, n2, n1,
                                         in, // Float
                                         out, // fftwf_complex
                                         FFTW3_PLANNING);
    fftwf_execute(p);
    fftwf_destroy_plan(p);
    return out;
}

void fft_mul(fftwf_complex * restrict C,
             fftwf_complex * restrict A,
             fftwf_complex * restrict B,
             const size_t n1, const size_t n2, const size_t n3)
{
    size_t N = (n1+3)/2*n2*n3;
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


void fft_mul_conj(fftwf_complex * restrict C,
                  fftwf_complex * restrict A,
                  fftwf_complex * restrict B,
                  const size_t n1, const size_t n2, const size_t n3)
/* Multiply and conjugate the elements in the array A
 * i.e. C = conj(A)*B
 * All inputs should have the same size [n1 x n2 x n3]
 * */
{
    size_t N = (n1+3)/2*n2*n3;
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
    size_t n = (M+3)/2*N*P;
    fftwf_complex * C = fftwf_malloc(n*sizeof(fftwf_complex));
    fft_mul(C, A, B, M, N, P);
    fftwf_free(B);

    float * out = fftwf_malloc(M*N*P*sizeof(float));

    fftwf_plan p = fftwf_plan_dft_c2r_3d(P, N, M,
                                         C, out,
                                         FFTW3_PLANNING);
    fftwf_execute(p);
    fftwf_destroy_plan(p);
    fftwf_free(C);

    const size_t MNP = M*N*P;

#pragma omp parallel for shared(out)
    for(size_t kk = 0; kk < MNP; kk++)
    {
        out[kk]/=(MNP);
    }

    return out;
}

float * fft_convolve_cc_conj_f2(fftwf_complex * A, fftwf_complex * B,
                                 const int M, const int N, const int P)
{

    size_t n = (M+3)/2*N*P;
    fftwf_complex * C = fftwf_malloc(n*sizeof(fftwf_complex));

    fft_mul_conj(C, A, B, M, N, P);

    // fft_mul(C, A, B, M, N, P); /* For symmetric PSFs this would give the same */
    fftwf_free(B);

    float * out = fftwf_malloc(M*N*P*sizeof(float));

    fftwf_plan p = fftwf_plan_dft_c2r_3d(P, N, M,
                                         C, out,
                                         FFTW3_PLANNING);
    fftwf_execute(p);
    fftwf_destroy_plan(p);
    fftwf_free(C);
    const float MNP_float = (float) M*N*P;
    const size_t MNP = M*N*P;
#pragma omp parallel for shared(out)
    for(size_t kk = 0; kk < MNP; kk++)
    {
        out[kk] /= MNP_float;
    }
    return out;
}



float * fft_convolve_cc(fftwf_complex * A, fftwf_complex * B,
                         const int M, const int N, const int P)
{
    size_t n = (M+3)/2*N*P;
    fftwf_complex * C = fftwf_malloc(n*sizeof(fftwf_complex));
    fft_mul(C, A, B, M, N, P);

    float * out = fftwf_malloc(M*N*P*sizeof(float));

    fftwf_plan p = fftwf_plan_dft_c2r_3d(P, N, M,
                                         C, out,
                                         FFTW3_PLANNING);
    fftwf_execute(p);
    fftwf_destroy_plan(p);
    fftwf_free(C);

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
    size_t n = (M+3)/2*N*P;
    fftwf_complex * C = fftwf_malloc(n*sizeof(fftwf_complex));
    fft_mul_conj(C, A, B, M, N, P);

    float * out = fftwf_malloc(M*N*P*sizeof(float));

    fftwf_plan p = fftwf_plan_dft_c2r_3d(P, N, M,
                                         C, out,
                                         FFTW3_PLANNING);
    fftwf_execute(p);
    fftwf_destroy_plan(p);
    fftwf_free(C);

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

    if(FFTW3_PLANNING == FFTW_ESTIMATE)
    {
        if(verbosity > 0)
        {
            printf("No training needed for FFTW_ESTIMATE\n");
        }
        fprintf(log, "No training needed for FFTW_ESTIMATE\n");
        return;
    }

    if(verbosity > 1){
        printf("fftw3 training ... \n"); fflush(stdout);
    }
    if(log != stdout)
    {
        fprintf(log, "--- fftw3 training ---\n");
    }

    size_t MNP = M*N*P;
    fftwf_complex * C = fftwf_malloc(MNP*sizeof(fftwf_complex));
    float * R = fftwf_malloc(MNP*sizeof(float));

    fftwf_plan p0 = fftwf_plan_dft_c2r_3d(P, N, M,
                                          C, R, FFTW3_PLANNING  | FFTW_WISDOM_ONLY);

    if(p0 == NULL)
    {
        if(verbosity > 0)
        {
            printf("> generating fftw3 c2r plan\n");
        }
        fprintf(log, "> generating fftw3 c2r plan\n");
        fftwf_plan p1 = fftwf_plan_dft_c2r_3d(P, N, M,
                                              C, R, FFTW3_PLANNING);
        fftwf_execute(p1);
        fftwf_destroy_plan(p1);
        updatedWisdom = 1;
    } else {
        if(verbosity > 1)
        {
            printf("\t Using cached fftw3 c2r plan\n");
        }
        if(log != stdout)
        {
            fprintf(log, "Using cached fftw3 c2r plan\n");
        }
    }
    fftwf_destroy_plan(p0);

    p0 = fftwf_plan_dft_r2c_3d(P, N, M,
                               R, C, FFTW3_PLANNING | FFTW_WISDOM_ONLY);

    if(p0 == NULL)
    {
        if(verbosity > 0){
            printf("> generating fftw3 r2c plan \n");
        }
        if(log != stdout)
        {
            fprintf(log, "> generating fftw3 r2c plan\n");
        }
        fftwf_plan p2 = fftwf_plan_dft_r2c_3d(P, N, M,
                                              R, C, FFTW3_PLANNING);
        fftwf_execute(p2);
        fftwf_destroy_plan(p2);
        updatedWisdom = 1;
    } else {
        if(verbosity > 1)
        {
            printf("\t Using cached fftw3 r2c plan\n");
        }
        fprintf(log, "Using cached fftw3 r2c plan\n");
    }
    fftwf_destroy_plan(p0);

    fftwf_free(R);
    fftwf_free(C);

    if(updatedWisdom)
    {
        char * swf = get_swf_file_name(nThreads);
        if(swf == NULL)
        { assert(0); }
        else {
            fprintf(log, "Exporting fftw wisdom to %s\n", swf);
            int ret = fftwf_export_wisdom_to_filename(swf);
            if(verbosity > 0)
            {
                if(ret != 0)
                {
                    printf("Exported fftw wisdom to %s\n", swf);
                } else {
                    printf("ERROR; Failed to write fftw wisdom to %s\n", swf);
                }

            }
            free(swf);
        }
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
    float * A = fftwf_malloc(M*N*P*sizeof(float));
    for(int kk = 0; kk<M*N*P; kk++)
    { A[kk] = (float) rand() / (float) RAND_MAX; }
    fim_stats(A, M*N*P);
    float * B = fftwf_malloc(M*N*P*sizeof(float));
    memcpy(B, A, M*N*P*sizeof(float));
    float * B_flipall = fftwf_malloc(M*N*P*sizeof(float));
    fim_flipall(B_flipall, B, M, N, P);


    fftwf_complex * FA = fft(A, M, N, P);
    fftwf_complex * FB = fft(B, M, N, P);
    fftwf_complex * FB_flipall = fft(B_flipall, M, N, P);

    float * Y1 = fft_convolve_cc(FA, FB_flipall, M, N, P);
    float * Y2 = fft_convolve_cc_conj(FA, FB, M, N, P);

    float mse = fim_mse(Y1, Y2, M*N*P);
    printf("mse=%f ", mse);
    if(mse < 1e-5)
    { printf("ok!\n"); } else
    { printf("BAD :(\n"); }

    fftwf_free(A); fftwf_free(FA);
    fftwf_free(B); fftwf_free(FB);
    fftwf_free(B_flipall); fftwf_free(FB_flipall);

    fftwf_free(Y1);
    fftwf_free(Y2);

    myfftw_stop();
}


double * fft_bench_1d(int64_t from, int64_t to, int niter)
{
    assert(to >= from);
    assert(from > 0);
    assert(niter > 0);

    struct timespec tictoc_start, tictoc_end;

    double * t = malloc(sizeof(double) * (to + 1 - from));
    memset(t, 0, sizeof(double)*(to+1-from));
    for(int kk = 0; kk<niter; kk++)
    {
        for(int64_t s = from; s <= to; s++)
        {
            float * fft_in = fftwf_malloc(s*sizeof(float));
            for(int64_t kk = 0; kk<s; kk++)
            {
                fft_in[kk] = (float)  ((double) rand() / (double) RAND_MAX);
            }
            fftwf_complex * fft_out = fftwf_malloc(s * sizeof(fftwf_complex));
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
            fftwf_free(fft_out);
            fftwf_free(fft_in);
        }
    }
    return t;
}

void fft_ut(void)
{
    tictoc
        tic
        fft_ut_wisdom_name();
    fft_ut_flipall_conj();
    toc(fft_ut took)

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
    fftwf_cleanup();
    return;
}
