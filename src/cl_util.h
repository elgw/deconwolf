#pragma once

/* Only a single clu should be used at a time.
 * One clu is only set up for a specific fft size.
 *
 * For vkFFT:
 * - Todo: Forward and backward FFT (inplace and/or out of place)
 *
 * For clFFT:
 * - Fix so that the inplace transforms work! (default is out-of-place, right?)
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "dw_util.h"

/* VkFFT targets version 1.20 */
#define CL_TARGET_OPENCL_VERSION 120

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif


#define CLU_KEEP_ALL 0
#define CLU_DROP_ALL 1
#define CLU_KEEP_2ND 2


#ifdef VKFFT_BACKEND /* Should be defined in the makefile */
// Unfortunately this does not turn off all warnings

#include <vkFFT.h>

#ifndef VKFFT
#define VKFFT
#endif
#else
#include <clFFT.h>
#endif

#include "fim.h"

/* Pad first dimension */
#define PAD_FIRST_DIM 1

const char* clGetErrorString(int errorCode);

typedef struct{
    cl_kernel kernel;
    cl_program program;
} clu_kernel_t;


typedef struct{
    int verbose;
    cl_context context;
    cl_platform_id platform_id;
    cl_device_id device_id;
    cl_command_queue command_queue;

    size_t M; size_t N; size_t P;
    /* For complex data */
    clu_kernel_t * kern_mul;
    clu_kernel_t * kern_mul_conj;
    clu_kernel_t * kern_mul_inplace;
    clu_kernel_t * kern_mul_conj_inplace;
    /* For real data */
    clu_kernel_t * kern_real_mul_inplace;
    clu_kernel_t * kern_real_positivity;
    /* mem buffer with 1 float. Used for positivity threshold and for alpha */
    cl_mem float_gpu;
    clu_kernel_t * kern_shb_update; /* Find next guess with shb */
    clu_kernel_t * kern_preprocess_image;
    clu_kernel_t * idiv_kernel;
    clu_kernel_t * update_y_kernel;

    size_t nb_allocated;
    size_t n_release;
    size_t n_alloc;

    /* Prefer in-place transformations over out-of place?
     * in-place should use less memory, but is it slower? */
    int prefer_inplace;
#ifdef VKFFT
    VkFFTApplication vkfft_app;
#else
    int clFFT_loaded;
    clfftSetupData fftSetup;
    size_t clfft_buffer_size;
    cl_mem clfft_buffer; // Only allocated if clfft_buffer_size > 0
    clfftPlanHandle r2h_plan;
    clfftPlanHandle r2h_inplace_plan;
    clfftPlanHandle h2r_plan;
    clfftPlanHandle h2r_inplace_plan;
#endif
} clu_env_t;

/* A fimcl object can be one of these types */
typedef enum {
    fimcl_real, /* Real (padding on is applied internally) */
    fimcl_hermitian /* complex */
}
    fimcl_type ;

typedef struct{
    /* Size of real array before padding */
    size_t M;
    size_t N;
    size_t P;

    cl_mem buf;
    size_t buf_size_nf; // Number of floats, not bytes

    fimcl_type type;

    clu_env_t * clu;
    /* Check this before data is used, initialized to NULL */
    cl_event wait_ev;
} fimcl_t;

/*******************************************************
 *     Setup and teardown of environment
 *******************************************************/

/** Create an environment with OpenCL and clFFT
 * See also clu_destroy
 */
clu_env_t * clu_new(int verbose, int cl_device);

/* Prepare to do FFTs of real arrays of size wM x wN x wP
 * the size M x N x P is the size of the original image
 * without padding for Berteros method
 */
void clu_prepare_kernels(clu_env_t * clu,
                         size_t wM, size_t wN, size_t wP,
                         size_t M, size_t N, size_t P);

/* Tear down what is crated with clu_new */
void clu_destroy(clu_env_t * );



/* Allocate a new float image on the GPU.
 * Use data in X unless it is NULL
 */
fimcl_t * fimcl_new(clu_env_t * clu, fimcl_type type,
                    const float * X, size_t M, size_t N, size_t P);

void fimcl_free(fimcl_t * );



/* Get back data from GPU.
 * Blocking (sync before and after) */
float * fimcl_download(fimcl_t *);

/* C = A*B, where C will be allocated or possibly reused from A
 * flags can be
 * - CLU_KEEP_ALL -- A and B are left untouched
 * - CLU_DROP_ALL -- A and B are freed.
 *   Inplace transformation can be used if A->fullsize=1
 * - CLU_DROP_2ND, B will be dropped (reused)
 * Returns an object which has transformed == 0, i.e. real data
 */
fimcl_t * fimcl_convolve(fimcl_t * A, fimcl_t * B, int flags);
fimcl_t * fimcl_convolve_conj(fimcl_t * A, fimcl_t * B, int flags);

/* Assumes that the input is real.
 * The output will be complex hermitian.
 * A = fft(A) */
fimcl_t * fimcl_fft(fimcl_t * A);
void fimcl_fft_inplace(fimcl_t * X);
fimcl_t * fimcl_ifft(fimcl_t * A);
void fimcl_ifft_inplace(fimcl_t * A);

/* B = copy(A)
 * non-blocking, sync object before using */
fimcl_t * fimcl_copy(fimcl_t * );

/* fZ = fX.*fY
 * set conj to 1 for complex conjugate
 * wait event on fZ */
void fimcl_complex_mul(fimcl_t * fX, fimcl_t * fY, fimcl_t * fZ, int conj);

/* fY = fX .* fY */
void fimcl_complex_mul_inplace(fimcl_t * fX, fimcl_t * fY, int conj);

/* Y = X.*Y */
void fimcl_real_mul_inplace(fimcl_t * X, fimcl_t * Y);

/*  P[idx] = x[idx] + alpha[0]*(x[idx]-xp[idx]) */
void fimcl_shb_update(fimcl_t * P, fimcl_t * X, fimcl_t * XP, float alpha);

/* X[kk] < val ? X[kk] = val : 0 */
void fimcl_positivity(fimcl_t * X, float val);

/* Wait until the object is available.
 * To be used after fft, ifft, complex_mul etc ...*/
void fimcl_sync(fimcl_t * X);

/* Try to resolve the error code and then exit */
void clu_exit_error(cl_int err,
                    const char * file,
                    const char * function,
                    int line,
                    int clfft);




/* host to host deconvolution via OpenCl
 * in case dropy is set, Y is freed (or reused for the output)
 * Uses in-place ffts
 */
float * clu_convolve(clu_env_t * clu,
                     int dropy,
                     float * X, float * Y,
                     size_t M, size_t N, size_t P);

#ifndef VKFFT
const char * get_clfft_error_string(clfftStatus error);
#endif

/* Load a program either from a file (if file_name != NULL)
 * or from a string (if program_code != NULL)
 * if both program_code and file_name is supplied, code is loaded
 * from the file.
 * program_size must be set when program_code is supplied
 * Returns NULL on failure
 */
clu_kernel_t * clu_kernel_new(clu_env_t * env,
                              const char * file_name,
                              const char * program_code,
                              size_t program_size,
                              const char * kernel_name);

/* As clu_kernel_new but also forward arguments to the compiler
 * useful for -D PI=3.14 etc. */
clu_kernel_t * clu_kernel_newa(clu_env_t * env,
                               const char * file,
                               const char * program_code,
                               size_t program_size,
                               const char * kernel_name,
                               const char * argument_string);

void clu_kernel_destroy(clu_kernel_t * kern);


/* replacement for the bussy-waiting clWaitForEvents
   check fimcl_sync for usage.
*/
cl_uint clu_wait_for_event(cl_event clev, size_t ns);

/* Print some of the device information available */
void clu_print_device_info(FILE *, cl_device_id dev_id);

/* Return the smallest integer >= N that has no
 * prime factors > 13. These are the only sizes supported
 * by clFFT
 */
size_t clu_next_fft_size(size_t N);

#ifndef VKFFT
/* Generate a plan from real to complex hermitian */
clfftPlanHandle  gen_r2h_plan(clu_env_t * clu,
                              size_t M, size_t N, size_t P);
clfftPlanHandle  gen_r2h_inplace_plan(clu_env_t * clu,
                                      size_t M, size_t N, size_t P);
/* Generate a plan from complex hermitian to real */
clfftPlanHandle  gen_h2r_plan(clu_env_t * clu,
                              size_t M, size_t N, size_t P);
clfftPlanHandle  gen_h2r_inplace_plan(clu_env_t * clu,
                                      size_t M, size_t N, size_t P);


/* increase the size of the clFFT buffer to size. Don't do anything if
 * the current buffer is already large enough */
cl_int clu_increase_clfft_buffer(clu_env_t * clu, size_t size);
#endif

void clu_benchmark_transfer(clu_env_t * clu);

float fimcl_error_idiv(fimcl_t * forward, fimcl_t * image);
void fimcl_update_y(fimcl_t * gy, fimcl_t * image);

void fimcl_preprocess( fimcl_t * fft_image, fimcl_t * fft_PSF, float value);

/* Pad the data and make it ready for in-place transformations
 * This is called by fimcl_new when type is set to fimcl_real_inplace
 */
float * pad_for_inplace(const float * X,
                        size_t M, size_t N, size_t P);

/* Reverse of unpad_from_inplace, used with fimcl_download when
 * X->type == fimcl_real_inplace */
float * unpad_from_inplace(const float * PX,
                           size_t M, size_t N, size_t P);
