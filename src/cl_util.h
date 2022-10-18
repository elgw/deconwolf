#ifndef __cl_util_h__
#define __cl_util_h__

/* This is tailored for a single FFT size. Doing it that way we can
 * do much initialization in bulk at the beginning.
 * Notes:
 * -
*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#define CLFFT_REQUEST_LIB_NOMEMALLOC
#define CL_TARGET_OPENCL_VERSION 300
#include <CL/cl.h>
#include <clFFT.h>

#define CLU_KEEP_ALL 0
#define CLU_DROP_ALL 1
#define CLU_KEEP_2ND 2

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

    int clFFT_loaded;
    clfftSetupData fftSetup;
    size_t clfft_buffer_size;
    cl_mem clfft_buffer; // Only allocated if clfft_buffer_size > 0
    clfftPlanHandle r2h_plan;
    clfftPlanHandle r2h_inplace_plan;
    clfftPlanHandle h2r_plan;
    clfftPlanHandle h2r_inplace_plan;
    size_t M; size_t N; size_t P;
    /* For complex data */
    clu_kernel_t kern_mul;
    clu_kernel_t kern_mul_conj;
    clu_kernel_t kern_mul_inplace;
    clu_kernel_t kern_mul_conj_inplace;
    /* For real data */
    cl_mem real_size;
    clu_kernel_t kern_real_mul_inplace;
    clu_kernel_t kern_error_idiv;
    clu_kernel_t kern_real_positivity;
    /* mem buffer with 1 float. Used for positivity threshold and for alpha */
    cl_mem float_gpu;
    clu_kernel_t kern_shb_update; /* Find next guess with shb */

    size_t nb_allocated;
    size_t n_release;
    size_t n_alloc;
} clu_env_t;


typedef struct{
    size_t M;
    size_t N;
    size_t P;
    cl_mem buf;
    size_t buf_size_nf; // Number of floats, not bytes
    int transformed; /* 0=real or 1=hermitian */
    int fullsize; /* If the buffer is large enough for in-place fft */
    int padded; /* Ready for inplace transform ? */
    clu_env_t * clu;
    cl_event wait_ev; // Check this before data is used
} fimcl_t;

/* Allocate a new float image on the GPU.
 * Use data in X unless it is NULL
 * if fullsize is set to 1, this object can be used for
 * inplace transformations
 * non-blocking. call fimcl_sync on the new object to manually sync
 */
fimcl_t * fimcl_new(clu_env_t * clu, int ffted, int fullsize,
                    const float * X, size_t M, size_t N, size_t P);

void fimcl_free(fimcl_t * );

/* Get back data from GPU. Blocking. */
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
 * non-blocking */
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


#define check_CL(x) if(x != CL_SUCCESS)                             \
    {                                                               \
        clu_exit_error(x, __FILE__, __FUNCTION__, __LINE__, 0);     \
    }                                                               \

#define check_clFFT(x) if(x != CL_SUCCESS)                          \
    {                                                               \
        clu_exit_error(x, __FILE__, __FUNCTION__, __LINE__, 1);     \
    }                                                               \


/* Try to resolve the error code and then exit */
void clu_exit_error(cl_int err,
                    const char * file,
                    const char * function,
                    int line,
                    int clfft);


/* Create an environment with OpenCL and clFFT */
clu_env_t * clu_new(int verbose);

/* Prepare to do FFTs */
void clu_prepare_fft(clu_env_t * clu,
                     size_t M, size_t N, size_t P);

/* Tear down what is crated with clu_new */
void clu_destroy(clu_env_t * );


/* host to host deconvolution via OpenCl
 * in case dropy is set, Y is freed (or reused for the output)
 * Uses in-place ffts
 */
float * clu_convolve(clu_env_t * clu,
                     int dropy,
                     float * X, float * Y,
                     size_t M, size_t N, size_t P);

const char * get_clfft_error_string(clfftStatus error);

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

void clu_kernel_destroy(clu_kernel_t kern);


/* replacement for clWaitForEvents which isn't bussy waiting.

   Unfortunately it looks like there is no way to signal a thread
   that a job is done.
   This is one way to avoid bussy waiting. See
   CL_QUEUE_THROTTLE_LOW_KHR for another alternative.
   https://github.com/intel/compute-runtime/issues/363
*/
cl_uint clu_wait_for_event(cl_event clev, size_t ns);

/* Print some of the device information available */
void clu_print_device_info(FILE *, cl_device_id dev_id);

/* Return the smallest integer >= N that has no
 * prime factors > 13. These are the only sizes supported
 * by clFFT
 */
size_t clu_next_fft_size(size_t N);

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

void clu_benchmark_transfer(clu_env_t * clu);
#endif
