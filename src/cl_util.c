#include "cl_util.h"
#include "clext.h"

#include "kernels/cl_complex_mul.h"
#include "kernels/cl_complex_mul_conj.h"
#include "kernels/cl_complex_mul_inplace.h"
#include "kernels/cl_complex_mul_conj_inplace.h"

#include "kernels/cl_real_mul_inplace.h"
#include "kernels/cl_positivity.h"

#include "kernels/cl_shb_update.h"
#include "kernels/cl_preprocess_image.h"
#include "kernels/cl_idiv_kernel.h"
#include "kernels/cl_update_y_kernel.h"

#define check_CL(x) if(x != CL_SUCCESS)                         \
    {                                                           \
        clu_exit_error(x, __FILE__, __FUNCTION__, __LINE__, 0); \
    }                                                           \

#define check_clFFT(x) if(x != CL_SUCCESS)                      \
    {                                                           \
        clu_exit_error(x, __FILE__, __FUNCTION__, __LINE__, 1); \
    }                                                           \


//#define here(x) printf("%s %s %d\n", __FILE__, __FUNCTION__, __LINE__); fflush(stdout);
#define here(x) ;

static char * read_program(const char * fname, size_t * size);
static double clockdiff(struct timespec* start, struct timespec * finish);

/* Number of complex values required for Hermitian representation */
static size_t fimcl_ncx(fimcl_t * X);

static size_t fimcl_hM(fimcl_t * X);
static size_t fimcl_hN(fimcl_t * X);
static size_t fimcl_hP(fimcl_t * X);

/* Number of complex numbers along the first dimension after FFT with
 * hermitian symmetry taken into considerations, based on the number
 * of real numbers.
 */
static size_t M_r2h(size_t M);

/* Number of floats in the first dimension required for in-place
 * transformation
 */
static size_t padsize(size_t M)
{
    if(M % 2 == 0)
    {
        return M+2;
    }
    return M+1;
}

/** Instead of bare clReleaseMemObject
 * since clReleaseMemObject only decrease the reference count,
 * we check that there is only one reference available, else
 * the object will not be freed and we have to figure out
 * what else uses the mem object  (or decrease the reference count again)
 * We must also check that clReleaseProgram is used and as a worst case
 * also call clReleaseContext
 * missing clReleaseEvent ? -- check all cl_event ...(wait_ev)
 * https://forums.developer.nvidia.com/t/clreleasememobject-not-working-the-clreleasememobject-function-doesnt-release-the-memory/15963/5
 */
static int clu_release(clu_env_t * C, cl_mem buf)
{
    cl_uint ref_count;
    check_CL(
             clGetMemObjectInfo(buf, CL_MEM_REFERENCE_COUNT,
                                sizeof(ref_count), &ref_count, NULL));
    if(ref_count != 1)
    {
        fprintf(stderr, "clu_release WARNING\n"
                " - Object has a reference count > 1\n");
    }
    check_CL(clReleaseMemObject(buf));
    C->n_release++;
}

void clu_exit_error(cl_int err,
                    const char * file,
                    const char * function,
                    int line,
                    int clfft)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Sorry! There was an unrecoverable error!\n"
            "   File: %s\n"
            "   Function: %s at line %d\n",
            file, function, line);

    const char * err_cl = clGetErrorString(err);
    printf("   OpenCl error=%s\n", err_cl);

#ifndef VKFFT_BACKEND
    if(clfft)
    {
        if(strcmp(err_cl, "CL_UNKNOWN_ERROR") == 0)
        {
            const char * err_clfft = get_clfft_error_string(err);
            fprintf(stderr, "clFFT error=%s\n", err_clfft);
        }
    }
#endif

    if(err == CL_MEM_OBJECT_ALLOCATION_FAILURE)
    {
        fprintf(stderr, "\n");
        fprintf(stderr,
                "   CL_MEM_OBJECT_ALLOCATION_FAILURE indicates that there was not\n"
                "   enough memory on the GPU to continue. Try with a smaller image\n"
                "   and look up the option --tilesize\n");
    }
    fprintf(stderr, "\n");
    fprintf(stderr,
            "   If you are sure that OpenCL works on this machine\n"
            "   and that it is a problem only related to deconwolf,\n"
            "   check open issues or create a new one at\n"
            "   https://github.com/elgw/deconwolf/issues\n");
    exit(EXIT_FAILURE);
}

/* Convert real dimensions to complex hermitian dimensions */
static size_t M_r2h(size_t M)
{
    return (1+M/2);
}

static size_t fimcl_hM(fimcl_t * X)
{
    return M_r2h(X->M);
}
static size_t fimcl_hN(fimcl_t * X)
{
    return X->N;
}

static size_t fimcl_hP(fimcl_t * X)
{
    return X->P;
}


static size_t fimcl_ncx(fimcl_t * X)
{
    return fimcl_hM(X)*fimcl_hN(X)*fimcl_hP(X);
}


float * fimcl_download(fimcl_t * gX)
{
    here();
    assert(gX != NULL);
    fimcl_sync(gX);
    if(gX->type == fimcl_real)
    {
        size_t MNP = gX->M*gX->N*gX->P;
        if(gX->clu->verbose > 1)
        {
            printf("Downloading real data %zu x %zu x %zu (%zu floats)\n",
                   gX->M, gX->N, gX->P, MNP);
        }

        size_t nfloat = padsize(gX->M)*gX->N*gX->P;
        float * PX = fim_malloc(nfloat*sizeof(float));

        if(PX == NULL)
        {
            fprintf(stderr, "Failed to allocate memory in %s %s %d\n", __FILE__,
                    __FUNCTION__,
                    __LINE__);
            fprintf(stderr, "Tried to allocate %zu x %zu x %zu (%zu floats)\n",
                    gX->M, gX->N, gX->P, MNP);
            exit(EXIT_FAILURE);
        }

        check_CL( clEnqueueReadBuffer( gX->clu->command_queue,
                                       gX->buf,
                                       CL_TRUE, // blocking
                                       0, // offset
                                       nfloat * sizeof( float ), // size
                                       PX,
                                       0,
                                       NULL,
                                       &gX->wait_ev ));
        fimcl_sync(gX);

        float * X = unpad_from_inplace(PX, gX->M, gX->N, gX->P);
        fim_free(PX);
        return X;
    }

    if(gX->type == fimcl_hermitian)
    {
        size_t cMNP = fimcl_ncx(gX);
        if(gX->clu->verbose > 1)
        {
            printf("Downloading transformed data %zu x %zu x %zu (%zu floats)\n",
                   fimcl_hM(gX), fimcl_hN(gX), fimcl_hP(gX), cMNP*2);
        }

        float * X = fim_malloc(cMNP*2*sizeof(float));
        if(X == NULL)
        {
            fprintf(stderr, "Failed to allocate memory in %s/%s/%d\n", __FILE__,
                    __FUNCTION__,
                    __LINE__);
            exit(EXIT_FAILURE);
        }

        check_CL( clEnqueueReadBuffer( gX->clu->command_queue,
                                       gX->buf,
                                       CL_TRUE, // blocking
                                       0, // offset
                                       cMNP *2* sizeof( float ), // size
                                       X,
                                       0,
                                       NULL,
                                       &gX->wait_ev ));
        fimcl_sync(gX);

        return X;
    }
    assert(0);
    return NULL;
}

fimcl_t * fimcl_new(clu_env_t * clu, fimcl_type type,
                    const float * X, size_t M, size_t N, size_t P)
{

    if(clu->verbose > 2)
    {
        printf("fimcl_new, type=");
        switch(type)
        {
        case fimcl_real:
            printf("real\n");
            break;
        case fimcl_hermitian:
            printf("fimcl_hermitian\n");
            break;
        }
    }

    fimcl_t * Y = calloc(1, sizeof(fimcl_t));
    Y->M = M;
    Y->N = N;
    Y->P = P;
    Y->type = type;
    Y->clu = clu;
    Y->wait_ev = NULL;

    cl_int ret;
    if(clu->verbose > 2)
    {
        printf("Creating cl buffer\n");
    }


    Y->buf_size_nf = fimcl_ncx(Y)*2;

    Y->buf = clCreateBuffer(clu->context,
                            CL_MEM_READ_WRITE,
                            Y->buf_size_nf*sizeof(float),
                            NULL, &ret );
    clu->nb_allocated += Y->buf_size_nf*sizeof(float);
    clu->n_alloc++;

    check_CL(ret);

    if(X == NULL)
    {
        goto done;
    }

    if(clu->verbose > 2)
    {
        printf("Uploading data to GPU\n");
    }


    if(type == fimcl_real)
    {
        //printf("fimcl_ncx(Y)*2: %zu\n", fimcl_ncx(Y)*2);
        float * PX = pad_for_inplace(X, M, N, P);
        check_CL( clEnqueueWriteBuffer( clu->command_queue,
                                        Y->buf,
                                        CL_TRUE, // blocking_write
                                        0,
                                        fimcl_ncx(Y)*2*sizeof( float ),
                                        PX,
                                        0, // num_events_in_wait_list
                                        NULL,
                                        &Y->wait_ev ) );
        fimcl_sync(Y);

        fim_free(PX);
    }

    if(type == fimcl_hermitian)
    {
        fprintf(stderr, "TODO: Unable to upload already ffted data in %s, line %d\n",
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    if(clu->verbose > 2)
    {
        printf("Upload complete\n");
    }

 done: ;
    return Y;
}

fimcl_t * fimcl_copy(fimcl_t * G)
{
    if(G->clu->verbose > 1)
    {
        printf("fimcl_copy\n");
    }

    fimcl_t * H = fimcl_new(G->clu, G->type, NULL, G->M, G->N, G->P);

    check_CL(clEnqueueCopyBuffer(
                                 G->clu->command_queue,//cl_command_queue command_queue,
                                 G->buf, //cl_mem src_buffer,
                                 H->buf, //cl_mem dst_buffer,
                                 0, //size_t src_offset,
                                 0, //size_t dst_offset,
                                 padsize(G->M)*G->N*G->P*sizeof(float),//size_t size,
                                 0, //cl_uint num_events_in_wait_list,
                                 NULL, //const cl_event* event_wait_list,
                                 &H->wait_ev)); //cl_event* event);
    fimcl_sync(H);

    H->buf_size_nf = G->buf_size_nf;
    assert(H->buf_size_nf == padsize(G->M)*G->N*G->P);
    return H;
}

void fimcl_free(fimcl_t * G)
{
    if(G == NULL)
    {
        fprintf(stderr, "Warning: fimcl_free(NULL);\n");
        return;
    }

    fimcl_sync(G);
    clu_release(G->clu, G->buf);
    free(G);
}


static fimcl_t * _fimcl_convolve(fimcl_t * X, fimcl_t * Y, int mode, int conj)
{
    assert(X != NULL);
    assert(Y != NULL);


    if(X->clu->verbose > 1)
    {
        printf("fimcl_convolve mode=%d, conj=%d\n", mode, conj);
    }

    if(X->type == fimcl_hermitian && Y->type == fimcl_hermitian)
    {

        if(mode == CLU_KEEP_ALL)
        {
            fimcl_t * Z = fimcl_new(X->clu, fimcl_hermitian,
                                    NULL, X->M, X->N, X->P);
            fimcl_complex_mul(X, Y, Z, conj);
            fimcl_ifft(Z);
            return Z;
        }

        if(mode == CLU_DROP_ALL)
        {
            fimcl_complex_mul_inplace(Y, Y, conj);
            fimcl_free(Y);

            fimcl_t * fX = NULL;

            if(X->clu->prefer_inplace)
            {
                fimcl_ifft_inplace(X);
                fX = X;
                X = NULL;
            } else {
                fX = fimcl_ifft(X);
                fimcl_free(X);
            }

            return fX;
        }



        if(mode == CLU_KEEP_2ND)
        {
            here();
            assert(Y->buf_size_nf > 0);
            assert(X->buf_size_nf > 0);
            fimcl_complex_mul_inplace(Y, X, conj);


            fimcl_t * fX = NULL;
            fX = fimcl_ifft(X);
            fimcl_sync(fX);
            assert(fX->buf_size_nf > 0);
            fimcl_free(X);

            return fX;
        }
    }

    fprintf(stderr, "Can't convolve the objects\n");
    exit(EXIT_FAILURE);

    return NULL;
    assert(NULL);
    return NULL;
}

fimcl_t * fimcl_convolve(fimcl_t * X, fimcl_t * Y, int mode)
{
    return _fimcl_convolve(X, Y, mode, 0);
}
fimcl_t * fimcl_convolve_conj(fimcl_t * X, fimcl_t * Y, int mode)
{
    return _fimcl_convolve(X, Y, mode, 1);
}

void fimcl_sync(fimcl_t * X)
{
    if(X->wait_ev == NULL)
    {
        return;
    }
    check_CL( clu_wait_for_event(X->wait_ev, 10000));
    clReleaseEvent(X->wait_ev);
    X->wait_ev = NULL;

    clFinish(X->clu->command_queue);
    return;
}

void fimcl_preprocess(fimcl_t * fft_image, fimcl_t * fft_PSF, float value)
{
    assert(fft_image->type == fimcl_hermitian);
    assert(fft_PSF->type == fimcl_hermitian);

    printf("fimcl_preprocess, value = %f\n", value);

    assert(fft_PSF->M == fft_image->M);
    cl_kernel kernel = fft_image->clu->kern_preprocess_image.kernel;

    /* Set sizes */
    size_t localWorkSize; /* Use largest possible */
    check_CL ( clGetKernelWorkGroupInfo(kernel,
                                        fft_image->clu->device_id,
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &localWorkSize,
                                        NULL) );

    size_t nel = fimcl_ncx(fft_image); // number of complex values
    size_t numWorkGroups = (nel + (localWorkSize -1) ) / localWorkSize;
    size_t globalWorkSize = localWorkSize * numWorkGroups;


    check_CL( clSetKernelArg(kernel,
                             0, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &fft_image->buf) ); // argument value

    check_CL( clSetKernelArg(kernel,
                             1, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &fft_PSF->buf) ); // argument value

    check_CL( clEnqueueWriteBuffer( fft_image->clu->command_queue,
                                    fft_image->clu->float_gpu,
                                    CL_TRUE,
                                    0,
                                    sizeof(float),
                                    &value,
                                    0, NULL, NULL));

    check_CL( clSetKernelArg(kernel,
                             2,
                             sizeof(cl_mem),
                             (void *) &fft_image->clu->float_gpu) );

    fimcl_sync(fft_image);
    fimcl_sync(fft_PSF);;

    check_CL( clEnqueueNDRangeKernel(fft_image->clu->command_queue,
                                     kernel,
                                     1, //3,
                                     NULL,
                                     &globalWorkSize, //global_work_size,
                                     &localWorkSize,
                                     0,
                                     NULL,
                                     &fft_image->wait_ev) );
    fimcl_sync(fft_image);
}


void fimcl_complex_mul(fimcl_t * X, fimcl_t * Y, fimcl_t * Z, int conj)
{
    assert(X->type == fimcl_hermitian);
    assert(Y->type == fimcl_hermitian);

    if(X->clu->verbose > 1)
    {
        printf("fimcl_complex_mul, conj=%d\n", conj);
    }

    cl_kernel kernel = X->clu->kern_mul.kernel;
    if(conj == 1)
    {
        kernel = X->clu->kern_mul_conj.kernel;
    }


    /* Set sizes */
    size_t localWorkSize; /* Use largest possible */
    check_CL ( clGetKernelWorkGroupInfo(kernel,
                                        X->clu->device_id,
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &localWorkSize,
                                        NULL) );

    size_t nel = fimcl_ncx(X); // number of complex values
    size_t numWorkGroups = (nel + (localWorkSize -1) ) / localWorkSize;
    size_t globalWorkSize = localWorkSize * numWorkGroups;


    check_CL( clSetKernelArg(kernel,
                             0, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &X->buf) ); // argument value

    check_CL( clSetKernelArg(kernel,
                             1, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &Y->buf) ); // argument value

    check_CL( clSetKernelArg(kernel,
                             2, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &Z->buf) ); // argument value


    fimcl_sync(X);
    fimcl_sync(Y);
    fimcl_sync(Z);
    check_CL( clEnqueueNDRangeKernel(X->clu->command_queue,
                                     kernel,
                                     1, //3,
                                     NULL,
                                     &globalWorkSize, //global_work_size,
                                     &localWorkSize,
                                     0,
                                     NULL,
                                     &Z->wait_ev) );
    fimcl_sync(Z);
    return;
}

void fimcl_complex_mul_old(fimcl_t * X, fimcl_t * Y, fimcl_t * Z, int conj)
{
    size_t global_work_offset[] = {0, 0, 0};
    size_t local_work_size[] = {1,1,1};

    cl_kernel kernel = X->clu->kern_mul.kernel;
    if(conj == 1)
    {
        kernel = X->clu->kern_mul_conj.kernel;
    }

    size_t cMNP = fimcl_ncx(X);

    check_CL( clSetKernelArg(kernel,
                             0, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &X->buf) ); // argument value

    check_CL( clSetKernelArg(kernel,
                             1, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &Y->buf) ); // argument value

    check_CL( clSetKernelArg(kernel,
                             2, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &Z->buf) ); // argument value

    fimcl_sync(X);
    fimcl_sync(Y);
    fimcl_sync(Z);
    check_CL( clEnqueueNDRangeKernel(X->clu->command_queue,
                                     kernel,
                                     1, //3,
                                     global_work_offset,
                                     &cMNP, //global_work_size,
                                     local_work_size,
                                     0,
                                     NULL,
                                     &Z->wait_ev) );
    fimcl_sync(Z);
    return;
}


/* Y = X.*Y */
void fimcl_real_mul_inplace(fimcl_t * X, fimcl_t * Y)
{
    assert(X->type == fimcl_real && Y->type == fimcl_real);

    assert(X->M == Y->M);
    assert(X->N == Y->N);
    assert(X->P == Y->P);


    cl_kernel kernel = X->clu->kern_real_mul_inplace.kernel;

    /* Set sizes */
    size_t localWorkSize; /* Use largest possible */
    check_CL ( clGetKernelWorkGroupInfo(kernel,
                                        X->clu->device_id,
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &localWorkSize,
                                        NULL) );

    size_t nel = padsize(X->M) * X->N * X->P;
    size_t numWorkGroups = (nel + (localWorkSize -1) ) / localWorkSize;
    size_t globalWorkSize = localWorkSize * numWorkGroups;

    check_CL( clSetKernelArg(kernel,
                             0, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &X->buf) ); // argument value

    check_CL( clSetKernelArg(kernel,
                             1, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &Y->buf) ); // argument value

    check_CL( clEnqueueNDRangeKernel(X->clu->command_queue,
                                     kernel,
                                     1, //3,
                                     NULL,
                                     &globalWorkSize, //global_work_size,
                                     &localWorkSize,
                                     0,
                                     NULL,
                                     &Y->wait_ev) );

    fimcl_sync(Y);
    return;
}

/* P = x + \alpha*((x-xp))
 * It actually does not matter if the arrays are padded along the first dimension
 * as long as we supply the extended size there is just a little more work to do.
 */
void fimcl_shb_update(fimcl_t * P, fimcl_t * X, fimcl_t * XP, float alpha)
{
    if( (P->type == fimcl_real)
        + (X->type == fimcl_real)
        + (XP->type == fimcl_real) < 3)
    {
        fprintf(stderr,
                "\n"
                "ERROR: fimcl_shb_update only support non-padded data, "
                "i.e. only fimcl_real\n");
        //       exit(EXIT_FAILURE);
    }

    cl_kernel kernel = X->clu->kern_shb_update.kernel;

    /* Set sizes */
    size_t localWorkSize; /* Use largest possible */
    check_CL ( clGetKernelWorkGroupInfo(kernel,
                                        X->clu->device_id,
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &localWorkSize,
                                        NULL) );

    size_t nel = padsize(X->M) * X->N * X->P;
    size_t numWorkGroups = (nel + (localWorkSize -1) ) / localWorkSize;
    size_t globalWorkSize = localWorkSize * numWorkGroups;

    check_CL( clSetKernelArg(kernel,
                             0, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &P->buf) ); // argument value

    check_CL( clSetKernelArg(kernel,
                             1, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &X->buf) ); // argument value

    check_CL( clSetKernelArg(kernel,
                             2, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &XP->buf) ); // argument value

    check_CL( clEnqueueWriteBuffer( X->clu->command_queue,
                                    X->clu->float_gpu,
                                    CL_TRUE,
                                    0,
                                    sizeof(float),
                                    &alpha,
                                    0, NULL, NULL));

    check_CL( clSetKernelArg(kernel,
                             3,
                             sizeof(cl_mem),
                             (void *) &X->clu->float_gpu) );


    check_CL( clEnqueueNDRangeKernel(X->clu->command_queue,
                                     kernel,
                                     1, //3,
                                     NULL,
                                     &globalWorkSize, //global_work_size,
                                     &localWorkSize,
                                     0,
                                     NULL,
                                     &P->wait_ev) );

    fimcl_sync(P);
    return;
}


void fimcl_positivity(fimcl_t * X, float val)
{
    assert(X != NULL);
    cl_kernel kernel = X->clu->kern_real_positivity.kernel;

    /* Set sizes */
    size_t localWorkSize; /* Use largest possible */
    check_CL ( clGetKernelWorkGroupInfo(kernel,
                                        X->clu->device_id,
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &localWorkSize,
                                        NULL) );

    size_t nel = padsize(X->M)*X->N*X->P;
    size_t numWorkGroups = (nel + (localWorkSize -1) ) / localWorkSize;
    size_t globalWorkSize = localWorkSize * numWorkGroups;

    check_CL( clSetKernelArg(kernel,
                             0, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &X->buf) ); // argument value


    check_CL( clEnqueueWriteBuffer( X->clu->command_queue,
                                    X->clu->float_gpu,
                                    CL_TRUE,
                                    0,
                                    sizeof(float),
                                    &val,
                                    0, NULL, NULL));

    check_CL( clSetKernelArg(kernel,
                             1,
                             sizeof(cl_mem),
                             (void *) &X->clu->float_gpu) );

    check_CL( clEnqueueNDRangeKernel(X->clu->command_queue,
                                     kernel,
                                     1, //3,
                                     NULL,
                                     &globalWorkSize, //global_work_size,
                                     &localWorkSize,
                                     0,
                                     NULL,
                                     &X->wait_ev) );

    fimcl_sync(X);
    here();
    return;
}

void fimcl_complex_mul_inplace(fimcl_t * X, fimcl_t * Y, int conj)
{
    assert(X->type == fimcl_hermitian);
    assert(Y->type == fimcl_hermitian);

    if(X->clu->verbose > 1)
    {
        printf("fimcl_complex_mul_inplace, conj=%d\n", conj);
    }

    cl_kernel kernel = X->clu->kern_mul_inplace.kernel;
    if(conj == 1)
    {
        kernel = X->clu->kern_mul_conj_inplace.kernel;
    }

    /* Set sizes */
    size_t localWorkSize; /* Use largest possible */
    check_CL ( clGetKernelWorkGroupInfo(kernel,
                                        X->clu->device_id,
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &localWorkSize,
                                        NULL) );

    size_t nel = fimcl_ncx(X);
    size_t numWorkGroups = (nel + (localWorkSize -1) ) / localWorkSize;
    size_t globalWorkSize = localWorkSize * numWorkGroups;

    fimcl_sync(X);
    fimcl_sync(Y);


    check_CL( clSetKernelArg(kernel,
                             0, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &X->buf) ); // argument value

    check_CL( clSetKernelArg(kernel,
                             1, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &Y->buf) ); // argument value


    check_CL( clEnqueueNDRangeKernel(X->clu->command_queue,
                                     kernel,
                                     1, //3,
                                     NULL,
                                     &globalWorkSize, //global_work_size,
                                     &localWorkSize,
                                     0,
                                     NULL,
                                     &Y->wait_ev) );

    fimcl_sync(Y);
    return;
}


fimcl_t * fimcl_ifft(fimcl_t * fX)
{

#ifdef VKFFT
    fimcl_t * X = fimcl_copy(fX);
    fimcl_ifft_inplace(X);
    return X;
#else

    if(fX->clu->verbose > 1)
    {
        printf("fimcl_ifft\n");
    }
    assert(fX->type == fimcl_hermitian);

    fimcl_t * X = calloc(1, sizeof(fimcl_t));
    X->M = fX->M;
    X->N = fX->N;
    X->P = fX->P;
    X->type = fimcl_real;
    X->wait_ev = CL_SUCCESS;
    X->clu = fX->clu;

    cl_int ret;
    X->buf = clCreateBuffer(fX->clu->context,
                            CL_MEM_READ_WRITE,
                            padsize(fX->M)*fX->N*fX->P *  sizeof(float),
                            NULL, &ret );
    X->buf_size_nf = padsize(fX->M)*fX->N*fX->P;
    fX->clu->nb_allocated += fX->M*fX->N*fX->P*sizeof(float);
    fX->clu->n_alloc++;
    check_CL(ret);

    /* And back again */
    clFinish(X->clu->command_queue);


    check_clFFT(clfftEnqueueTransform(fX->clu->h2r_plan, // XXX
                                      CLFFT_BACKWARD,
                                      1, // numQueuesAndEcents
                                      &fX->clu->command_queue, // commQueues
                                      0, // numWaitEvents
                                      NULL, // waitEvents
                                      &X->wait_ev, // cl_event * outEvents
                                      &fX->buf, // Input buffer
                                      &X->buf, // output buffer
                                      fX->clu->clfft_buffer)); // temp buffer

    clFinish(X->clu->command_queue);
    fimcl_sync(X);
    here();
    return X;
#endif
}

void fimcl_ifft_inplace(fimcl_t * X)
{
    assert(X->type == fimcl_hermitian);

    if(X->clu->verbose > 1)
    {
        printf("fimcl_ifft_inplace\n");
    }


    fimcl_sync(X);

#ifdef VKFFT
    VkFFTLaunchParams launchParams = {};

    launchParams.commandQueue = &X->clu->command_queue;
    /* The main buffer is called buffer and it always has to be
       provided, either during the plan creation or when the plan is
       executed. All calculations are performed in this buffer and it
       is always overwritten. */

    launchParams.buffer = &X->buf;
    /* Buffer size does not have to be supplied with OpenCL */
    int inverse = 1;
    VkFFTResult resFFT = VkFFTAppend(&X->clu->vkfft_app,
                                     inverse, &launchParams);
    if(resFFT != 0)
    {
        fprintf(stderr, "VkFFT failed with error %d\n", resFFT);
        exit(EXIT_FAILURE);
    }


#else
    check_clFFT(clfftEnqueueTransform(X->clu->h2r_inplace_plan,
                                      CLFFT_BACKWARD,
                                      1, // numQueuesAndEvents
                                      &X->clu->command_queue, // commQueues
                                      0, // numWaitEvents
                                      NULL, // waitEvents
                                      &X->wait_ev, // cl_event * outEvents
                                      &X->buf, // Input buffer
                                      NULL, // output buffer
                                      NULL)); // temp buffer
#endif
    clFinish(X->clu->command_queue);

    fimcl_sync(X);
    X->type = fimcl_real;
    return;
}


fimcl_t * fimcl_fft(fimcl_t * X)
{

#ifdef VKFFT
    fimcl_t * Y = fimcl_copy(X);
    fimcl_fft_inplace(Y);
    return Y;
#else

    if(X->clu->verbose > 1)
    {
        printf("fimcl_fft\n");
    }
    if(X->type != fimcl_real)
    {
        fprintf(stderr,
                "\n"
                "ERROR: fimcl_fft was not called with a fimcl_real object\n");
        exit(EXIT_FAILURE);
    }
    fimcl_sync(X);

    /* Number of complex elements in the FFT */
    size_t cMNP = fimcl_ncx(X);

    fimcl_t * fX = calloc(1, sizeof(fimcl_t));
    fX->M = X->M;
    fX->N = X->N;
    fX->P = X->P;
    fX->type = fimcl_hermitian;

    fX->clu = X->clu;

    cl_int ret;
    fX->buf_size_nf = cMNP * 2;
    fX->buf = clCreateBuffer(X->clu->context,
                             CL_MEM_READ_WRITE,
                             fX->buf_size_nf *  sizeof(float),
                             NULL, &ret );
    X->clu->nb_allocated += fX->buf_size_nf *  sizeof(float);
    X->clu->n_alloc++;
    check_CL(ret);




    check_clFFT(clfftEnqueueTransform(X->clu->r2h_plan,
                                      CLFFT_FORWARD, // direction
                                      1, // numQueuesAndEcents
                                      &X->clu->command_queue, // commQueues
                                      0, // numWaitEvents
                                      NULL, // waitEvents
                                      &fX->wait_ev, // cl_event * outEvents
                                      &X->buf, // Input buffer
                                      &fX->buf, // output buffer
                                      X->clu->clfft_buffer)); // temp buffer

    fimcl_sync(fX);
    return fX;
#endif
}

void fimcl_fft_inplace(fimcl_t * X)
{
    assert(X != NULL);
    assert(X->type == fimcl_real);

    if(X->clu->verbose > 1)
    {
        printf("fimcl_fft_inplace\n");
    }

    /* Number of complex elements in the FFT */
    size_t cMNP = fimcl_ncx(X);

    if(X->buf_size_nf < cMNP * 2)
    {
        // Set in fimcl_new and fimcl_copy
        fprintf(stderr,
                "fimcl_fft_inplace: the object does not have a "
                "large enough buffer for in-place fft\n"
                "current size: %zu required size: %zu\n",
                X->buf_size_nf, cMNP*2);

        exit(EXIT_FAILURE);
    }

#ifdef VKFFT
    VkFFTLaunchParams launchParams = {};

    launchParams.commandQueue = &X->clu->command_queue;
    /* The main buffer is called buffer and it always has to be
       provided, either during the plan creation or when the plan is
       executed. All calculations are performed in this buffer and it
       is always overwritten. */

    launchParams.buffer = &X->buf;
    /* Buffer size does not have to be supplied with OpenCL */
    int inverse = 0;
    if(X->clu->verbose > 1)
    {
        printf("VkFFTAppend (for in-place forward transform)\n");
    }
    VkFFTResult resFFT = VkFFTAppend(&X->clu->vkfft_app,
                                     inverse, &launchParams);
    if(resFFT != 0)
    {
        fprintf(stderr, "VkFFT failed with error %d\n", resFFT);
        exit(EXIT_FAILURE);
    }

    clFinish(X->clu->command_queue);

#else
    fimcl_sync(X);
    check_clFFT(clfftEnqueueTransform(X->clu->r2h_inplace_plan,
                                      CLFFT_FORWARD, // direction
                                      1, // numQueuesAndEcents
                                      &X->clu->command_queue, // commQueues
                                      0, // numWaitEvents
                                      NULL, // waitEvents
                                      &X->wait_ev, // cl_event * outEvents
                                      &X->buf, // Input buffer
                                      NULL, // output buffer
                                      X->clu->clfft_buffer)); // temp buffer
    fimcl_sync(X);
#endif

    X->type = fimcl_hermitian;
    return;
}


float * clu_convolve(clu_env_t * clu,
                     int dropy,
                     float * X, float * Y,
                     size_t M, size_t N, size_t P)
{
    size_t cMNP = M_r2h(M)*N*P;
    if(clu->verbose > 1)
    {
        size_t buf_mem = 2*cMNP*2*sizeof(float)/1000000;
        printf("clu_convolve: GPU memory needed for storage: %zu MB\n"
               "              Approximate total need: %zu\n",
               buf_mem, 15*buf_mem/10);
    }
    // Note full size specified for in-place transformations
    struct timespec t0, t1, t2;
    dw_gettime(&t0);
    fimcl_t * gX = fimcl_new(clu, fimcl_real, X, M, N, P);
    dw_gettime(&t1);
    fimcl_t * gY = fimcl_new(clu, fimcl_real, Y, M, N, P);
    dw_gettime(&t2);

    float  dt0 = clockdiff(&t0, &t1);
    float  dt1 = clockdiff(&t1, &t2);
    if(clu->verbose > 2)
    {
        printf("## Buffer upload times: %f + %f = %f\n", dt0, dt1, dt0+dt1);
    }

    /* Convolve and free gX and gY as soon as possible */
    fimcl_t * gZ = fimcl_convolve(gX, gY, CLU_DROP_ALL);

    if(dropy)
    {
        fim_free(Y);
    }
    float * Z = fimcl_download(gZ);
    fimcl_free(gZ);
    return Z;
}

static double clockdiff(struct timespec* start, struct timespec * finish)
{
    double elapsed = (finish->tv_sec - start->tv_sec);
    elapsed += (finish->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}

char * get_cl_device_info_string(cl_device_id dev_id, cl_device_info info_type)
{
    size_t string_size = 16*1024;
    size_t string_size_used = 0;
    char * string = calloc(string_size, 1);
    // Note: Can return CL_INVALID_VALUE if the string_size is not large enough.
    check_CL(clGetDeviceInfo(dev_id,
                             info_type,
                             string_size,
                             string,
                             &string_size_used));

    return string;

}

void clu_print_device_info(FILE * fid, cl_device_id dev_id)
{
    cl_device_type dev_type;
    check_CL( clGetDeviceInfo(dev_id, // 	cl_device_id device,
                              CL_DEVICE_TYPE,
                              sizeof(cl_device_type),
                              &dev_type,
                              NULL) );

    fprintf(fid, "CL_DEVICE_TYPE=");
    switch(dev_type)
    {
    case CL_DEVICE_TYPE_CPU:
        fprintf(fid, "CL_DEVICE_TYPE_CPU\n");
        break;
    case CL_DEVICE_TYPE_GPU:
        fprintf(fid, "CL_DEVICE_TYPE_GPU\n");
        break;
    case CL_DEVICE_TYPE_ACCELERATOR:
        fprintf(fid, "CL_DEVICE_TYPE_ACCELERATOR\n");
        break;
    case CL_DEVICE_TYPE_DEFAULT:
        fprintf(fid, "CL_DEVICE_TYPE_DEFAULT\n");
        break;
    }

    cl_ulong dev_memory;
    check_CL( clGetDeviceInfo(dev_id, // 	cl_device_id device,
                              CL_DEVICE_GLOBAL_MEM_SIZE, // cl_device_info param_name,
                              sizeof(cl_ulong), //size_t param_value_size,
                              &dev_memory,
                              NULL));

    fprintf(fid, "CL_DEVICE_GLOBAL_MEM_SIZE = %" PRIu64 " (%" PRIu64 " MiB)\n",
            dev_memory,
            dev_memory/1000000);

    char * string = NULL;
    string = get_cl_device_info_string(dev_id, CL_DEVICE_NAME);
    fprintf(fid, "CL_DEVICE_NAME = %s\n", string);
    free(string);
    string = get_cl_device_info_string(dev_id, CL_DEVICE_VENDOR);
    fprintf(fid, "CL_DEVICE_VENDOR = %s\n", string);
    free(string);
    string = get_cl_device_info_string(dev_id, CL_DRIVER_VERSION);
    fprintf(fid, "CL_DRIVER_VERSION = %s\n", string);
    free(string);
    string = get_cl_device_info_string(dev_id, CL_DEVICE_EXTENSIONS);
    fprintf(fid, "CL_DEVICE_EXTENSIONS = %s\n", string);
    free(string);


    return;
}

static char * read_program(const char * fname, size_t * size)
{
    FILE * fid = fopen(fname, "r");
    if (!fid) {
        fprintf(stderr, "Failed to load kernel from %s.\n", fname);
        exit(EXIT_FAILURE);
    }

    fseek(fid, 0, SEEK_END);
    long fsize = ftell(fid);
    fseek(fid, 0, SEEK_SET);

    char * source_str = (char*) malloc(fsize); // XXX
    size[0] = fread( source_str, 1, fsize, fid);
    fclose(fid);
    return source_str;
}

/** Callback from clCreateContext */
void cl_crash_fun (
                   const char *errinfo,
                   __attribute__((unused)) const void *private_info,
                   __attribute__((unused)) size_t cb,
                   __attribute__((unused)) void *user_data
                   )
{
    printf("On no... bad new from OpenCL:\n");
    printf("errinfo: %s\n", errinfo);
    printf("\n");
    //printf("private_info: %s\n", private_info);
    fflush(stdout);
    return;
}

clu_env_t * clu_new(int verbose, int cl_device)
{

    clu_env_t * env = calloc(1, sizeof(clu_env_t));
    env->verbose = verbose;


#ifndef VKFFT
    /* Might be of use in future releases */
    setenv("CLFFT_REQUEST_LIB_NOMEMALLOC", "1", 1);

    assert(2*M_r2h(128) == padsize(128));
    assert(2*M_r2h(129) == padsize(129));

    char * clFFT_cache_path = malloc(1024);
    char * dir_home = getenv("HOME");
    sprintf(clFFT_cache_path, "%s/.config/deconwolf/clFFT/", dir_home);
    setenv("CLFFT_CACHE_PATH", clFFT_cache_path, 1);
    ensuredir(clFFT_cache_path);
    free(clFFT_cache_path);
    env->clfft_buffer = NULL;
    env->clfft_buffer_size = 0;
    env->clFFT_loaded = 0;
#endif



    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;

    check_CL( clGetPlatformIDs(1, &env->platform_id, &ret_num_platforms) );

    if(env->verbose > 1)
    {
        printf("Found %d CL platforms\n", ret_num_platforms);
    }
    if(ret_num_platforms < 1)
    {
        fprintf(stderr, "Unable to find a CL platform\n");
        return NULL;
    }

    /* Get the number of devices */
    check_CL( clGetDeviceIDs(env->platform_id,
                             CL_DEVICE_TYPE_ALL,
                             0,
                             NULL,
                             &ret_num_devices));
    if(ret_num_devices < 1)
    {
        fprintf(stderr, "Could not find any OpenCL devices\n");
        return NULL;
    }
    cl_device_id * devices = malloc(ret_num_devices*sizeof(cl_device_id));
    /* Get the devices */
    check_CL( clGetDeviceIDs(env->platform_id,
                             CL_DEVICE_TYPE_ALL,
                             ret_num_devices,
                             devices,
                             NULL));


    if(env->verbose > 1)
    {
        printf("Found %d CL devices\n", ret_num_devices);
        printf("Will use device %d (first = 0)\n", cl_device);
    }
    fflush(stdout);

    env->device_id = devices[cl_device];

    if(env->verbose > 1)
    {
        for(cl_uint kk = 0; kk<ret_num_devices; kk++)
	{
            printf("CL device #%d\n", kk);
            clu_print_device_info(stdout, devices[kk]);
	}
    }

    cl_int ret = CL_SUCCESS;
    env->context = clCreateContext( NULL,
                                    1,
                                    &env->device_id,
                                    cl_crash_fun,
                                    NULL,
                                    &ret);
    check_CL(ret);

    // Create a command queue
    env->command_queue = clCreateCommandQueue(env->context,
                                              env->device_id,
                                              0,
                                              &ret);
    check_CL(ret);
    if(env->verbose > 1)
    {
        printf("CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS=%u\n",
               CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS);
    }

    if(env->verbose > 1)
    {
#ifdef VKFFT
        printf("Using VkFFT version %d\n", VkFFTGetVersion());
#else
        cl_uint major, minor, patch;
        clfftGetVersion(&major, &minor, &patch);
        printf("clfft version: %d.%d.%d\n", major, minor, patch);
#endif
    }

    env->prefer_inplace = 1;
    return env;
}

void clu_prepare_kernels(clu_env_t * clu,
                         size_t M, size_t N, size_t P,
                         size_t wM, size_t wN, size_t wP)
{
    if(clu->verbose > 1)
    {
        printf("Preparing for convolutions of size %zu x %zu x %zu\n",
               wM, wN, wP);
    }

#ifndef VKFFT
    assert(clu->clFFT_loaded == 0);

    /* Setup clFFT. */
    check_CL( clfftInitSetupData(&clu->fftSetup));
    check_CL( clfftSetup(&clu->fftSetup));

    if(clu->verbose > 1)
    {
        cl_uint major, minor, patch;
        check_clFFT(clfftGetVersion(&major, &minor, &patch));
        printf("clFFT version %d.%d.%d\n", major, minor, patch);
    }

    /* Create clFFT plans for this particular size */
    clu->r2h_plan = gen_r2h_plan(clu, wM, wN, wP);
    clu->r2h_inplace_plan = gen_r2h_inplace_plan(clu, wM, wN, wP);

    clu->h2r_plan = gen_h2r_plan(clu, wM, wN, wP);
    clu->h2r_inplace_plan = gen_h2r_inplace_plan(clu, wM, wN, wP);
    clu->clFFT_loaded = 1;
#endif

    char * argstring = malloc(1024);

    /*
     *  Create the OpenCL kernels that we will use
     */

    /* Kernels for complex floats
     * These kernels only know the number of complex values
     */
    sprintf(argstring, "-D cMNP=%zu",
            M_r2h(wM)*wN*wP);
    clu->kern_mul =
        *clu_kernel_newa(clu,
                         NULL, //"cl_complex_mul.c",
                         (const char *) cl_complex_mul,
                         cl_complex_mul_len,
                         "cl_complex_mul",
                         argstring);

    clu->kern_mul_conj =
        *clu_kernel_newa(clu,
                         NULL, //"cl_complex_mul_conj.c",
                         (const char *) cl_complex_mul_conj,
                         cl_complex_mul_conj_len,
                         "cl_complex_mul_conj",
                         argstring);
    clu->kern_mul_inplace =
        *clu_kernel_newa(clu,
                         NULL,
                         (const char *) cl_complex_mul_inplace,
                         cl_complex_mul_inplace_len,
                         "cl_complex_mul_inplace",
                         argstring);
    clu->kern_mul_conj_inplace =
        *clu_kernel_newa(clu,
                         NULL,
                         (const char *) cl_complex_mul_conj_inplace,
                         cl_complex_mul_conj_inplace_len,
                         "cl_complex_mul_conj_inplace",
                         argstring);
    clu->kern_preprocess_image =
        * clu_kernel_newa(clu,
                          NULL,
                          (const char *) src_kernels_cl_preprocess_image_c,
                          src_kernels_cl_preprocess_image_c_len,
                          "cl_preprocess_image",
                          argstring);

    /* Kernels for real float */
    sprintf(argstring, "-D wMNP=%zu", padsize(wM)*wN*wP);
    clu->kern_real_mul_inplace =
        * clu_kernel_newa(clu,
                          NULL,
                          (const char *) src_kernels_cl_real_mul_inplace_c,
                          src_kernels_cl_real_mul_inplace_c_len,
                          "cl_real_mul_inplace",
                          argstring);
    clu->kern_real_positivity =
        * clu_kernel_newa(clu,
                          NULL,
                          (const char *) src_kernels_cl_positivity_c,
                          src_kernels_cl_positivity_c_len,
                          "cl_positivity",
                          argstring);
    clu->kern_shb_update =
        * clu_kernel_newa(clu,
                          NULL,
                          (const char *) src_kernels_cl_shb_update_c,
                          src_kernels_cl_shb_update_c_len,
                          "cl_shb_update",
                          argstring);

    /* CL_DEVICE_MAX_CONSTANT_ARGS is at least 8 so 7 is safe
       without checking capabilities */

    /* Number of voxels to process per work item
       Please check so that it matches when the kernel is called */
    int nPerItem = 16;

    sprintf(argstring,
            "-DNELEMENTS=%zu "
            "-DM0=%zu "
            "-DM=%zu -DN=%zu -DP=%zu "
            "-DwM=%zu -DwN=%zu -DwP=%zu "
            "-DnPerItem=%d",
            padsize(wM)*wN*wP,
            M,
            padsize(M), N, P,
            padsize(wM), wN, wP,
            nPerItem);

    clu->idiv_kernel = clu_kernel_newa(clu,
                                       "kernels/cl_idiv_kernel.c",
                                       (const char *) cl_idiv_kernel,
                                       cl_idiv_kernel_len,
                                       "idiv_kernel",
                                       argstring);

    sprintf(argstring,
            "-DNELEMENTS=%zu "
            "-DM0=%zu "
            "-DM=%zu -DN=%zu -DP=%zu "
            "-DwM=%zu -DwN=%zu -DwP=%zu",
            padsize(wM)*wN*wP,
            M,
            padsize(M), N, P,
            padsize(wM), wN, wP);

    clu->update_y_kernel = clu_kernel_newa(clu,
                                           "kernels/cl_update_y_kernel.c",
                                           (const char *) cl_update_y_kernel,
                                           cl_update_y_kernel_len,
                                           "update_y_kernel",
                                           argstring);

    free(argstring);

    float value = 0;
    cl_int status;
    clu->float_gpu = clCreateBuffer(clu->context,
                                    CL_MEM_READ_ONLY,
                                    sizeof(float),
                                    NULL, &status);
    check_CL(status);
    clu->n_alloc++;

    check_CL( clEnqueueWriteBuffer( clu->command_queue,
                                    clu->float_gpu,
                                    CL_TRUE,
                                    0,
                                    sizeof(float),
                                    &value,
                                    0, NULL, NULL));



#ifdef VKFFT

    char * vkfft_cache_file_name = NULL;
    char * dir_home = getenv("HOME");
    if(dir_home == NULL)
    {
        if(clu->verbose > 0)
        {
            printf("Warning: Will write the VkFFT configuration in the current folder.\n");
            if(clu->verbose > 1) {
                printf("Reason: Can not determine a suitable folder under Windows.\n");
            }
        }
        vkfft_cache_file_name = malloc(1024);
        sprintf(vkfft_cache_file_name, "VkFFT_kernelCache_%zux%zux%zu.binary",
                wM, wN, wP);
    } else {
        vkfft_cache_file_name = malloc(strlen(dir_home) + 1024);
        sprintf(vkfft_cache_file_name, "%s/.config/deconwolf/VkFFT_kernelCache_%zux%zux%zu.binary",
                dir_home,
                wM, wN, wP);
    }

    if(clu->verbose > 1)
    {
        printf("vkFFT cache file: %s\n", vkfft_cache_file_name);
    }

    int load_vkfft_from_file = 0;
    if(dw_isfile(vkfft_cache_file_name))
    {
        load_vkfft_from_file = 1;
    }

    struct timespec t0, t1;


    dw_gettime(&t0);
    VkFFTConfiguration configuration = {};

    //Device management + code submission
    configuration.device = &clu->device_id;
    configuration.platform = &clu->platform_id;
    configuration.context = &clu->context;
    configuration.performR2C = 1; /*  Perform R2C/C2R decomposition */
    configuration.FFTdim = 3; //FFT dimension
    configuration.size[0] = wM; //FFT size -- number of complex numbers in this case
    configuration.size[1] = wN;
    configuration.size[2] = wP;

    //configuration.isOutputFormatted = 1; // padded for output size
    configuration.normalize = 1; // Normalize inverse transform (0 - off, 1- on
    // printf("sizeof(VkFFTApplication))=%zu\n", sizeof(VkFFTApplication));;
    if( load_vkfft_from_file )
    {
        configuration.loadApplicationFromString = 1;
    } else {
        configuration.saveApplicationToString = 1;
    }

    if(clu->verbose > 1)
    {
        printf("Initializing VkFFT for size %" PRIu64 " x %" PRIu64 " x %" PRIu64 "\n",
               configuration.size[0],
               configuration.size[1],
               configuration.size[2]);
    }
    if(configuration.loadApplicationFromString)
    {
        FILE* kernelCache;
        uint64_t str_len;
        kernelCache = fopen(vkfft_cache_file_name, "rb");
        if(kernelCache == NULL)
        {
            fprintf(stderr, "Unable to load %s\n", vkfft_cache_file_name);
            configuration.loadApplicationString = 0;
            goto skipLoadApplicationFromString;
        }
        fseek(kernelCache, 0, SEEK_END);
        str_len = ftell(kernelCache);
        fseek(kernelCache, 0, SEEK_SET);
        configuration.loadApplicationString = malloc(str_len);
        size_t nread = fread(configuration.loadApplicationString,
                             1, str_len, kernelCache);
        if(nread != str_len)
        {
            fprintf(stderr, "Unable to read from %s\n", vkfft_cache_file_name);
            fprintf(stderr, "The file might be corrupt or generated by an "
                    "incompatible version. Try delete it and run again\n");
            fprintf(stderr, "nread=%zu str_len=%zu\n", nread, str_len);
            exit(EXIT_FAILURE);
        }
        fclose(kernelCache);
    skipLoadApplicationFromString:
        ;
    }
    if(clu->verbose > 1)
    {
        printf("VkFFT configuration:\n");
        printf("loadApplicationFromString: %d\n", configuration.loadApplicationFromString);
        printf("saveApplicationToString: %d\n", configuration.saveApplicationToString);
    }

    VkFFTResult resFFT = initializeVkFFT(&clu->vkfft_app, configuration);
    if(resFFT != 0)
    {
        fprintf(stderr, "\n"
                "initializeVkFFT failed with error code %d\n", resFFT);
        fprintf(stderr,
                "%s\n", getVkFFTErrorString(resFFT));
        fprintf(stderr,
                "\n"
                "Unfortunately deconwolf can not continue from this point.\n"
                "Please see the documentation about possible solutions or open\n"
                "an issue on https://github.com/elgw/deconwolf/issues if you\n"
                "think that this is a bug.\n");
        fprintf(stderr,
                "\n"
                "In most cases errors at this point indicate that there\n"
                "was not enough memory on the GPU\n");
        exit(EXIT_FAILURE);
    }
    if(clu->verbose > 1)
    {
        printf("VkFFT was initialized\n");
    }
    if (configuration.loadApplicationFromString) {
        free(configuration.loadApplicationString);
    }
    if (configuration.saveApplicationToString) {
        FILE* kernelCache;
        kernelCache = fopen(vkfft_cache_file_name, "wb");
        if(kernelCache == NULL)
        {
            fprintf(stderr, "Unable to open %s for writing\n", vkfft_cache_file_name);
        } else {
            fwrite(clu->vkfft_app.saveApplicationString,
                   clu->vkfft_app.applicationStringSize, 1, kernelCache);
            fclose(kernelCache);
        }
    }
    dw_gettime(&t1);
    if(clu->verbose > 2)
    {
        // 0.6 s for 512x256x128
        printf("VkFFT intialization took %f s\n", timespec_diff(&t1, &t0));
    }

    free(vkfft_cache_file_name);
#endif

    return;
}

void clu_destroy(clu_env_t * clu)
{
    if(clu->verbose > 1)
    {
        printf("Closing the OpenCL environment\n");
    }

#ifdef VKFFT
    deleteVkFFT(&clu->vkfft_app);
#else
    /* Release clFFT library. */
    if(clu->clFFT_loaded)
    {
        /* Teardown */
        check_clFFT(clfftDestroyPlan( &clu->r2h_plan ));
        check_clFFT(clfftDestroyPlan( &clu->h2r_plan ));
        check_clFFT(clfftDestroyPlan( &clu->r2h_inplace_plan ));
        check_clFFT(clfftDestroyPlan( &clu->h2r_inplace_plan ));
        check_clFFT(clfftTeardown());

        /* Free memory */
        if(clu->clfft_buffer_size > 0)
        {
            clu_release(clu, clu->clfft_buffer);
        }
    }
#endif

    if(clu->float_gpu != NULL)
    {
        clu_release(clu, clu->float_gpu);
    }

    /* Clear up the OpenCL stuff */
    check_CL(clFlush(clu->command_queue));
    /* For complex data */
    clu_kernel_destroy(clu->kern_mul);
    clu_kernel_destroy(clu->kern_mul_conj);
    clu_kernel_destroy(clu->kern_mul_inplace);
    clu_kernel_destroy(clu->kern_mul_conj_inplace);
    /* For real data */

    clu_kernel_destroy(clu->kern_real_mul_inplace);
    clu_kernel_destroy(clu->kern_real_positivity);

    clu_kernel_destroy(clu->kern_shb_update);
    clu_kernel_destroy(clu->kern_preprocess_image);

    clu_kernel_destroy(*clu->idiv_kernel);

    clu_kernel_destroy(*clu->update_y_kernel);

    check_CL(clFinish(clu->command_queue));
    check_CL(clReleaseCommandQueue(clu->command_queue) );
    check_CL(clReleaseContext(clu->context));


    if(clu->n_alloc != clu->n_release)
    {
        fprintf(stderr, "clu_destroy warning: allocated: %zu, released: %zu\n",
                clu->n_alloc, clu->n_release);
        fprintf(stderr, "This indicates a potential memory leak on the GPU.\n");
    }

    free(clu);
    return;
}


clu_kernel_t * clu_kernel_new(clu_env_t * env,
                              const char * file,
                              const char * program_code,
                              size_t program_size,
                              const char * kernel_name)
{
    return clu_kernel_newa(env, file, program_code, program_size, kernel_name, NULL);
}


clu_kernel_t * clu_kernel_newa(clu_env_t * env,
                               const char * file,
                               const char * program_code,
                               size_t program_size,
                               const char * kernel_name,
                               const char * argument_string)
{
    clu_kernel_t * clk = calloc(1, sizeof(clu_kernel_t));
    cl_int ret = CL_SUCCESS;

    // Create a program from the kernel source
    size_t source_size = 0;
    const char * source_str = NULL;
    char * file_str = NULL;

    if(program_code == NULL)
    {
        char * file_str = read_program(file, &source_size);
        source_str = (const char *) file_str;
    } else {
        if(program_code == NULL)
        {
            return NULL;
        }
        source_str = program_code;
        source_size = program_size;
    }

    assert(source_str != NULL);
    //printf("1 clCreateProgram\n");
    clk->program = clCreateProgramWithSource(env->context, 1,
                                             (const char **) &source_str,
                                             (const size_t *) &source_size,
                                             &ret);
    check_CL(ret);
    if(file_str != NULL)
    {
        free(file_str);
    }

    cl_uint err = clBuildProgram(clk->program,
                                 1,
                                 &env->device_id,
                                 argument_string,
                                 NULL,
                                 NULL);

    if (err != CL_SUCCESS) {
        char *buff_erro;
        cl_int errcode;
        size_t build_log_len;
        errcode = clGetProgramBuildInfo(clk->program,
                                        env->device_id,
                                        CL_PROGRAM_BUILD_LOG,
                                        0,
                                        NULL,
                                        &build_log_len);
        if (errcode) {
            printf("clGetProgramBuildInfo failed at line %d\n", __LINE__);
            exit(-1);
        }

        buff_erro = malloc(build_log_len);
        if (!buff_erro) {
            printf("malloc failed at line %d\n", __LINE__);
            exit(-2);
        }

        errcode = clGetProgramBuildInfo(clk->program,
                                        env->device_id,
                                        CL_PROGRAM_BUILD_LOG,
                                        build_log_len,
                                        buff_erro,
                                        NULL);
        check_CL(errcode);

        fprintf(stderr, "\n");
        fprintf(stderr, "Sorry, a fatal error occurred and dw can't continue\n");
        fprintf(stderr, "File: %s, Function: %s Line: %d\n",
                __FILE__, __FUNCTION__, __LINE__);
        fprintf(stderr, "This happened when trying to compile the kernel %s\n",
                kernel_name);
        fprintf(stderr,"Build log: \n%s\n", buff_erro);
        free(buff_erro);

        exit(EXIT_FAILURE);
    }


    // Create the OpenCL kernel
    clk->kernel = clCreateKernel(clk->program, kernel_name, &ret);
    check_CL(ret);
    return clk;
}


void clu_kernel_destroy(clu_kernel_t kern)
{
    check_CL( clReleaseKernel(kern.kernel));
    check_CL( clReleaseProgram(kern.program));
    //printf("1 clReleaseProgram\n");
    return;
}


cl_uint clu_wait_for_event(cl_event clev, size_t ns)
{
    cl_uint ret = CL_SUCCESS;
    size_t param_value_size = sizeof(cl_int);
    size_t param_value_size_ret = 0;
    struct timespec sleep_req = {};
    sleep_req.tv_nsec = ns;
    struct timespec sleep_rem = {};
    cl_uint status = CL_QUEUED;

    /* No sleeping at all if it is already done */
    ret = clGetEventInfo(clev,
                         CL_EVENT_COMMAND_EXECUTION_STATUS,
                         param_value_size,
                         &status,
                         &param_value_size_ret);

    if(status == CL_COMPLETE || ret != CL_SUCCESS)
    {
        return ret;
    }

    /* If it wasn't done, wait for that */
    while(status != CL_COMPLETE)
    {
        // TODO: Bussy waiting on windows ...
#ifndef WINDOWS
        nanosleep(&sleep_req,
                  &sleep_rem);
#endif
        ret = clGetEventInfo(clev,
                             CL_EVENT_COMMAND_EXECUTION_STATUS,
                             param_value_size,
                             &status,
                             &param_value_size_ret);
        if(ret != CL_SUCCESS)
            break;
    }
    return ret;
}

size_t clu_next_fft_size(size_t N)
{
    const int has_factor_11 = 1;
    const int has_factor_13 = 1;

    while(1)
    {
        size_t t = N;
        while( t > 0 && (t % 2) == 0)
        {
            t/=2;
        }
        while( t > 0 && (t % 3) == 0)
        {
            t/=3;
        }
        while( t > 0 && (t % 5) == 0)
        {
            t/=5;
        }
        while( t > 0 && (t % 7) == 0)
        {
            t/=7;
        }
        if(has_factor_11)
        {
            while( t > 0 && (t % 11) == 0)
            {
                t/=11;
            }
        }
        if(has_factor_13)
        {
            while( t > 0 && (t % 13) == 0)
            {
                t/=13;
            }
        }
        if(t == 1)
        {
            return N;
        }
        N++;
    }
    return N;
}


#ifndef VKFFT
const char * get_clfft_error_string(clfftStatus error)
{
    switch (error) {
    case CLFFT_BUGCHECK:
        return "CLFFT_BUGCHECK";
    case CLFFT_NOTIMPLEMENTED:
        return "CLFFT_NOTIMPLEMENTED";
    case CLFFT_TRANSPOSED_NOTIMPLEMENTED:
        return "CLFFT_TRANSPOSED_NOTIMPLEMENTED";
    case CLFFT_FILE_NOT_FOUND:
        return "CLFFT_FILE_NOT_FOUND";
    case CLFFT_FILE_CREATE_FAILURE:
        return "CLFFT_FILE_CREATE_FAILURE";
    case CLFFT_VERSION_MISMATCH:
        return "CLFFT_VERSION_MISMATCH";
    case CLFFT_INVALID_PLAN:
        return "CLFFT_INVALID_PLAN";
    case CLFFT_DEVICE_NO_DOUBLE:
        return "CLFFT_DEVICE_NO_DOUBLE";
    case CLFFT_DEVICE_MISMATCH:
        return "CLFFT_DEVICE_MISMATCH";
    default:
        return "UNKNOWN";
    }
}
#endif

#ifndef VKFFT
clfftPlanHandle  gen_r2h_plan(clu_env_t * clu,
                              size_t M, size_t N, size_t P)
{
    size_t cM = M_r2h(M);
    size_t cN = N;

    clfftPlanHandle planHandle;
    size_t real_size[3] = {M, N, P};


    //printf("Real size: %zu x %zu x %zu\n", M, N, P);

    /* Start with a default plan and then adjust it */
    check_clFFT(clfftCreateDefaultPlan(&planHandle,
                                       clu->context,
                                       CLFFT_3D,
                                       real_size));

    check_clFFT(clfftSetPlanPrecision(planHandle, CLFFT_SINGLE));
    check_clFFT(clfftSetLayout(planHandle,
                               CLFFT_REAL,
                               CLFFT_HERMITIAN_INTERLEAVED));
    check_clFFT(clfftSetResultLocation(planHandle,
                                       CLFFT_OUTOFPLACE));

#if PAD_FIRST_DIM
    size_t inStride[3] =  {1, padsize(M), padsize(M)*N};
#else
    size_t inStride[3] =  {1, M, M*N};
#endif
    size_t outStride[3] = {1, cM, cM*cN};

    check_clFFT(clfftSetPlanInStride(planHandle,  CLFFT_3D, inStride));
    check_clFFT(clfftSetPlanOutStride(planHandle, CLFFT_3D, outStride));

    /* Bake the plan -- i.e. generate OpenCL kernel */
    check_clFFT(clfftBakePlan(planHandle, 1,
                              &clu->command_queue,
                              NULL, NULL));
    check_CL(clFinish(clu->command_queue));

#ifndef VKFFT
    size_t clfft_buffersize = 0;
    check_clFFT(clfftGetTmpBufSize(planHandle, &clfft_buffersize));
    if(clu->verbose > 2)
    {
        if(clfft_buffersize > 0 && clu->verbose > 2)
        {
            printf("Will require a buffer of size: %zu\n", clfft_buffersize);
        }
    }
    clu_increase_clfft_buffer(clu, clfft_buffersize);
#endif

    return planHandle;
}
#endif


#ifndef VKFFT
/* TEST: Will inplace be much faster when the first dimension is padded? */
clfftPlanHandle  gen_r2h_inplace_plan(clu_env_t * clu,
                                      size_t M, size_t N, size_t P)
{
    size_t cM = M_r2h(M);
    size_t cN = N;

    clfftPlanHandle planHandle;
    size_t real_size[3] = {M, N, P};

    /* Start with a default plan and then adjust it */
    check_clFFT(clfftCreateDefaultPlan(&planHandle,
                                       clu->context,
                                       CLFFT_3D,
                                       real_size));

    check_clFFT(clfftSetPlanPrecision(planHandle, CLFFT_SINGLE));
    check_clFFT(clfftSetLayout(planHandle,
                               CLFFT_REAL,
                               CLFFT_HERMITIAN_INTERLEAVED));

    check_clFFT(clfftSetResultLocation(planHandle,
                                       CLFFT_INPLACE));

#if PAD_FIRST_DIM
    size_t inStride[3] = {1, padsize(M), padsize(M)*N};
#else
    printf("gen_r2h_inplace_plan -- bad values\n");
    size_t inStride[3] = {1, M, M*N};
#endif

    //size_t inStride[3] = {1, M, M*N};
    size_t outStride[3] = {1, cM, cM*cN};

    check_clFFT(clfftSetPlanInStride(planHandle,  CLFFT_3D, inStride));
    check_clFFT(clfftSetPlanOutStride(planHandle, CLFFT_3D, outStride));

    /* Bake the plan -- i.e. generate OpenCL kernel */
    check_clFFT(clfftBakePlan(planHandle, 1,
                              &clu->command_queue,
                              NULL, NULL));

    check_CL(clFinish(clu->command_queue));

    size_t clfft_buffersize = 0;
    check_clFFT(clfftGetTmpBufSize(planHandle, &clfft_buffersize));

    if(clu->verbose > 2)
    {
        if(clfft_buffersize > 0)
        {
            printf("gen_r2h_inplace_plan: Will require a buffer of size: %zu\n", clfft_buffersize);
        } else {
            printf("gen_r2h_inplace_plan: No buffer required\n");
        }
    }

    clu_increase_clfft_buffer(clu, clfft_buffersize);

    return planHandle;
}
#endif

#ifndef VKFFT
cl_int clu_increase_clfft_buffer(clu_env_t * clu, size_t req_size)
{
    cl_int ret = CL_SUCCESS;

    if(req_size <= clu->clfft_buffer_size)
    {
        return ret;
    }

    if(clu->clfft_buffer_size > 0)
    {
        /* If there was already a smaller buffer allocated */
        clu_release(clu, clu->clfft_buffer);
    }
    clu->clfft_buffer = clCreateBuffer(clu->context,
                                       CL_MEM_READ_WRITE,
                                       req_size,
                                       NULL, &ret );
    clu->nb_allocated += req_size;
    clu->n_alloc++;
    check_CL(ret);
    clu->clfft_buffer_size = req_size;

    return ret;
}
#endif

#ifndef VKFFT
clfftPlanHandle  gen_h2r_plan(clu_env_t * clu,
                              size_t M, size_t N, size_t P)
{
    size_t cM = M_r2h(M);
    size_t cN = N;

    /* From Hermitian to Real of size MxNxP */
    clfftPlanHandle planHandle;

    size_t real_size[3] = {M, N, P};

    /* Start with a default plan and then adjust it */
    check_clFFT(clfftCreateDefaultPlan(&planHandle,
                                       clu->context,
                                       CLFFT_3D,
                                       real_size));

    check_clFFT(clfftSetPlanPrecision(planHandle, CLFFT_SINGLE));
    check_clFFT(clfftSetLayout(planHandle,
                               CLFFT_HERMITIAN_INTERLEAVED,
                               CLFFT_REAL));
    check_clFFT(clfftSetResultLocation(planHandle,
                                       CLFFT_OUTOFPLACE));

    size_t inStride[3] = {1, cM, cM*cN};
#if PAD_FIRST_DIM
    size_t outStride[3] = {1, padsize(M), padsize(M)*N};
#else
    size_t outStride[3] = {1, M, M*N};
#endif

    check_clFFT(clfftSetPlanInStride(planHandle,  CLFFT_3D, inStride));
    check_clFFT(clfftSetPlanOutStride(planHandle, CLFFT_3D, outStride));

    /* Bake the plan -- i.e. generate OpenCL kernel */
    check_clFFT(clfftBakePlan(planHandle, 1, // XXX
                              &clu->command_queue,
                              NULL, NULL));
    check_CL(clFinish(clu->command_queue));

    size_t clfft_buffersize = 0;
    check_clFFT(clfftGetTmpBufSize(planHandle, &clfft_buffersize));
    if(clfft_buffersize > 0 && clu->verbose > 2)
    {
        printf("Will require a buffer of size: %zu\n", clfft_buffersize);
    }
    clu_increase_clfft_buffer(clu, clfft_buffersize);
    return planHandle;
}
#endif

#ifndef VKFFT
clfftPlanHandle  gen_h2r_inplace_plan(clu_env_t * clu,
                                      size_t M, size_t N, size_t P)
{
    size_t cM = M_r2h(M);
    size_t cN = N;

    /* From Hermitian to Real of size MxNxP */
    clfftPlanHandle planHandle;

    size_t real_size[3] = {M, N, P};

    /* Start with a default plan and then adjust it */
    check_clFFT(clfftCreateDefaultPlan(&planHandle,
                                       clu->context,
                                       CLFFT_3D,
                                       real_size));

    check_clFFT(clfftSetPlanPrecision(planHandle, CLFFT_SINGLE));
    check_clFFT(clfftSetLayout(planHandle,
                               CLFFT_HERMITIAN_INTERLEAVED,
                               CLFFT_REAL));
    check_clFFT(clfftSetResultLocation(planHandle,
                                       CLFFT_INPLACE));

    size_t inStride[3] = {1, cM, cM*cN};
#if PAD_FIRST_DIM
    size_t outStride[3] = {1, padsize(M), padsize(M)*N};
#else
    printf("gen_h2r_inplace_plan got bad values\n");
    size_t outStride[3] = {1, M, M*N};
#endif

    check_clFFT(clfftSetPlanInStride(planHandle,  CLFFT_3D, inStride));
    check_clFFT(clfftSetPlanOutStride(planHandle, CLFFT_3D, outStride));

    /* Bake the plan -- i.e. generate OpenCL kernel */
    check_clFFT(clfftBakePlan(planHandle, 1, // XXX
                              &clu->command_queue,
                              NULL, NULL));
    check_CL(clFinish(clu->command_queue));

    size_t clfft_buffersize = 0;
    check_clFFT(clfftGetTmpBufSize(planHandle, &clfft_buffersize));
    if(clu->verbose > 1)
    {
        if(clfft_buffersize > 0)
        {
            printf("gen_h2r_inplace_plan: requires a buffer of size: %zu\n",
                   clfft_buffersize);
        } else {
            printf("gen_h2r_inplace_plan: no buffer required\n");
        }
    }
    clu_increase_clfft_buffer(clu, clfft_buffersize);
    return planHandle;
}
#endif


void clu_benchmark_transfer(clu_env_t * clu)
{

    size_t buf_size_nf = 1024*1024*1024;
    /* Got the same performance for the WRITE_ONLY and READ_ONLY */
    cl_mem_flags mem_flag = CL_MEM_READ_WRITE;

    printf(" -> Testing transfer rates using a buffer of %zu B\n", buf_size_nf*sizeof(float));
    uint8_t * buf = fim_malloc(buf_size_nf*sizeof(float));
    uint8_t * buf_copy = fim_malloc(buf_size_nf*sizeof(float));

    for(size_t kk = 0; kk<buf_size_nf; kk++)
    {
        buf[kk] = rand() % 256;
    }
    struct timespec t0, t1;
    cl_int ret;
    cl_mem buf_gpu = clCreateBuffer(clu->context,
                                    mem_flag,
                                    buf_size_nf*sizeof(float),
                                    NULL, &ret );
    clu->nb_allocated += buf_size_nf*sizeof(float);
    clu->n_alloc++;
    check_CL(ret);
    dw_gettime(&t0);
    check_CL( clEnqueueWriteBuffer( clu->command_queue,
                                    buf_gpu,
                                    CL_TRUE, // block
                                    0,
                                    buf_size_nf*sizeof(float), buf, 0,
                                    NULL, NULL));
    dw_gettime(&t1);
    float dt = clockdiff(&t0, &t1);
    printf("    Upload speed: %f GB/s\n", sizeof(float) * 1.0/dt);

    dw_gettime(&t0);
    check_CL( clEnqueueReadBuffer( clu->command_queue,
                                   buf_gpu,
                                   CL_TRUE, // block
                                   0,
                                   buf_size_nf*sizeof(float),
                                   buf_copy,
                                   0,
                                   NULL,
                                   NULL));
    dw_gettime(&t1);
    dt = clockdiff(&t0, &t1);
    printf("    Download speed: %f GB/s\n", sizeof(float)*1.0/dt);
    for(size_t kk = 0; kk<buf_size_nf; kk++)
    {
        if(buf[kk] != buf_copy[kk])
        {
            fprintf(stderr, "ERROR: downloaded data does not match the original data! buf[%zu] = %u buf_copy[%zu= %u]\n", kk, buf[kk], kk, buf_copy[kk]);
        }
    }

    fim_free(buf);
    fim_free(buf_copy);
    clu_release(clu, buf_gpu);
    return;
}


void fimcl_update_y(fimcl_t * gy, fimcl_t * image)
{
    here();

    /* Configure sizes */
    size_t nel = padsize(gy->M)*gy->N*gy->P;
    size_t localWorkSize; /* Use largest possible */
    cl_kernel kernel = gy->clu->update_y_kernel->kernel;

    /* This takes no time, probably just a look-up so no need to
     * cache it */
    check_CL(
             clGetKernelWorkGroupInfo(kernel,
                                      gy->clu->device_id,
                                      CL_KERNEL_WORK_GROUP_SIZE,
                                      sizeof(size_t), &localWorkSize, NULL)
             );

    /* The global size needs to be a multiple of the localWorkSize
     * i.e. it will be larger than the number of elements */
    size_t numWorkGroups = (nel + (localWorkSize -1) ) / localWorkSize;
    size_t globalWorkSize = localWorkSize * numWorkGroups;
    if(0)
    {
        printf("   globalWorkSize: %zu\n", globalWorkSize);
        printf("   localWorkSize: %zu\n", localWorkSize);
        printf("   numWorkGroups: %zu\n", numWorkGroups);
    }
    assert(globalWorkSize >= image->M*image->N*image->P);

    fimcl_sync(gy);
    //fimcl_sync(image); // never changed

    /* The image is assumed to be smaller or the same size
     * as the forward projection (which might be extended). */
    assert(image->M <= gy->M);
    assert(image->N <= gy->N);
    assert(image->P <= gy->P);

    check_CL( clSetKernelArg(kernel,
                             0, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &gy->buf) ); // argument value

    check_CL( clSetKernelArg(kernel,
                             1, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &image->buf) ); // argument value

    check_CL( clEnqueueNDRangeKernel(gy->clu->command_queue,
                                     kernel,
                                     1,
                                     NULL,
                                     &globalWorkSize, //global_work_size,
                                     &localWorkSize,
                                     0,
                                     NULL,
                                     &gy->wait_ev) );

    fimcl_sync(gy);

    return;
}


float fimcl_error_idiv(fimcl_t * forward, fimcl_t * image)
{
    int timers = 0;
    struct timespec tstart, tend, t0, t1;

    if(timers){
        printf("\n");
        dw_gettime(&tstart);
    }

    /* Number of voxels to process per work item
     */
    size_t nPerItem = 16;

    /* Performance considerations:
     * Re-use the partial-sums buffer? (avoid create and destroy at each iter)
     * Further reduction so that we don't have to download that much data?
     * Other reduction kernel?
     * */
    cl_int status = CL_SUCCESS;

    /* Configure sizes */
    size_t nel = padsize(forward->M)*forward->N*forward->P;
    size_t localWorkSize; /* Use largest possible */

    cl_kernel kernel = forward->clu->idiv_kernel->kernel;
    status = clGetKernelWorkGroupInfo(kernel,
                                      forward->clu->device_id,
                                      CL_KERNEL_WORK_GROUP_SIZE,
                                      sizeof(size_t), &localWorkSize, NULL);

    /* The global size needs to be a multiple of the localWorkSize
     * i.e. it will be larger than the number of elements */
    size_t numWorkGroups = (nel/nPerItem + (localWorkSize -1) ) / localWorkSize;
    numWorkGroups < 1 ? numWorkGroups = 1 : 0;

    size_t globalWorkSize = localWorkSize * numWorkGroups;
    if(0)
    {
        printf("   globalWorkSize: %zu\n", globalWorkSize);
        printf("   localWorkSize: %zu\n", localWorkSize); // Typically 256, hardware dependent
        printf("   numWorkGroups: %zu\n", numWorkGroups); // global/local
    }

    fimcl_sync(forward);
    fimcl_sync(image);

    /* The image is assumed to be smaller or the same size
     * as the forward projection (which might be extended). */
    assert(image->M <= forward->M);
    assert(image->N <= forward->N);
    assert(image->P <= forward->P);

    /* Create a buffer for the partial sums */
    if(timers){dw_gettime(&t0);}
    cl_mem partial_sums_gpu = clCreateBuffer(forward->clu->context,
                                             CL_MEM_WRITE_ONLY,
                                             numWorkGroups * sizeof(float),
                                             NULL,
                                             &status );
    forward->clu->n_alloc++;

    if(timers){
        dw_gettime(&t1);
        printf("Create buffer: %f \n", timespec_diff(&t1, &t0));
    }

    check_CL( clSetKernelArg(kernel,
                             0, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &forward->buf) ); // argument value

    check_CL( clSetKernelArg(kernel,
                             1, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &image->buf) ); // argument value

    check_CL( clSetKernelArg(kernel,
                             2, localWorkSize*sizeof(float), NULL) );

    check_CL( clSetKernelArg(kernel,
                             3, sizeof(cl_mem), &partial_sums_gpu));
    dw_gettime(&t0);
    check_CL( clEnqueueNDRangeKernel(forward->clu->command_queue,
                                     kernel,
                                     1,
                                     NULL,
                                     &globalWorkSize, //global_work_size,
                                     &localWorkSize,
                                     0,
                                     NULL,
                                     &forward->wait_ev) );

    fimcl_sync(forward);
    if(timers)
    {
        dw_gettime(&t1);
        printf("Reduce: %f \n", timespec_diff(&t1, &t0));
    }

    /* Download the result */
    if(timers){dw_gettime(&t0);}
    float * partial_sums = malloc(numWorkGroups*sizeof(float));
    status = clEnqueueReadBuffer(forward->clu->command_queue,
                                 partial_sums_gpu,
                                 CL_TRUE,
                                 0,
                                 numWorkGroups*sizeof(float),
                                 partial_sums,
                                 0,
                                 NULL,
                                 NULL);
    clFinish(forward->clu->command_queue);
    if(timers)
    {
        dw_gettime(&t1);
        printf("Download: %f \n", timespec_diff(&t1, &t0));
    }

    if(timers){dw_gettime(&t0);}
    clu_release(forward->clu, partial_sums_gpu);

    if(timers)
    {
        dw_gettime(&t1);
        printf("Release buffer: %f \n", timespec_diff(&t1, &t0));
    }

    float sum_gpu = 0;
#pragma omp parallel for reduction(+ : sum_gpu)
    for(size_t kk = 0 ; kk<numWorkGroups; kk++)
    {
        sum_gpu += partial_sums[kk];
    }
    free(partial_sums);

    if(timers)
    {
        dw_gettime(&tend);
        printf("Total: %f \n", timespec_diff(&tend, &tstart));
    }

    return sum_gpu / (float) (image->M*image->N*image->P);
}


float * pad_for_inplace(const float * X, size_t M, size_t N, size_t P)
{
    const size_t nchunk = N*P;
    const size_t chunk_size = M*sizeof(float);
    float * PX = fim_malloc(padsize(M)*N*P*sizeof(float));

    assert(PX != NULL);
    /* P values at a time */
    for(size_t c = 0; c < nchunk; c++)
    {
        memcpy(PX+c*padsize(M), X + c*M, chunk_size);
    }

    return PX;
}

float * unpad_from_inplace(const float * PX, size_t M, size_t N, size_t P)
{
    const size_t nchunk = N*P;
    const size_t chunk_size = M*sizeof(float);
    float * X = fim_malloc(M*N*P*sizeof(float));
    assert(PX != NULL);
    /* P values at a time */
    for(size_t c = 0; c < nchunk; c++)
    {
        memcpy(X + c*M, PX+c*padsize(M), chunk_size);
    }

    return X;
}
