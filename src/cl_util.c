#include "cl_util.h"
#include "clext.h"

#include "kernels/cl_complex_mul.h" // corresponds to cl_complex_mul
#include "kernels/cl_complex_mul_conj.h" // corresponds to cl_complex_mul_conj
#include "kernels/cl_complex_mul_inplace.h"
#include "kernels/cl_complex_mul_conj_inplace.h"
#include "kernels/cl_error_idiv.h"
#include "kernels/cl_real_mul_inplace.h"
#include "kernels/cl_positivity.h"
#include "kernels/cl_shb_update.h"


static char * read_program(const char * fname, size_t * size);
static double clockdiff(struct timespec* start, struct timespec * finish);

/* Number of floats */
static size_t fimcl_nreal(fimcl_t * X);
/* Number of complex values required for Hermitian representation */
static size_t fimcl_ncx(fimcl_t * X);

static size_t fimcl_hM(fimcl_t * X);
static size_t fimcl_hN(fimcl_t * X);
static size_t fimcl_hP(fimcl_t * X);

static size_t M_r2h(size_t M);
static size_t N_r2h(size_t N);
static size_t P_r2h(size_t P);

void clu_exit_error(cl_int err,
                    const char * file,
                    const char * function,
                    int line,
                    int clfft)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "ERROR!\n");
    fprintf(stderr, "There was an unrecoverable problem in %s (%s) at line %d\n",
            file, function, line);

    const char * err_cl = clGetErrorString(err);
    printf("OpenCl error=%s\n", err_cl);

    if(clfft)
    {
        if(strcmp(err_cl, "CL_UNKNOWN_ERROR") == 0)
        {
            const char * err_clfft = get_clfft_error_string(err);
            printf("clFFT error=%s\n", err_clfft);
        }
    }

    exit(EXIT_FAILURE);
}

#if 1
static size_t M_r2h(size_t M)
{
    return (1+M/2);
}
static size_t N_r2h(size_t N)
{
    return N;
}
static size_t P_r2h(size_t P)
{
    return P;
}
#endif

#if 0
static size_t M_r2h(size_t M)
{
    return M;
}
static size_t N_r2h(size_t N)
{
    return N;
}
static size_t P_r2h(size_t P)
{
    return (1+P/2);
}
#endif


static size_t fimcl_hM(fimcl_t * X)
{
    return M_r2h(X->M);
}
static size_t fimcl_hN(fimcl_t * X)
{
    return N_r2h(X->N);
}

static size_t fimcl_hP(fimcl_t * X)
{
    return P_r2h(X->P);
}


static size_t fimcl_ncx(fimcl_t * X)
{
    return fimcl_hM(X)*fimcl_hN(X)*fimcl_hP(X);
}

static size_t fimcl_nreal(fimcl_t * X)
{
    return X->M * X->N * X->P;
}

float * fimcl_download(fimcl_t * gX)
{
    assert(gX != NULL);
    fimcl_sync(gX);
    if(gX->transformed == 0)
    {
        size_t MNP = gX->M*gX->N*gX->P;
        if(gX->clu->verbose > 1)
        {
            printf("Downloading real data %zu x %zu x %zu (%zu floats)\n",
                   gX->M, gX->N, gX->P, MNP);
        }

        float * X = calloc(MNP, sizeof(float));
        if(X == NULL)
        {
            fprintf(stderr, "Failed to allocate memory in %s %s %d\n", __FILE__,
                    __FUNCTION__,
                    __LINE__);
            fprintf(stderr, "Tried to allocate %zu x %zu x %zu (%zu floats)\n",
                    gX->M, gX->N, gX->P, MNP);
            exit(EXIT_FAILURE);
        }
        assert(X != NULL);
        check_CL( clEnqueueReadBuffer( gX->clu->command_queue,
                                       gX->buf,
                                       CL_TRUE, // blocking
                                       0, // offset
                                       MNP * sizeof( float ), // size
                                       X,
                                       0,
                                       NULL,
                                       &gX->wait_ev ));
        fimcl_sync(gX); // todo remove
        return X;
    }
    if(gX->transformed == 1)
    {
        size_t cMNP = fimcl_ncx(gX);
        if(gX->clu->verbose > 1)
        {
            printf("Downloading transformed data %zu x %zu x %zu (%zu floats)\n",
                   fimcl_hM(gX), fimcl_hN(gX), fimcl_hP(gX), cMNP*2);
        }

        float * X = malloc(cMNP*2*sizeof(float));
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
        fimcl_sync(gX); // todo remove
        return X;
    }
    assert(0);
    return NULL;
}

fimcl_t * fimcl_new(clu_env_t * clu, int ffted, int fullsize,
                    const float * X, size_t M, size_t N, size_t P)
{
    // Consider CL_MEM_COPY_HOST_PTR on creation
    if(clu->verbose > 2)
    {
        printf("fimcl_new, ffted = %d, fullsize = %d\n", ffted, fullsize);;
    }
    fimcl_t * Y = calloc(1, sizeof(fimcl_t));
    Y->M = M;
    Y->N = N;
    Y->P = P;
    Y->transformed = ffted;
    Y->clu = clu;
    Y->wait_ev = CL_SUCCESS;

    cl_int ret;
    if(ffted == 0 && fullsize == 0)
    {
        Y->buf = clCreateBuffer(clu->context,
                                CL_MEM_READ_WRITE,
                                M*N*P* sizeof(float),
                                NULL, &ret );
        Y->buf_size_nf = M*N*P;
        clu->nb_allocated += M*N*P*sizeof(float);
    } else {
        Y->buf_size_nf = fimcl_ncx(Y)*2;
        Y->buf = clCreateBuffer(clu->context,
                                CL_MEM_READ_WRITE,
                                Y->buf_size_nf*sizeof(float),
                                NULL, &ret );
        clu->nb_allocated += Y->buf_size_nf*sizeof(float);
        Y->fullsize = 1;
    }
    check_CL(ret);

    if(X == NULL)
    {
        goto done;
    }

    if(ffted == 0)
    {

        if(fullsize)
        {
            /* Until we get the padding correct we fill the whole memory region,
             * which is larger than the image data, with zeros  */
            if(clu->verbose > 2)
            {
                printf("Clearing the clBuffer before uploading\n");
            }
            float zero_pattern[1] = {0.0};
            check_CL ( clEnqueueFillBuffer (clu->command_queue, //  command_queue ,
                                            Y->buf, // cl_mem  buffer ,
                                            &zero_pattern, // const void  *pattern ,
                                            sizeof(float), // size_t  pattern_size ,
                                            0, // size_t  offset ,
                                            Y->buf_size_nf*sizeof(float), // size_t  size ,
                                            0, // cl_uint  num_events_in_wait_list ,
                                            NULL, // const cl_event  *event_wait_list ,
                                            &Y->wait_ev)); // cl_event  *event
        }
        check_CL( clEnqueueWriteBuffer( clu->command_queue,
                                        Y->buf,
                                        CL_TRUE, // blocking_write
                                        0,
                                        M*N*P* sizeof( float ),
                                        X,
                                        0, // num_events_in_wait_list
                                        NULL,
                                        &Y->wait_ev ) );
        fimcl_sync(Y);
    }
    if(ffted == 1)
    {
        fprintf(stderr, "TODO: Unable to upload already ffted data in %s, line %d\n",
                __FILE__, __LINE__);
        exit(EXIT_FAILURE);
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
    fimcl_t * H = NULL;
    if(G->transformed == 0)
    {
        H = fimcl_new(G->clu, 0, G->fullsize, NULL, G->M, G->N, G->P);

        fimcl_sync(G);
        fimcl_sync(H);
        check_CL(clEnqueueCopyBuffer(
                                     G->clu->command_queue,//cl_command_queue command_queue,
                                     G->buf, //cl_mem src_buffer,
                                     H->buf, //cl_mem dst_buffer,
                                     0, //size_t src_offset,
                                     0, //size_t dst_offset,
                                     G->M*G->N*G->P*sizeof(float),//size_t size,
                                     0, //cl_uint num_events_in_wait_list,
                                     NULL, //const cl_event* event_wait_list,
                                     &H->wait_ev)); //cl_event* event);
        fimcl_sync(H);
    } else {
        fprintf(stderr, "fimcl_copy for transformed objects is not implemented \n");
        exit(EXIT_FAILURE);
    }
    return H;
}

void fimcl_free(fimcl_t * G)
{
    assert(G != NULL);
    fimcl_sync(G);
    check_CL(clReleaseMemObject(G->buf));
    G->clu->n_release++;
    free(G);
}


static fimcl_t * _fimcl_convolve(fimcl_t * X, fimcl_t * Y, int mode, int conj)
{
    struct timespec tk0, tk1, tfft0, tfft1, tifft0, tifft1;
    if(X->clu->verbose > 1)
    {
        printf("fimcl_convolve\n");
    }
    clFinish(X->clu->command_queue);
    if(X->transformed == 0)
    {
        if(X->clu->verbose > 1)
        {
            printf("fimcl_convolve - convolving first argument\n");
        }

        clock_gettime(CLOCK_MONOTONIC, &tfft0);
        if(X->fullsize)
        {
            fimcl_fft_inplace(X);
        } else {
            fimcl_fft(X);
        }
        clock_gettime(CLOCK_MONOTONIC, &tfft1);
    }
    clFinish(X->clu->command_queue);
    if(Y->transformed == 0)
    {
        if(X->clu->verbose > 1)
        {
            printf("fimcl_convolve - convolving 2nd argument\n");
        }
        clock_gettime(CLOCK_MONOTONIC, &tfft0);
        if(Y->fullsize)
        {
            fimcl_fft_inplace(Y);
        } else {
            fimcl_fft(Y);
        }
        clock_gettime(CLOCK_MONOTONIC, &tfft1);
    }
    clFinish(X->clu->command_queue);
    assert(X->transformed == 1);
    assert(Y->transformed == 1);

    if(mode == CLU_KEEP_ALL)
    {
        fimcl_t * Z = fimcl_new(X->clu, 1, 0, NULL, X->M, X->N, X->P);
        assert(Z->transformed == 1);

        clock_gettime(CLOCK_MONOTONIC, &tk1);
        fimcl_complex_mul(X, Y, Z, conj);

        clock_gettime(CLOCK_MONOTONIC, &tk1);
        fimcl_ifft(Z);

        return Z;
    }

    if(mode == CLU_DROP_ALL)
    {
        clock_gettime(CLOCK_MONOTONIC, &tk0);
        fimcl_complex_mul_inplace(Y, Y, conj);
        clock_gettime(CLOCK_MONOTONIC, &tk1);
        fimcl_free(Y);

        fimcl_t * fX = NULL;
        if(X->fullsize == 1)
        {
            clock_gettime(CLOCK_MONOTONIC, &tifft0);
            fimcl_ifft_inplace(X);
            fX = X;
            clock_gettime(CLOCK_MONOTONIC, &tifft1);
        } else {
            fX = fimcl_ifft(X);
            fimcl_free(X);
        }

        double dfft = clockdiff(&tfft0, &tfft1);
        double dtk = clockdiff(&tk0, &tk1);
        double difft = clockdiff(&tifft0, &tifft1);

        if(X->clu->verbose > 1)
        {
            printf("- One fft took %f s\n", dfft);
            printf("- Multiplication took %f s\n", dtk);
            printf("- Inverse fft took %f s\n", difft);
        }
        clFinish(X->clu->command_queue);
        return fX;
    }

    if(mode == CLU_KEEP_2ND)
    {
        clock_gettime(CLOCK_MONOTONIC, &tk0);
        fimcl_complex_mul_inplace(Y, X, conj);
        clock_gettime(CLOCK_MONOTONIC, &tk1);

        fimcl_t * fX = NULL;
        if(X->fullsize == 1)
        {
            clock_gettime(CLOCK_MONOTONIC, &tifft0);
            fimcl_ifft_inplace(X);
            fX = X;
            clock_gettime(CLOCK_MONOTONIC, &tifft1);
        } else {
            fX = fimcl_ifft(X);
            fimcl_sync(fX);
            fimcl_free(X);
        }

        double dtk = clockdiff(&tk0, &tk1);
        double difft = clockdiff(&tifft0, &tifft1);
        double dfft = clockdiff(&tfft0, &tfft1);
        if(fX->clu->verbose > 1)
        {
            printf("- One fft took %f s\n", dfft);
            printf("- Multiplication took %f s\n", dtk);
            printf("- Inverse fft took %f s\n", difft);
        }

        return fX;
    }

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
    check_CL( clu_wait_for_event(X->wait_ev, 10000));
    clFinish(X->clu->command_queue);
    return;
}

void fimcl_complex_mul(fimcl_t * X, fimcl_t * Y, fimcl_t * Z, int conj)
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
    assert(X->M == Y->M);
    assert(X->N == Y->N);
    assert(X->P == Y->P);
    assert(X->transformed == 0);
    assert(Y->transformed == 0);

    cl_kernel kernel = X->clu->kern_real_mul_inplace.kernel;

    /* Set sizes */
    size_t localWorkSize; /* Use largest possible */
    check_CL ( clGetKernelWorkGroupInfo(kernel,
                                        X->clu->device_id,
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &localWorkSize,
                                        NULL) );

    size_t nel = X->M * X->N * X->P;
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
                             2,
                             sizeof(cl_mem),
                             (void *) &X->clu->real_size) );

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


}

void fimcl_shb_update(fimcl_t * P, fimcl_t * X, fimcl_t * XP, float alpha)
{

    cl_kernel kernel = X->clu->kern_shb_update.kernel;

    /* Set sizes */
    size_t localWorkSize; /* Use largest possible */
    check_CL ( clGetKernelWorkGroupInfo(kernel,
                                        X->clu->device_id,
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &localWorkSize,
                                        NULL) );

    size_t nel = X->M * X->N * X->P;
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

    check_CL( clSetKernelArg(kernel,
                             3,
                             sizeof(cl_mem),
                             (void *) &X->clu->real_size) );

    check_CL( clEnqueueWriteBuffer( X->clu->command_queue,
                                    X->clu->float_gpu,
                                    CL_TRUE,
                                    0,
                                    sizeof(float),
                                    &alpha,
                                    0, NULL, NULL));

    check_CL( clSetKernelArg(kernel,
                             4,
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


}


void fimcl_positivity(fimcl_t * X, float val)
{
    assert(X != NULL);
    assert(X->transformed == 0);
    cl_kernel kernel = X->clu->kern_real_positivity.kernel;

    /* Set sizes */
    size_t localWorkSize; /* Use largest possible */
    check_CL ( clGetKernelWorkGroupInfo(kernel,
                                        X->clu->device_id,
                                        CL_KERNEL_WORK_GROUP_SIZE,
                                        sizeof(size_t),
                                        &localWorkSize,
                                        NULL) );

    size_t nel = fimcl_nreal(X); //X->M * X->N * X->P;
    size_t numWorkGroups = (nel + (localWorkSize -1) ) / localWorkSize;
    size_t globalWorkSize = localWorkSize * numWorkGroups;

    check_CL( clSetKernelArg(kernel,
                             0, // argument index
                             sizeof(cl_mem), // argument size
                             (void *) &X->buf) ); // argument value

    check_CL( clSetKernelArg(kernel,
                             1,
                             sizeof(cl_mem),
                             (void *) &X->clu->real_size) );

    check_CL( clEnqueueWriteBuffer( X->clu->command_queue,
                                    X->clu->float_gpu,
                                    CL_TRUE,
                                    0,
                                    sizeof(float),
                                    &val,
                                    0, NULL, NULL));

    check_CL( clSetKernelArg(kernel,
                             2,
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
}

void fimcl_complex_mul_inplace(fimcl_t * X, fimcl_t * Y, int conj)
{
    size_t global_work_offset[] = {0, 0, 0};
    size_t local_work_size[] = {1,1,1};

    cl_kernel kernel = X->clu->kern_mul_inplace.kernel;
    if(conj == 1)
    {
        kernel = X->clu->kern_mul_conj_inplace.kernel;
    }

    size_t cMNP = fimcl_ncx(X);

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
                                     global_work_offset,
                                     &cMNP, //global_work_size,
                                     local_work_size,
                                     0,
                                     NULL,
                                     &Y->wait_ev) );

    fimcl_sync(Y);
    return;
}


fimcl_t * fimcl_ifft(fimcl_t * fX)
{
    if(fX->clu->verbose > 1)
    {
        printf("fimcl_ifft\n");
    }
    assert(fX->transformed == 1);

    fimcl_t * X = calloc(1, sizeof(fimcl_t));
    X->M = fX->M;
    X->N = fX->N;
    X->P = fX->P;
    X->transformed = 0;
    X->wait_ev = CL_SUCCESS;
    X->fullsize = 0;
    X->clu = fX->clu;

    cl_int ret;
    X->buf = clCreateBuffer(fX->clu->context,
                            CL_MEM_READ_WRITE,
                            fX->M*fX->N*fX->P *  sizeof(float),
                            NULL, &ret );
    fX->clu->nb_allocated += fX->M*fX->N*fX->P *  sizeof(float);
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

    return X;
}

void fimcl_ifft_inplace(fimcl_t * X)
{
    if(X->clu->verbose > 1)
    {
        printf("fimcl_ifft_inplace\n");
    }
    assert(X->transformed == 1);

    fimcl_sync(X);
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
    clFinish(X->clu->command_queue);
    X->transformed = 0;
    fimcl_sync(X);
    return;
}


fimcl_t * fimcl_fft(fimcl_t * X)
{
    if(X->clu->verbose > 1)
    {
        printf("fimcl_fft\n");
    }
    assert(X->transformed == 0);
    fimcl_sync(X);

    /* Number of complex elements in the FFT */
    size_t cMNP = fimcl_ncx(X);

    fimcl_t * fX = calloc(1, sizeof(fimcl_t));
    fX->M = X->M;
    fX->N = X->N;
    fX->P = X->P;
    fX->transformed = 1;
    fX->clu = X->clu;

    cl_int ret;
    fX->buf_size_nf = cMNP * 2;
    fX->buf = clCreateBuffer(X->clu->context,
                             CL_MEM_READ_WRITE,
                             fX->buf_size_nf *  sizeof(float),
                             NULL, &ret );
    X->clu->nb_allocated += fX->buf_size_nf *  sizeof(float);
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
}

void fimcl_fft_inplace(fimcl_t * X)
{
    if(X->clu->verbose > 1)
    {
        printf("fimcl_fft_inplace\n");
    }
    assert(X->transformed == 0);

    /* Number of complex elements in the FFT */
    size_t cMNP = fimcl_ncx(X);

    if(X->buf_size_nf < cMNP * 2)
    {
        fprintf(stderr,
                "fimcl_fft_inplace: the object does not have a "
                "large enough buffer for in-place fft\n"
                "current size: %zu required size: %zu\n",
                X->buf_size_nf, cMNP*2);

        exit(EXIT_FAILURE);
    }

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
    X->transformed = 1;
    return;
}


float * clu_convolve(clu_env_t * clu,
                     int dropy,
                     float * X, float * Y,
                     size_t M, size_t N, size_t P)
{
    size_t cMNP = M_r2h(M)*N_r2h(N)*P_r2h(P);
    if(clu->verbose > 1)
    {
        size_t buf_mem = 2*cMNP*2*sizeof(float)/1000000;
        printf("clu_convolve: GPU memory needed for storage: %zu MB\n"
               "              Approximate total need: %zu\n",
               buf_mem, 15*buf_mem/10);
    }
    // Note full size specified for in-place transformations
    struct timespec t0, t1, t2;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    fimcl_t * gX = fimcl_new(clu, 0, 1, X, M, N, P);
    clock_gettime(CLOCK_MONOTONIC, &t1);
    fimcl_t * gY = fimcl_new(clu, 0, 1, Y, M, N, P);
    clock_gettime(CLOCK_MONOTONIC, &t2);

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
        free(Y);
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

    fprintf(fid, "CL_DEVICE_GLOBAL_MEM_SIZE = %lu (%lu MiB)\n",
            dev_memory,
            dev_memory/1000000);

    size_t driver_version_buff_size = 256;
    size_t driver_version_buff_size_used = 0;
    char * driver_version_buff = malloc(driver_version_buff_size);
    check_CL(clGetDeviceInfo(dev_id, // 	cl_device_id device,
                             CL_DRIVER_VERSION,
                             driver_version_buff_size,
                             driver_version_buff,
                             &driver_version_buff_size_used));

    fprintf(fid, "CL_DRIVER_VERSION = %s\n", driver_version_buff);
    free(driver_version_buff);
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

clu_env_t * clu_new(int verbose)
{
    clu_env_t * env = calloc(1, sizeof(clu_env_t));
    env->verbose = verbose;
    env->clfft_buffer = NULL;
    env->clfft_buffer_size = 0;
    env->clFFT_loaded = 0;

    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;

    check_CL( clGetPlatformIDs(1, &env->platform_id, &ret_num_platforms) );

    if(env->verbose > 1)
    {
        printf("Found %d CL platforms\n", ret_num_platforms);
    }
    assert(ret_num_platforms > 0);


    check_CL( clGetDeviceIDs(env->platform_id,
                             CL_DEVICE_TYPE_ALL,
                             1,
                             &env->device_id,
                             &ret_num_devices));

    if(env->verbose > 1)
    {
        printf("Found %d CL devices\n", ret_num_devices);
    }

    if(env->verbose > 1)
    {
        clu_print_device_info(stdout, env->device_id);
    }


    cl_int ret = CL_SUCCESS;
    env->context = clCreateContext( NULL,
                                    1,
                                    &env->device_id,
                                    NULL,
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

    return env;
}

void clu_prepare_fft(clu_env_t * clu, size_t M, size_t N, size_t P)
{
    if(clu->verbose > 1)
    {
        printf("Preparing for convolutions of size %zu x %zu x %zu\n",
               M, N, P);
    }
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
    clu->r2h_plan = gen_r2h_plan(clu, M, N, P);
    clu->r2h_inplace_plan = gen_r2h_inplace_plan(clu, M, N, P);

    clu->h2r_plan = gen_h2r_plan(clu, M, N, P);
    clu->h2r_inplace_plan = gen_h2r_inplace_plan(clu, M, N, P);

    /* Create the OpenCL kernels that we will use*/
    clu->kern_mul =
        *clu_kernel_new(clu,
                        NULL, //"cl_complex_mul.c",
                        (const char *) cl_complex_mul,
                        cl_complex_mul_len,
                        "cl_complex_mul");

    clu->kern_mul_conj =
        *clu_kernel_new(clu,
                        NULL, //"cl_complex_mul_conj.c",
                        (const char *) cl_complex_mul_conj,
                        cl_complex_mul_conj_len,
                        "cl_complex_mul_conj");
    clu->kern_mul_inplace =
        *clu_kernel_new(clu,
                        NULL,
                        (const char *) cl_complex_mul_inplace,
                        cl_complex_mul_inplace_len,
                        "cl_complex_mul_inplace");
    clu->kern_mul_conj_inplace =
        *clu_kernel_new(clu,
                        NULL,
                        (const char *) cl_complex_mul_conj_inplace,
                        cl_complex_mul_conj_inplace_len,
                        "cl_complex_mul_conj_inplace");

    clu->kern_real_mul_inplace =
        * clu_kernel_new(clu,
                         NULL,
                         (const char *) src_kernels_cl_real_mul_inplace_c,
                         src_kernels_cl_real_mul_inplace_c_len,
                         "cl_real_mul_inplace");
    clu->kern_real_positivity =
        * clu_kernel_new(clu,
                       NULL,
                       (const char *) src_kernels_cl_positivity_c,
                       src_kernels_cl_positivity_c_len,
                       "cl_positivity");

    clu->kern_shb_update =
        * clu_kernel_new(clu,
                         NULL,
                         (const char *) src_kernels_cl_shb_update_c,
                         src_kernels_cl_shb_update_c_len,
                         "cl_shb_update");


    /* Set up size-dependent data */
    size_t * size = malloc(4*sizeof(size_t));
    size[0] = M*N*P;
    size[1] = M;
    size[2] = M;
    size[3] = P;
    cl_int status;
    clu->real_size = clCreateBuffer(clu->context,
                                    CL_MEM_READ_ONLY,
                                    4*sizeof(size_t),
                                    NULL, &status);
    check_CL(status);
    check_CL( clEnqueueWriteBuffer( clu->command_queue,
                                    clu->real_size,
                                    CL_TRUE,
                                    0,
                                    4*sizeof(size_t),
                                    size,
                                    0, NULL, NULL));
    clu->n_alloc++;
    free(size);

    float value = 0;
    clu->float_gpu = clCreateBuffer(clu->context,
                                 CL_MEM_READ_ONLY,
                                 sizeof(float),
                                 NULL, &status);
    check_CL(status);
    check_CL( clEnqueueWriteBuffer( clu->command_queue,
                                    clu->float_gpu,
                                    CL_TRUE,
                                    0,
                                    sizeof(float),
                                    &value,
                                    0, NULL, NULL));
    clu->n_alloc++;


    clu->clFFT_loaded = 1;
    return;
}

void clu_destroy(clu_env_t * clu)
{
    if(clu->verbose > 1)
    {
        printf("Closing the OpenCL environment\n");
    }

    /* Release clFFT library. */
    if(clu->clFFT_loaded)
    {
        /* Teardown */
        check_clFFT(clfftDestroyPlan( &clu->r2h_plan ));
        check_clFFT(clfftDestroyPlan( &clu->h2r_plan ));
        check_clFFT(clfftTeardown());

        /* Free memory */
        if(clu->clfft_buffer_size > 0)
        {
            check_CL(clReleaseMemObject(clu->clfft_buffer));
        }
        if(clu->real_size != NULL)
        {
            check_CL(clReleaseMemObject(clu->real_size));
        }
        if(clu->float_gpu != NULL)
        {
            check_CL(clReleaseMemObject(clu->float_gpu));
        }
    }

    /* Clear up the OpenCL stuff */
    check_CL(clFlush(clu->command_queue));
    clu_kernel_destroy(clu->kern_mul);
    clu_kernel_destroy(clu->kern_mul_conj);
    check_CL(clReleaseCommandQueue(clu->command_queue) );
    check_CL(clReleaseContext(clu->context));

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

        fprintf(stderr,"Build log: \n%s\n", buff_erro); //Be careful with  the fprint
        free(buff_erro);
        fprintf(stderr,"clBuildProgram failed\n");
        return NULL;
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
    return;
}


cl_uint clu_wait_for_event(cl_event clev, size_t ns)
{
    /* Always use a clFinish afterwards since ...
       "Using clGetEventInfo to determine if a command identified by
       event has finished execution
       (i.e. CL_EVENT_COMMAND_EXECUTION_STATUS returns CL_COMPLETE) is
       not a synchronization point. There are no guarantees that the
       memory objects being modified by command associated with event
       will be visible to other enqueued commands".
    */

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
        nanosleep(&sleep_req,
                  &sleep_rem);
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
        while( t > 0 && (t % 11) == 0)
        {
            t/=11;
        }
        while( t > 0 && (t % 13) == 0)
        {
            t/=13;
        }
        if(t == 1)
        {
            return N;
        }
        N++;
    }
    return N;
}


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

clfftPlanHandle  gen_r2h_plan(clu_env_t * clu,
                              size_t M, size_t N, size_t P)
{
    size_t cM = M_r2h(M);
    size_t cN = N_r2h(N);

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
                                       CLFFT_OUTOFPLACE));


    size_t inStride[3] =  {1, M, M*N};
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
        if(clfft_buffersize > 0 && clu->verbose > 2)
        {
            printf("Will require a buffer of size: %zu\n", clfft_buffersize);
        }
    }
    clu_increase_clfft_buffer(clu, clfft_buffersize);

    return planHandle;
}

clfftPlanHandle  gen_r2h_inplace_plan(clu_env_t * clu,
                                      size_t M, size_t N, size_t P)
{
    size_t cM = M_r2h(M);
    size_t cN = N_r2h(N);;

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

    size_t inStride[3] = {1, M, M*N};
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

cl_int clu_increase_clfft_buffer(clu_env_t * clu, size_t req_size)
{
    cl_int ret = CL_SUCCESS;

    if(req_size <= clu->clfft_buffer_size)
    {
        return ret;
    }
    if(clu->clfft_buffer_size > 0)
    {
        check_CL(clReleaseMemObject( clu->clfft_buffer));
    }
    clu->clfft_buffer = clCreateBuffer(clu->context,
                                       CL_MEM_READ_WRITE,
                                       req_size,
                                       NULL, &ret );
    clu->nb_allocated += req_size;
    check_CL(ret);
    clu->clfft_buffer_size = req_size;

    return ret;
}

clfftPlanHandle  gen_h2r_plan(clu_env_t * clu,
                              size_t M, size_t N, size_t P)
{
    size_t cM = M_r2h(M);
    size_t cN = N_r2h(N);;

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
    size_t outStride[3] = {1, M, M*N};

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

clfftPlanHandle  gen_h2r_inplace_plan(clu_env_t * clu,
                                      size_t M, size_t N, size_t P)
{
    size_t cM = M_r2h(M);
    size_t cN = N_r2h(N);;

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
    size_t outStride[3] = {1, M, M*N};

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

void clu_benchmark_transfer(clu_env_t * clu)
{

    size_t buf_size_nf = 1024*1024*1024;
    /* Got the same performance for the WRITE_ONLY and READ_ONLY */
    cl_mem_flags mem_flag = CL_MEM_READ_WRITE;

    printf(" -> Testing transfer rates using a buffer of %zu B\n", buf_size_nf*sizeof(float));
    uint8_t * buf = calloc(buf_size_nf, sizeof(float));
    uint8_t * buf_copy = calloc(buf_size_nf, sizeof(float));

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
    check_CL(ret);
    clock_gettime(CLOCK_MONOTONIC, &t0);
    check_CL( clEnqueueWriteBuffer( clu->command_queue,
                                    buf_gpu,
                                    CL_TRUE, // block
                                    0,
                                    buf_size_nf*sizeof(float), buf, 0,
                                    NULL, NULL));
    clock_gettime(CLOCK_MONOTONIC, &t1);
    float dt = clockdiff(&t0, &t1);
    printf("    Upload speed: %f GB/s\n", sizeof(float) * 1.0/dt);

    clock_gettime(CLOCK_MONOTONIC, &t0);
    check_CL( clEnqueueReadBuffer( clu->command_queue,
                                   buf_gpu,
                                   CL_TRUE, // block
                                   0,
                                   buf_size_nf*sizeof(float),
                                   buf_copy,
                                   0,
                                   NULL,
                                   NULL));
    clock_gettime(CLOCK_MONOTONIC, &t1);
    dt = clockdiff(&t0, &t1);
    printf("    Download speed: %f GB/s\n", sizeof(float)*1.0/dt);
    for(size_t kk = 0; kk<buf_size_nf; kk++)
    {
        if(buf[kk] != buf_copy[kk])
        {
            fprintf(stderr, "ERROR: downloaded data does not match the original data! buf[%zu] = %u buf_copy[%zu= %u]\n", kk, buf[kk], kk, buf_copy[kk]);
        }
    }

    free(buf);
    free(buf_copy);
    check_CL(clReleaseMemObject(buf_gpu));
    return;
}
