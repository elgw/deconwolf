#include "cl_util.h"
#include "clext.h"

#include "kernels/cl_complex_mul.h" // corresponds to cl_complex_mul
#include "kernels/cl_complex_mul_conj.h" // corresponds to cl_complex_mul_conj

static char * read_program(const char * fname, size_t * size);
static double clockdiff(struct timespec* start, struct timespec * finish);

void clu_exit_error(cl_int err,
                    const char * file,
                    const char * function,
                    int line,
                    int clfft)
{
    fprintf(stderr, "\n");
    fprintf(stderr, "ERROR!\n");
    fprintf(stderr, "There was an unrecoverable problem in %s/%s at line %d\n",
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


float * fimcl_download(fimcl_t * gX)
{
    size_t MNP = gX->M*gX->N*gX->P;
    float * X = malloc(MNP*sizeof(float));
    assert(X != NULL);
    assert(gX->transformed == 0);

    check_CL( clEnqueueReadBuffer( gX->clu->command_queue,
                                   gX->buf,
                                   CL_TRUE,
                                   0,
                                   MNP * sizeof( float ),
                                   X,
                                   0,
                                   NULL,
                                   &gX->wait_ev ));

    fimcl_sync(gX);
    return X;
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
                                M*N*P *  sizeof(float),
                                NULL, &ret );
        Y->buf_size = M*N*P;
    } else {
        Y->buf = clCreateBuffer(clu->context,
                                CL_MEM_READ_WRITE,
                                (1+M/2)*N*P *2*  sizeof(float),
                                NULL, &ret );
        Y->buf_size = (1+M/2)*N*P*2;
        Y->fullsize = 1;
    }
    check_CL(ret);

    if(X == NULL)
    {
        goto done;
    }

    if(ffted == 0)
    {

        check_CL( clEnqueueWriteBuffer( clu->command_queue,
                                        Y->buf,
                                        CL_FALSE, // blocking_write
                                        0,
                                        M*N*P* sizeof( float ),
                                        X,
                                        0, // num_events_in_wait_list
                                        NULL,
                                        &Y->wait_ev ) );
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
        printf("fimcl_convolve\n");
    }
    fimcl_t * H = NULL;
    if(G->transformed == 0)
    {
        H = fimcl_new(G->clu, 0, G->fullsize, NULL, G->M, G->N, G->P);
        fimcl_sync(G);
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
    } else {
        fprintf(stderr, "fimcl_copy for transformed objects is not implemented \n");
        exit(EXIT_FAILURE);
    }
    return H;
}

void fimcl_free(fimcl_t * G)
{
    clReleaseMemObject(G->buf);
    free(G);
}

fimcl_t * fimcl_convolve(fimcl_t * X, fimcl_t * Y, int mode)
{
    struct timespec tk0, tk1, tfft0, tfft1, tifft0, tifft1;
    if(X->clu->verbose > 1)
    {
        printf("fimcl_convolve\n");
    }

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

    assert(X->transformed == 1);
    assert(Y->transformed == 1);

    if(mode == CLU_KEEP_ALL)
    {
        fimcl_t * Z = fimcl_new(X->clu, 1, 0, NULL, X->M, X->N, X->P);
        assert(Z->transformed == 1);

        clock_gettime(CLOCK_MONOTONIC, &tk1);
        fimcl_complex_mul(X, Y, Z, 0);

        clock_gettime(CLOCK_MONOTONIC, &tk1);
        fimcl_ifft(Z);

        return Z;
    }

    if(mode == CLU_DROP_ALL)
    {
        clock_gettime(CLOCK_MONOTONIC, &tk0);
        fimcl_complex_mul(X, Y, X, 0);
        fimcl_sync(X);
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
        return fX;

    }

    if(mode == CLU_KEEP_2ND)
    {
        clock_gettime(CLOCK_MONOTONIC, &tk0);
        fimcl_complex_mul(X, Y, X, 0);
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

void fimcl_sync(fimcl_t * X)
{
    check_CL( clu_wait_for_event(X->wait_ev, 10000));
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

    size_t cMNP = (1+X->M/2)*X->N*X->P;

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
    check_CL(ret);

    /* And back again */
    check_clFFT(clfftEnqueueTransform(fX->clu->h2r_plan, // XXX
                                      CLFFT_BACKWARD,
                                      1, // numQueuesAndEcents
                                      &fX->clu->command_queue, // commQueues
                                      0, // numWaitEvents
                                      NULL, // waitEvents
                                      &fX->wait_ev, // cl_event * outEvents
                                      &fX->buf, // Input buffer
                                      &X->buf, // output buffer
                                      fX->clu->clfft_buffer)); // temp buffer
    fimcl_sync(fX);

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
                                      1, // numQueuesAndEcents
                                      &X->clu->command_queue, // commQueues
                                      0, // numWaitEvents
                                      NULL, // waitEvents
                                      &X->wait_ev, // cl_event * outEvents
                                      &X->buf, // Input buffer
                                      &X->buf, // output buffer
                                      NULL)); // temp buffer

    X->transformed = 0;
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

    /* Number of elements in real domain */
    size_t M = X->M;
    size_t N = X->N;
    size_t P = X->P;

    /* Number of complex elements in the FFT */
    size_t cM = 1 + M/2;
    size_t cN = N;
    size_t cP = P;
    size_t cMNP = cM*cN*cP;

    fimcl_t * fX = calloc(1, sizeof(fimcl_t));
    fX->M = X->M;
    fX->N = X->N;
    fX->P = X->P;
    fX->transformed = 1;
    fX->clu = X->clu;

    cl_int ret;
    fX->buf_size = cMNP * 2;
    fX->buf = clCreateBuffer(X->clu->context,
                             CL_MEM_READ_WRITE,
                             fX->buf_size *  sizeof(float),
                             NULL, &ret );
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
    return fX;
}

void fimcl_fft_inplace(fimcl_t * X)
{
    if(X->clu->verbose > 1)
    {
        printf("fimcl_fft_inplace\n");
    }
    assert(X->transformed == 0);


    /* Number of elements in real domain */
    size_t M = X->M;
    size_t N = X->N;
    size_t P = X->P;

    /* Number of complex elements in the FFT */
    size_t cM = 1 + M/2;
    size_t cN = N;
    size_t cP = P;
    size_t cMNP = cM*cN*cP;

    if(X->buf_size < cMNP * 2)
    {
        fprintf(stderr,
                "fimcl_fft_inplace: the object does not have a "
                "large enough buffer for in-place fft\n"
                "current size: %zu required size: %zu\n",
                X->buf_size, cMNP*2);

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
                                      &X->buf, // output buffer
                                      NULL)); // temp buffer

    X->transformed = 1;
    return;
}


float * clu_convolve(clu_env_t * clu,
                     int dropy,
                     float * X, float * Y,
                     size_t M, size_t N, size_t P)
{
    size_t cMNP = (1+M/2)*N*P;
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
    env->clfft_buffer_size = 0;
    env->clFFT_loaded = 0;
    cl_uint ret_num_devices;
    cl_uint ret_num_platforms;

    check_CL( clGetPlatformIDs(1, &env->platform_id, &ret_num_platforms) );


    if(env->verbose > 0)
    {
        printf("Found %d CL platforms\n", ret_num_platforms);
    }
    assert(ret_num_platforms > 0);


    check_CL( clGetDeviceIDs(env->platform_id,
                             CL_DEVICE_TYPE_ALL,
                             1,
                             &env->device_id,
                             &ret_num_devices));

    if(env->verbose > 0)
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
    if(env->verbose > 0)
    {
        printf("CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS=%u\n",
               CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS);
    }

    return env;
}

void clu_prepare_fft(clu_env_t * env, size_t M, size_t N, size_t P)
{
    assert(env->clFFT_loaded == 0);

    /* Setup clFFT. */
    check_CL( clfftInitSetupData(&env->fftSetup));
    check_CL( clfftSetup(&env->fftSetup));

    if(env->verbose > 1)
    {
        cl_uint major, minor, patch;
        check_clFFT(clfftGetVersion(&major, &minor, &patch));
        printf("clFFT version %d.%d.%d\n", major, minor, patch);
    }

    /* Create clFFT plans for this particular size */
    env->r2h_plan = gen_r2h_plan(env, M, N, P);
    env->r2h_inplace_plan = gen_r2h_inplace_plan(env, M, N, P);
    env->h2r_plan = gen_h2r_plan(env, M, N, P);
    env->h2r_inplace_plan = gen_h2r_inplace_plan(env, M, N, P);

    /* Create the OpenCL kernels that we will use*/
    env->kern_mul = *clu_kernel_new(env,
                                    NULL, //"cl_complex_mul.c",
                                    (const char *) cl_complex_mul,
                                    cl_complex_mul_len,
                                    "cl_complex_mul");
    env->kern_mul_conj = *clu_kernel_new(env,
                                    NULL, //"cl_complex_mul_conj.c",
                                    (const char *) cl_complex_mul_conj,
                                    cl_complex_mul_conj_len,
                                    "cl_complex_mul_conj");
    env->clFFT_loaded = 1;
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
        clfftTeardown( );

        /* Free memory */
        if(clu->clfft_buffer_size > 0)
        {
            clReleaseMemObject(clu->clfft_buffer);
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
                                 NULL,
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
        if (errcode) {
            printf("clGetProgramBuildInfo failed at line %d\n", __LINE__);
            exit(-3);
        }

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
    size_t cM = 1 + M/2;
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
                                       CLFFT_OUTOFPLACE));


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
    size_t cM = 1 + M/2;
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
    if(clu->verbose > 2)
    {
        printf("clu_increase_clfft_buffer: warning -- function disabled for testing\n");
    }
    return ret;

    if(req_size <= clu->clfft_buffer_size)
    {
        return ret;
    }
    clReleaseMemObject( clu->clfft_buffer);
    clu->clfft_buffer = clCreateBuffer(clu->context,
                                       CL_MEM_READ_WRITE,
                                       req_size,
                                       NULL, &ret );
    check_CL(ret);
    return ret;
}

clfftPlanHandle  gen_h2r_plan(clu_env_t * clu,
                              size_t M, size_t N, size_t P)
{
    size_t cM = 1 + M/2;
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
    return planHandle;
}

clfftPlanHandle  gen_h2r_inplace_plan(clu_env_t * clu,
                                      size_t M, size_t N, size_t P)
{
    size_t cM = 1 + M/2;
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
            printf("gen_h2r_inplace_plan: requires a buffer of size: %zu\n", clfft_buffersize);
        } else {
            printf("gen_h2r_inplace_plan: no buffer required\n");
        }
    }
    return planHandle;
}

void clu_benchmark_transfer(clu_env_t * clu)
{

    size_t buf_size = 1024*1024*1024;
    /* Got the same performance for the WRITE_ONLY and READ_ONLY */
    cl_mem_flags mem_flag = CL_MEM_READ_WRITE;

    printf(" -> Testing transfer rates using a buffer of %zu B\n", buf_size);
    uint8_t * buf = calloc(buf_size, 1);
    struct timespec t0, t1;
    cl_int ret;
    cl_mem buf_gpu = clCreateBuffer(clu->context,
                                 mem_flag,
                                 buf_size,
                                 NULL, &ret );
    check_CL(ret);
    clock_gettime(CLOCK_MONOTONIC, &t0);
    check_CL( clEnqueueWriteBuffer( clu->command_queue, buf_gpu, CL_TRUE, 0,
                                    buf_size, buf, 0,
                                    NULL, NULL));
    clock_gettime(CLOCK_MONOTONIC, &t1);
    float dt = clockdiff(&t0, &t1);
    printf("    Upload speed: %f GB/s\n", 1.0/dt);

    clock_gettime(CLOCK_MONOTONIC, &t0);
    check_CL( clEnqueueReadBuffer( clu->command_queue,
                                   buf_gpu,
                                   CL_TRUE,
                                   0,
                                   buf_size,
                                   buf,
                                   0,
                                   NULL,
                                   NULL));
    clock_gettime(CLOCK_MONOTONIC, &t1);
    dt = clockdiff(&t0, &t1);
    printf("    Download speed: %f GB/s\n", 1.0/dt);
    free(buf);
    clReleaseMemObject(buf_gpu);
    return;
}
