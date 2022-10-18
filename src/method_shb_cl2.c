#include "method_shb_cl.h"

//#define here(x) printf("%s %s %d\n", __FILE__, __FUNCTION__, __LINE__);
#define here(x) ;


#include "kernels/cl_idiv_kernel.h"
#include "kernels/cl_update_y_kernel.h"

clu_kernel_t * idiv_kernel;
clu_kernel_t * reduction_kernel;
clu_kernel_t * update_y_kernel;


void fimcl_to_tiff(fimcl_t * I_gpu, char * filename)
{
    float * I = fimcl_download(I_gpu);
    printf("Writing to %s\n", filename);
    fim_tiff_write_float(filename, I, NULL, I_gpu->M, I_gpu->N, I_gpu->P);
    fftwf_free(I);
}

void prepare_kernels(clu_env_t * clu,
                     size_t M, size_t N, size_t P,
                     size_t wM, size_t wN, size_t wP)
{
    char * argstring = malloc(1024);
    sprintf(argstring,
            "-D NELEMENTS=%zu "
            "-D M=%zu -D N=%zu -D P=%zu "
            "-D wM=%zu -D wN=%zu -D wP=%zu",
            wM*wN*wP, M, N, P, wM, wN, wP);
    idiv_kernel = clu_kernel_newa(clu,
                                  "kernels/cl_idiv_kernel.c",
                                  (const char *) cl_idiv_kernel,
                                  cl_idiv_kernel_len,
                                  "idiv_kernel",
                                  argstring);
    update_y_kernel = clu_kernel_newa(clu,
                                      "kernels/cl_update_y_kernel.c",
                                      (const char *) cl_update_y_kernel,
                                      cl_update_y_kernel_len,
                                      "update_y_kernel",
                                      argstring);

    free(argstring);
}

float fimcl_error_idiv(fimcl_t * forward, fimcl_t * image);

void gpu_update_y(fimcl_t * gy, fimcl_t * image)
{
    here();

    /* Configure sizes */
    size_t nel = gy->M*gy->N*gy->P;
    size_t localWorkSize; /* Use largest possible */
    cl_kernel kernel = update_y_kernel->kernel;
    // TOOD: does this take time?
    cl_int status = clGetKernelWorkGroupInfo(kernel,
                                             gy->clu->device_id,
                                             CL_KERNEL_WORK_GROUP_SIZE,
                                             sizeof(size_t), &localWorkSize, NULL);
    check_CL(status);
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
    cl_int status = CL_SUCCESS;

    /* Configure sizes */
    size_t nel = forward->M*forward->N*forward->P;
    size_t localWorkSize; /* Use largest possible */
    // TOOD: does this take time?
    status = clGetKernelWorkGroupInfo(idiv_kernel->kernel,
                                      forward->clu->device_id,
                                      CL_KERNEL_WORK_GROUP_SIZE,
                                      sizeof(size_t), &localWorkSize, NULL);
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

    fimcl_sync(forward);
    fimcl_sync(image);

    /* The image is assumed to be smaller or the same size
     * as the forward projection (which might be extended). */
    assert(image->M <= forward->M);
    assert(image->N <= forward->N);
    assert(image->P <= forward->P);

    /* Create a buffer for the partial sums */
    cl_mem partial_sums_gpu = clCreateBuffer(forward->clu->context,
                                             CL_MEM_WRITE_ONLY,
                                             numWorkGroups * sizeof(float),
                                             NULL,
                                             &status );


    cl_kernel kernel = idiv_kernel->kernel;
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

    /* Download the result */
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

    float sum_gpu = 0;
#pragma omp parallel for reduction(+ : sum_gpu)
    for(size_t kk = 0 ; kk<numWorkGroups; kk++)
    {
        sum_gpu += partial_sums[kk];
    }
    free(partial_sums);
    clReleaseMemObject(partial_sums_gpu);

    return sum_gpu / (float) (image->M*image->N*image->P);
}





static fimcl_t  * initial_guess_cl(clu_env_t * clu,
                                   const int64_t M, const int64_t N, const int64_t P,
                                   const int64_t wM, const int64_t wN, const int64_t wP)
{
    /* Create initial guess: the fft of an image that is 1 in MNP and 0 outside
     * M, N, P is the dimension of the microscopic image
     *
     * Possibly more stable to use the mean of the input image rather than 1
     */

    assert(wM >= M); assert(wN >= N); assert(wP >= P);

    float * one = fim_zeros(wM*wN*wP);

#pragma omp parallel for shared(one)
    for(int64_t cc = 0; cc < P; cc++) {
        for(int64_t bb = 0; bb < N; bb++) {
            for(int64_t aa = 0; aa < M; aa++) {
                one[aa + wM*bb + wM*wN*cc] = 1;
            }
        }
    }

    fimcl_t * gOne = fimcl_new(clu, 0, 0, one, wM, wN, wP);
    here();
    fimcl_sync(gOne);
    here();
    fftwf_free(one);
    fimcl_t * gfOne = fimcl_fft(gOne);

    return gfOne;
}


float iter_shb_cl2(fimcl_t ** _xp_gpu, // Output, f_(t+1)
                   fimcl_t * restrict im_gpu, // as above, on GPU
                   fimcl_t * restrict cK, // fft(psf)
                   fimcl_t * restrict pk_gpu, // p_k, estimation of the gradient
                   fimcl_t * W_gpu, // Bertero Weights
                   const dw_opts * s)
{
    fimcl_t * xp_gpu = _xp_gpu[0];

    here();
    fimcl_t * fft_pk_gpu = fimcl_fft(pk_gpu); // Pk -> fft_pk_gpu
    here();
    clFinish(fft_pk_gpu->clu->command_queue);
    here();
    putdot(s);
    fimcl_t * gy = fimcl_convolve(fft_pk_gpu, cK, CLU_KEEP_2ND);
    fft_pk_gpu = NULL;

    here();

    float error = fimcl_error_idiv(gy, im_gpu);

    putdot(s);

    gpu_update_y(gy, im_gpu);
    fimcl_t * _Y = gy;
    fimcl_sync(_Y);

    fimcl_t * Y = fimcl_fft(_Y);
    fimcl_free(_Y);
    clFinish(Y->clu->command_queue);

    fimcl_t * gx = fimcl_convolve_conj(Y, cK, CLU_KEEP_2ND);
    here();
    fimcl_real_mul_inplace(pk_gpu, gx); // gx = gpk .* gx
    here();
    fimcl_free(pk_gpu);
    here();
    if(W_gpu != NULL)
    {
        fimcl_real_mul_inplace(W_gpu, gx);
    }

    here();
    fimcl_positivity(gx, s->bg);
    here();
    fimcl_free(xp_gpu);
    _xp_gpu[0] = gx;
    return error;
}


float * deconvolve_shb_cl2(float * restrict im,
                           const int64_t M, const int64_t N, const int64_t P,
                           float * restrict psf,
                           const int64_t pM, const int64_t pN, const int64_t pP,
                           dw_opts * s)
{

    if(s->verbosity > 3)
    {
        printf("deconv_w debug info:\n");
        printf("Will deconvolve an image of size %" PRId64 "x%" PRId64 "x%" PRId64 "\n", M, N, P);
        printf("Using a PSF of size %" PRId64 "x%" PRId64 "x%" PRId64 "\n", pM, pN, pP);
        printf("Output will be of size %" PRId64 "x%" PRId64 "x%" PRId64 "\n", M, N, P);
        printf("psf will be freed\n");
    }


    if(s->verbosity > 1)
    {
        printf("Deconvolving\n");
    }

    if(s->nIter == 0)
    {
        fftw_free(psf);
        return fim_copy(im, M*N*P);
    }

    if(fim_maxAtOrigo(psf, pM, pN, pP) == 0)
    {
        if(s->verbosity > 0)
        {
            printf(" ! SHBCL: PSF is not centered!\n");
        }
        if(s->log != stdout)
        {
            fprintf(s->log, " ! SHBCL: PSF is not centered!\n");
        }
    }


    /* This is the work dimensions, i.e., dimensions
     * that will be used for all FFTs
     */

    int64_t wM = M + pM -1;
    int64_t wN = N + pN -1;
    int64_t wP = P + pP -1;

    if(s->borderQuality == 1)
    {
        wM = M + (pM+1)/2;
        wN = N + (pN+1)/2;
        wP = P + (pP+1)/2;
    }

    if(s->borderQuality == 0)
    {
        wM = int64_t_max(M, pM);
        wN = int64_t_max(N, pN);
        wP = int64_t_max(P, pP);
    }

    if(wP %2 == 1)
    {
        /*  The jobb size is not divisable by 2.
         * potentially it could be faster to extend the job by
         * one slice in Z, but for some sizes, e.g. 85->86 it gets slower
         * some benchmarking or prediction should guide this, not just a flag.
         */
        if(s->experimental1)
        {
            printf("      Adding one extra slice to the job for performance"
                   " this is still experimental. \n");
            // One slice of the "job" should be cleared every iterations.
            wP++;
        }
    }

    clu_env_t * clu = clu_new(s->verbosity);
    wM = clu_next_fft_size(wM);
    wN = clu_next_fft_size(wN);
    wP = clu_next_fft_size(wP);
    clu_prepare_fft(clu, wM, wN, wP);

    prepare_kernels(clu, M, N, P, wM, wN, wP);

    /* Total number of pixels */
    size_t wMNP = wM*wN*wP;

    if(s->verbosity > 0)
    { printf("image: [%" PRId64 "x%" PRId64 "x%" PRId64 "], psf: [%" PRId64 "x%" PRId64 "x%" PRId64 "], job: [%" PRId64 "x%" PRId64 "x%" PRId64 "]\n",
             M, N, P, pM, pN, pP, wM, wN, wP);
    }
    if(s->verbosity > 1)
    {
        // TODO
        // printf("Estimated peak memory usage: %.1f GB\n", wMNP*35.0/1e9);
    }
    fprintf(s->log, "image: [%" PRId64 "x%" PRId64 "x%" PRId64 "]\n"
            "psf: [%" PRId64 "x%" PRId64 "x%" PRId64 "]\n"
            "job: [%" PRId64 "x%" PRId64 "x%" PRId64 "] (%zu voxels)\n",
            M, N, P, pM, pN, pP, wM, wN, wP, wMNP);
    fflush(s->log);

    if(s->verbosity > 0)
    {
        printf("Iterating "); fflush(stdout);
    }

    // cK : "full size" fft of the PSF
    float * Z = fftwf_malloc(wMNP*sizeof(float));
    memset(Z, 0, wMNP*sizeof(float));
    /* Insert the psf into the bigger Z */
    fim_insert(Z, wM, wN, wP,
               psf, pM, pN, pP);

    free(psf);

    /* Shift the PSF so that the mid is at (0,0,0) */
    int64_t midM, midN, midP = -1;
    fim_argmax(Z, wM, wN, wP, &midM, &midN, &midP);
    //printf("max(PSF) = %f\n", fim_max(Z, wM*wN*wP));
    //printf("PSF(%ld, %ld, %ld) = %f\n", midM, midN, midP, Z[midM + wM*midN + wM*wN*midP]);
    fim_circshift(Z, wM, wN, wP, -midM, -midN, -midP);
    if(Z[0] != fim_max(Z, wM*wN*wP))
    {
        printf("Something went wrong here, the max of the PSF is in the wrong place\n");
        exit(1);
    }

    if(s->fulldump)
    {
        printf("Dumping to fullPSF.tif\n");
        fim_tiff_write_float("fullPSF.tif", Z, NULL, wM, wN, wP);
    }

    fimcl_t * im_gpu = fimcl_new(clu, 0, 0, im, M, N, P);

    //fftwf_complex * fftPSF = fft(Z, wM, wN, wP);
    fimcl_t * _clfftPSF = fimcl_new(clu, 0, 0,
                                    Z, wM, wN, wP);
    // fimcl_fft_inplace(clfftPSF);
    fimcl_t * clfftPSF = fimcl_fft(_clfftPSF);
    fimcl_free(_clfftPSF);

    fftwf_free(Z);

    putdot(s);

    float * W = NULL;
    fimcl_t * W_gpu = NULL;
    /* Sigma in Bertero's paper, introduced for Eq. 17 */
    if(s->borderQuality > 0)
    {
        //fftwf_complex * F_one = initial_guess(M, N, P, wM, wN, wP);
        here();
        fimcl_t * clfftOne = initial_guess_cl(clu, M, N, P, wM, wN, wP);
        here();
        //float * P1 = fft_convolve_cc_conj_f2(fftPSF, F_one, wM, wN, wP);
        fimcl_t * gclP1 = fimcl_convolve_conj(clfftOne, clfftPSF, CLU_KEEP_2ND);
        here();
        clfftOne = NULL; // freed already
        here();
        float * clP1 = fimcl_download(gclP1);
        here();
        //fim_tiff_write_float("P1.tif", P1, NULL, wM, wN, wP);
        //fim_tiff_write_float("P2.tif", clP1, NULL, wM, wN, wP);


        float sigma = 0.01;
#pragma omp parallel for shared(clP1)
        for(size_t kk = 0; kk<wMNP; kk++)
        {
            if(clP1[kk] > sigma)
            {
                clP1[kk] = 1/clP1[kk];
            } else {
                clP1[kk] = 0;
            }
        }
        W = clP1;
        W_gpu = fimcl_new(clu, 0, 0, W, wM, wN, wP);
        fftwf_free(W);
    }

    float sumg = fim_sum(im, M*N*P);

    /* x is the initial guess, initially the previous iteration,
     *  xp is
     *  set to be the same */

    float * x = fim_constant(wMNP, sumg / (float) wMNP);
    float * xp = fim_copy(x, wMNP);

    /* TODO: we only need two of these to be allocated at a time */
    fimcl_t * x_gpu = fimcl_new(clu, 0, 1, x, wM, wN, wP);


    fimcl_t * xp_gpu = fimcl_new(clu, 0, 1, xp, wM, wN, wP);
    fimcl_t * p_gpu = NULL;

    fftwf_free(x);
    fftwf_free(xp);

    dw_iterator_t * it = dw_iterator_new(s);
    while(dw_iterator_next(it) >= 0)
    {
        here();

        /* Eq. 10 in SHB paper */
        double alpha = ((float) it->iter-1.0)/((float) it->iter+2.0);
        //alpha /= 1.5;
        alpha < 0 ? alpha = 0: 0;

        here();
        /* To be interpreted as p^k in Eq. 7 of SHB
         *  p[kk] = x[kk] + alpha*(x[kk]-xp[kk]); */
        p_gpu = fimcl_new(clu, 0, 1, NULL, wM, wN, wP);

        fimcl_shb_update(p_gpu, x_gpu, xp_gpu, alpha);


        putdot(s);
        here();
        double err =
            iter_shb_cl2(
                         &xp_gpu, // xp is updated to the next guess
                                  //im,
                         im_gpu,
                         clfftPSF, // FFT of PSF
                         p_gpu, // Current guess
                         //p,
                         W_gpu, // Weights (to handle boundaries)
                         s);

        //fimcl_to_tiff(xp_gpu, "xp_gpu.tif");
        //exit(EXIT_FAILURE);

        here();
        //fimcl_free(p_gpu);
        //here();

        dw_iterator_set_error(it, err);
        {
            /* Swap so that the current is named x */
            fimcl_t * t = x_gpu;
            x_gpu = xp_gpu;
            xp_gpu = t;
        }

        putdot(s);
        dw_iterator_show(it, s);
        benchmark_write(s, it->iter, it->error, x, M, N, P, wM, wN, wP);
        here();
    } /* End of main loop */

    dw_iterator_free(it);
    here();
    fimcl_free(im_gpu);

    {
        /* Swap back so that x is the final iteration */
        fimcl_t * t = x_gpu;
        x_gpu = xp_gpu;
        xp_gpu = t;
    }
    here();

    if(s->verbosity > 0) {
        printf("\n");
    }
    here();
    fimcl_free(clfftPSF);

    if(s->fulldump)
    {
        printf("Dumping to fulldump.tif\n");
        fim_tiff_write("fulldump.tif", x, NULL, wM, wN, wP);
    }
    here();
    float * out_full = fimcl_download(x_gpu);
    fimcl_sync(x_gpu);
    fimcl_free(x_gpu);
    here();
    float * out = fim_subregion(out_full, wM, wN, wP, M, N, P);
    here();
    fftwf_free(out_full);
    here();

    clu_destroy(clu);
    return out;
}
