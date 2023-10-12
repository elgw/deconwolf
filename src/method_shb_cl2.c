#include "method_shb_cl2.h"

#define use_inplace_clfft 0


//#define here(x) printf("%s %s %d\n", __FILE__, __FUNCTION__, __LINE__);
#define here(x) ;


static void fimcl_to_tiff(fimcl_t * I_gpu, char * filename)
{
    float * I = fimcl_download(I_gpu);
    printf("Writing to %s\n", filename);
    fim_tiff_write_float(filename, I, NULL, I_gpu->M, I_gpu->N, I_gpu->P);
    free(I);
}

static fimcl_t  * fft_block_of_ones(clu_env_t * clu,
                                    const int64_t M, const int64_t N, const int64_t P,
                                    const int64_t wM, const int64_t wN, const int64_t wP)
{
    /* Return the fft of an image that is 1 in MNP and 0 outside
     * M, N, P is the dimension of the microscopic image
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

    fimcl_t * gOne = fimcl_new(clu, fimcl_real, one, wM, wN, wP);
    here();
    fimcl_sync(gOne);
    here();
    free(one);
    fimcl_t * gfOne = fimcl_fft(gOne);
    fimcl_free(gOne);
    return gfOne;
}

fimcl_t * create_initial_W(clu_env_t * clu,
                           size_t M, size_t N, size_t P,
                           size_t wM, size_t wN, size_t wP,
                           fimcl_t * fft_PSF_gpu)
{

    here();
    fimcl_t * fft_block_ones_gpu = fft_block_of_ones(clu, M, N, P, wM, wN, wP);
    here();
    fimcl_t * W_pre_gpu = fimcl_convolve_conj(fft_block_ones_gpu, fft_PSF_gpu, CLU_KEEP_2ND);
    fft_block_ones_gpu = NULL; // freed already
    here();

    //TODO: Don't download and upload, do it directly on the GPU.
    // took 0.2 s for a job of size 2240x2240x25 (total time with 15 iter 5.1 s)


    float * W_pre = fimcl_download(W_pre_gpu);
    fimcl_free(W_pre_gpu);
    here();

    const size_t wMNP = wM*wN*wP;

    /* Sigma in Bertero's paper, introduced for Eq. 17 */
    float sigma = 0.01;
#pragma omp parallel for shared(W_pre)
    for(size_t kk = 0; kk<wMNP; kk++)
    {
        if(W_pre[kk] > sigma)
        {
            W_pre[kk] = 1.0/W_pre[kk];
        } else {
            W_pre[kk] = 0;
        }
    }

    fimcl_t * W_gpu = fimcl_new(clu, fimcl_real, W_pre, wM, wN, wP);

    free(W_pre);
    return W_gpu;
}

float iter_shb_cl2(fimcl_t ** _xp_gpu, // Output, f_(t+1)
                   fimcl_t * restrict im_gpu, // as above, on GPU
                   fimcl_t * restrict fft_psf_gpu,
                   fimcl_t * restrict pk_gpu, // p_k, estimation of the gradient
                   fimcl_t * W_gpu, // Bertero Weights
                   const dw_opts * s)
{
    fimcl_t * xp_gpu = _xp_gpu[0];
    fimcl_free(xp_gpu);

    fimcl_t * fft_pk_gpu = fimcl_fft(pk_gpu); // Pk -> fft_pk_gpu

    here();

    putdot(s);

    // could be done inplace
    fimcl_t * gy = fimcl_convolve(fft_pk_gpu, fft_psf_gpu, CLU_KEEP_2ND);
    fft_pk_gpu = NULL;

    putdot(s);

    here();

    float error = fimcl_error_idiv(gy, im_gpu);
    putdot(s);

    fimcl_update_y(gy, im_gpu);
    putdot(s);

    fimcl_t * _Y = gy;

#if use_inplace_clfft
    fimcl_fft_inplace(_Y);
    fimcl_t * Y = _Y;
#else
    fimcl_t * Y = fimcl_fft(_Y);
    fimcl_free(_Y);
#endif

    // gx = Y * cK, Y is freed
    // could be done inplace
    fimcl_t * gx = fimcl_convolve_conj(Y, fft_psf_gpu, CLU_KEEP_2ND);
    Y = NULL;
    putdot(s);

    fimcl_real_mul_inplace(pk_gpu, gx); // gx = gpk .* gx
    fimcl_free(pk_gpu);

    if(W_gpu != NULL)
    {
        // gx = gx .* W_gpu
        fimcl_real_mul_inplace(W_gpu, gx);
    }

    // gx[kk] < s->bg ? gx[kk] = s->bg : 0;
    fimcl_positivity(gx, s->bg);

    _xp_gpu[0] = gx;
    return error;
}


float * deconvolve_shb_cl2(float * restrict im,
                           const int64_t M, const int64_t N, const int64_t P,
                           float * restrict psf,
                           const int64_t pM, const int64_t pN, const int64_t pP,
                           dw_opts * s)
{

    if(s->verbosity > 0)
    {
#if use_inplace_clfft
        printf("Deconvolving using shbcl2 (using inplace)\n");
#else
        printf("Deconvolving using shbcl2\n");
#endif
    }

    fprintf(s->log, "Deconvolving with shbcl2\n");
    struct timespec t_deconvolution_start, t_deconvolution_end;
    clock_gettime(CLOCK_REALTIME, &t_deconvolution_start);


    if(s->nIter == 0)
    {
        free(psf);
        return fim_copy(im, M*N*P);
    }

    if(s->bg_auto)
    {
        s->bg = fim_min(im, M*N*P);
        s->bg < 1e-2 ? s->bg = 1e-2 : 0;
        if(s->verbosity > 1)
        {
            printf("Setting the background level to %f\n", s->bg);
        }
    }

    if(fim_maxAtOrigo(psf, pM, pN, pP) == 0)
    {
        if(s->verbosity > 0)
        {
            printf(" ! SHBCL2: PSF is not centered!\n");
        }
        if(s->log != stdout)
        {
            fprintf(s->log, " ! SHBCL2: PSF is not centered!\n");
        }
    }

    /* These are the work dimensions, i.e., dimensions
     * that will be used for all FFTs. Will possibly be increased
     * later on if clFFT handle these sizes.
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
        /* Should we require that pM <= M, etc ? */
        wM = int64_t_max(M, pM);
        wN = int64_t_max(N, pN);
        wP = int64_t_max(P, pP);
    }

    wM = clu_next_fft_size(wM);
    wN = clu_next_fft_size(wN);
    wP = clu_next_fft_size(wP);

    if(s->borderQuality == 0)
    {
        if(wM != M || wN != N || wP != P)
        {
            // This is probably not good.
        }
    }

    /* Total number of pixels */
    const size_t wMNP = wM*wN*wP;

    if(s->verbosity > 0)
    {
        printf("image: [%" PRId64 "x%" PRId64 "x%" PRId64 "], "
               "psf: [%" PRId64 "x%" PRId64 "x%" PRId64 "], "
               "job: [%" PRId64 "x%" PRId64 "x%" PRId64 "]\n",
               M, N, P,
               pM, pN, pP,
               wM, wN, wP);
    }

    /* Set up the OpenCL environment */
    clu_env_t * clu = clu_new(s->verbosity, s->cl_device);
    /* Set up the kernels defined in this function */
    clu_prepare_kernels(clu, M, N, P, wM, wN, wP);

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


    /* Prepare the PSF */
    // cK : "full size" fft of the PSF
    float * Z = fim_malloc(wMNP*sizeof(float));
    memset(Z, 0, wMNP*sizeof(float));
    /* Insert the psf into the bigger Z */
    fim_insert(Z, wM, wN, wP,
               psf, pM, pN, pP);

    free(psf);
    psf = NULL;

    /* Shift the PSF so that the mid is at (0,0,0) */
    int64_t midM, midN, midP = -1;
    float psf_max = 0;
    fim_argmax_max(Z, wM, wN, wP, &midM, &midN, &midP, &psf_max);

    fim_circshift(Z, wM, wN, wP, -midM, -midN, -midP);
    if(Z[0] != psf_max)
    {
        printf("Something went wrong here, the max of the PSF is in the wrong place\n");
        exit(1);
    }

    if(s->fulldump)
    {
        printf("Dumping to fullPSF.tif\n");
        fim_tiff_write_float("fullPSF.tif", Z, NULL, wM, wN, wP);
    }

#if use_inplace_clfft
    fimcl_t * PSF_gpu = fimcl_new(clu, fimcl_real_inplace,
                                  Z, wM, wN, wP);
    fimcl_fft_inplace(PSF_gpu);
    fimcl_t * fft_PSF_gpu = PSF_gpu;
    PSF_gpu = NULL;
    free(Z);
#else
    fimcl_t * PSF_gpu = fimcl_new(clu, fimcl_real,
                                  Z, wM, wN, wP);

    fimcl_t * fft_PSF_gpu = fimcl_fft(PSF_gpu);
    fimcl_free(PSF_gpu);
    free(Z);
#endif

    /* Prepare the image */
    fimcl_t * im_gpu = NULL;

    if(s->psf_pass == 0)
    {
        im_gpu = fimcl_new(clu, fimcl_real, im, M, N, P);
    } else {
        printf("Experimental pre-processing path. Don't use. Work in progress.\n");
        /* Leaves some artifacts at the boundaries. FFT-based
           pre-filtering is a tricky path. */
        float * im_full = fim_malloc(wMNP*sizeof(float));
        memset(im_full, 0, wMNP*sizeof(float));
        /* Insert the psf into the bigger Z */
        fim_insert(im_full, wM, wN, wP,
                   im, M, N, P);
        fimcl_t * im_full_gpu = fimcl_new(clu, fimcl_real, im_full, wM, wN, wP);
        fimcl_t * fft_im_full_gpu = fimcl_fft(im_full_gpu);
        im_full_gpu = NULL;
        // magic thresholding

        fimcl_preprocess(fft_im_full_gpu, fft_PSF_gpu, s->psf_pass);
        im_full_gpu = fimcl_ifft(fft_im_full_gpu);
        float * im_filt = fimcl_download(im_full_gpu);
        float * im_filt_cropped = fim_subregion(im_filt, wM, wN, wP, M, N, P);
        free(im_filt);
        im_gpu = fimcl_new(clu, fimcl_real, im_filt_cropped, M, N, P);
        fimcl_to_tiff(im_gpu, "ifft_fft_input.tif");
        free(im_filt_cropped);
    }

    putdot(s);

    /* Create W_gpu that contains the weights used for the
       boundary handling  */

    fimcl_t * W_gpu = NULL;

    if(s->borderQuality > 0)
    {
        if(s->verbosity > 1)
        {
            printf("Creating weight map for boundary handling\n");
        }
        W_gpu = create_initial_W(clu, M, N, P, wM, wN, wP, fft_PSF_gpu);
    }


    /* Set up the initial guess, x_gpu as well as arrays for the
       previous iterations. p_gpu will be used later on. */

    fimcl_t * x_gpu = NULL;
    fimcl_t * xp_gpu = NULL;
    fimcl_t * p_gpu = NULL;
    {
        float sumg = fim_sum(im, M*N*P);
        float * x = fim_constant(wMNP, sumg / (float) wMNP);

#if use_inplace_clfft
        x_gpu = fimcl_new(clu, fimcl_real_inplace, x, wM, wN, wP);
#else
        x_gpu = fimcl_new(clu, fimcl_real, x, wM, wN, wP);
#endif
        xp_gpu = fimcl_copy(x_gpu);


        free(x);
    }

    if(s->verbosity > 0)
    {
        printf("Iterating "); fflush(stdout);
    }

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
#if use_inplace_clfft
        p_gpu = fimcl_new(clu, fimcl_real_inplace, NULL, wM, wN, wP);
#else
        p_gpu = fimcl_new(clu, fimcl_real, NULL, wM, wN, wP);
#endif


        fimcl_shb_update(p_gpu, x_gpu, xp_gpu, alpha);
        fimcl_positivity(p_gpu, s->bg);

        putdot(s);
        here();
        double err =
            iter_shb_cl2(
                &xp_gpu, // xp is updated to the next guess
                im_gpu,
                fft_PSF_gpu, // FFT of PSF
                p_gpu, // Current guess
                W_gpu, // Weights (to handle boundaries)
                s);

        here();

        dw_iterator_set_error(it, err);
        {
            /* Swap so that the current is named x */
            fimcl_t * t = x_gpu;
            x_gpu = xp_gpu;
            xp_gpu = t;
        }

        putdot(s);
        dw_iterator_show(it, s);
        // benchmark_write(s, it->iter, it->error, x, M, N, P, wM, wN, wP);
        here();
    } /* End of main loop */
    dw_iterator_free(it);
    here();
    fimcl_free(fft_PSF_gpu);
    fft_PSF_gpu = NULL;
    fimcl_free(im_gpu);
    im_gpu = NULL;

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

    float * out_full = fimcl_download(x_gpu);
    fimcl_free(x_gpu);
    if(s->fulldump)
    {
        printf("Dumping to fulldump.tif\n");
        fim_tiff_write("fulldump.tif", out_full, NULL, wM, wN, wP);
    }
    here();

    float * out = fim_subregion(out_full, wM, wN, wP, M, N, P);
    free(out_full);

    here();

    clu_destroy(clu);

    clock_gettime(CLOCK_REALTIME, &t_deconvolution_end);
    float dt_deconvolution = clockdiff(&t_deconvolution_end,
                                       &t_deconvolution_start);
    fprintf(s->log, "Deconvolution took %.2f s\n", dt_deconvolution);

    return out;
}
