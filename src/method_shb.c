#include "method_shb.h"


#include "fim.h"
#include "fft.h"

//#define here() printf("\n%s %s %d\n", __FILE__, __FUNCTION__, __LINE__);
#define here() ;

static float iter_shb(dw_fft * ff,
    float ** xp, // Output, f_(t+1)
    const float * restrict im, // Input image
    fftwf_complex * restrict cK, // fft(psf)
    float * restrict pk, // p_k, estimation of the gradient
    float * restrict W, // Bertero Weights
    const int64_t wM, const int64_t wN, const int64_t wP, // expanded size
    const int64_t M, const int64_t N, const int64_t P, // input image size
    __attribute__((unused)) const dw_opts * s)
{
    const size_t wMNP = wM*wN*wP;

    fftwf_complex * Pk = fft(ff, pk, wM, wN, wP);

    putdot(s);
    float * y = fft_convolve_cc_f2(ff, cK, Pk, wM, wN, wP); // Pk is freed

    float error = getError(y, im, M, N, P, wM, wN, wP, s->metric);
    putdot(s);


    const double mindiv = 1e-6; /* Smallest allowed divisor */
#pragma omp parallel for shared(y, im) collapse(2)
    for(size_t cc =0; cc < (size_t) wP; cc++){
        for(size_t bb = 0; bb < (size_t) wN; bb++){
            for(size_t aa = 0; aa < (size_t) wM; aa++){
                size_t yidx = aa + bb*wM + cc*wM*wN;
                size_t imidx = aa + bb*M + cc*M*N;
                if(aa < (size_t) M && bb < (size_t) N && cc < (size_t) P)
                {
                    /* abs and sign */
                    fabs(y[yidx]) < mindiv ? y[yidx] = copysign(mindiv, y[yidx]) : 0;
                    y[yidx] = im[imidx]/y[yidx];
                } else {
                    y[yidx]=0;
                }
            }
        }
    }

    here();
    fftwf_complex * Y = fft_and_free(ff, y, wM, wN, wP);
    here();
    float * x = fft_convolve_cc_conj_f2(ff, cK, Y, wM, wN, wP); // Y is freed
    here();
    /* Eq. 18 in Bertero */
    if(W != NULL)
    {
#pragma omp parallel for shared(x, pk, W)
        for(size_t cc = 0; cc<wMNP; cc++)
        {
            x[cc] *= pk[cc]*W[cc];
        }
    } else {
#pragma omp parallel for shared(x, pk)
        for(size_t cc = 0; cc<wMNP; cc++)
        {
            x[cc] *= pk[cc];
        }
    }
    fim_free(pk);
    here();
    xp[0] = x;
    return error;
}


float * deconvolve_shb(float * restrict im,
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

    ftif_t * ftif = fim_tiff_new(s->log, s->verbosity);

    if(s->verbosity > 1)
    {
        printf("Deconvolving\n");
    }


    if(s->nIter == 0)
    {
        fim_free(psf);
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
        fprintf(s->log, "Setting the background level to %f\n", s->bg);
    }


    if(fim_maxAtOrigo(psf, pM, pN, pP) == 0)
    {
        if(s->verbosity > 0)
        {
            printf(" ! SHB: PSF is not centered!\n");
        }
        if(s->log != stdout)
        {
            fprintf(s->log, " ! SHB: PSF is not centered!\n");
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

    /* Total number of pixels */
    size_t wMNP = wM*wN*wP;

    if(s->verbosity > 0)
    {
        printf("image: [%" PRId64 "x%" PRId64 "x%" PRId64 "], "
               "psf: [%" PRId64 "x%" PRId64 "x%" PRId64 "], "
               "job: [%" PRId64 "x%" PRId64 "x%" PRId64 "]\n",
             M, N, P, pM, pN, pP, wM, wN, wP);
    }
    if(s->verbosity > 1)
    {
        printf("Estimated peak memory usage: %.1f GB\n", wMNP*35.0/1e9);
    }
    fprintf(s->log, "image: [%" PRId64 "x%" PRId64 "x%" PRId64 "]\n"
            "psf: [%" PRId64 "x%" PRId64 "x%" PRId64 "]\n"
            "job: [%" PRId64 "x%" PRId64 "x%" PRId64 "] (%zu voxels)\n",
            M, N, P, pM, pN, pP, wM, wN, wP, wMNP);
    fflush(s->log);

    //myfftw_stop(); nope that was not the problem
    //myfftw_start(s->nThreads_FFT, s->verbosity, s->log);
    dw_fft * ff = dw_fft_new(s->nThreads_FFT, s->verbosity, s->log, wM, wN, wP, s->fftw3_planning);

    if(s->verbosity > 0)
    {
        printf("Iterating "); fflush(stdout);
    }

    // cK : "full size" fft of the PSF
    float * Z = fim_malloc(wMNP*sizeof(float));
    memset(Z, 0, wMNP*sizeof(float));
    /* Insert the psf into the bigger Z */
    fim_insert(Z, wM, wN, wP,
               psf, pM, pN, pP);

    fim_free(psf);

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
        fim_tiff_write_float(ftif, "fullPSF.tif", Z, NULL, wM, wN, wP);
    }

    here();
    fftwf_complex * cK = fft_and_free(ff, Z, wM, wN, wP);
    here();
    putdot(s);

    float * W = NULL;
    /* Sigma in Bertero's paper, introduced for Eq. 17 */
    if(s->borderQuality > 0)
    {
        fftwf_complex * F_one = initial_guess(ff, M, N, P, wM, wN, wP);
        here();
        float * P1 = fft_convolve_cc_conj_f2(ff, cK, F_one, wM, wN, wP);
        here();
        float sigma = 0.01; // Until 2021.11.25 used 0.001
#pragma omp parallel for shared(P1)
        for(size_t kk = 0; kk<wMNP; kk++)
        {
            if(P1[kk] > sigma)
            {
                P1[kk] = 1.0/P1[kk];
            } else {
                P1[kk] = 0;
            }
        }
        W = P1;
        if(0){
        printf("\n\nTesting\n\n");
        for(size_t kk = 0; kk<wMNP; kk++)
        {
            P1[kk] = 0;
        }
        for(int mm = 0; mm<M; mm++)
        {
            for(int nn = 0; nn<N; nn++)
            {
                for(int pp = 0; pp<P; pp++)
                {
                    P1[mm + wM*nn + wM*wN*pp] = 1;
                }
            }
        }
        }
    }

    float sumg = fim_sum(im, M*N*P);

    /*  x -- the initial guess
     *  xp -- the previous guess, initially set the same as x
     */

    float * x = NULL;
    float * xp = NULL;

    if(s->start_condition == DW_START_FLAT)
    {
         x = fim_constant(wMNP, sumg/wMNP);
         xp = fim_copy(x, wMNP);
    }

    if(s->start_condition == DW_START_IDENTITY)
    {
        x = fim_malloc(wMNP*sizeof(float));
        fim_insert(x, wM, wN, wP,
                   im, M, N, P);
        xp = fim_copy(x, wMNP);
    }

    if(s->start_condition == DW_START_LP)
    {
        x = fim_malloc(wMNP*sizeof(float));
        float * im_lp = fim_copy(im, M*N*P);
        fim_gsmooth(im_lp, M, N, P, 8);
        fim_insert(x, wM, wN, wP,
                   im_lp, M, N, P);
        xp = fim_copy(x, wMNP);
        fim_free(im_lp);
    }


    if( x == NULL || xp == NULL)
    {
        fprintf(stderr, "Internal error in %s %d\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }



    dw_iterator_t * it = dw_iterator_new(s);
    while(dw_iterator_next(it) >= 0)
    {
        here();
        if(s->iterdump > 0){
            if(it->iter % s->iterdump == 0)
            {
                float * temp = fim_subregion(x, wM, wN, wP, M, N, P);
                if(s->offset > 0)
                {
                    fim_add_scalar(temp, M*N*P, -s->offset);
                    fim_project_positive(im, M*N*P);
                }
                char * outname = gen_iterdump_name(s, it->iter);
                //printf(" Writing current guess to %s\n", outname);
                if(s->outFormat == 32)
                {
                    fim_tiff_write_float(ftif, outname, temp, NULL, M, N, P);
                } else {
                    fim_tiff_write(ftif, outname, temp, NULL, M, N, P);
                }
                free(outname);
                fim_free(temp);
            }
        }

        //float * p = fim_copy(x, wMNP);
        float * p = xp; /* We don't need xp more */

        /* Eq. 10 in SHB paper */
        double alpha = ((float) it->iter-1.0)/((float) it->iter+2.0);
        alpha < 0 ? alpha = 0: 0;
        alpha > s->alphamax ? alpha = s->alphamax : 0;


#pragma omp parallel for shared(p, x, xp)
            /* To be interpreted as p^k in Eq. 7 of SHB */
            for(size_t kk = 0; kk<wMNP; kk++)
            {
                p[kk] = x[kk] + alpha*(x[kk]-xp[kk]);
                p[kk] < s->bg ? p[kk] = s->bg : 0; // TODO
            }


        putdot(s);

        double err = iter_shb(ff,
            &xp, // xp is updated to the next guess
            im,
            cK, // FFT of PSF
            p, // Current guess
            //p,
            W, // Weights (to handle boundaries)
            wM, wN, wP, // Expanded size
            M, N, P, // Original size
            s);
        here();
	//        free(p); // free'ed in iter_shb
        here();
        dw_iterator_set_error(it, err);
        if(1){
            /* Swap so that the current is named x */
            float * t = x;
            x = xp;
            xp = t;
        }
        //free(p);
        /* Enforce a priori information about the lowest possible value */
        if(s->positivity)
        {
#pragma omp parallel for shared(x)
            for(size_t kk = 0; kk<wMNP; kk++)
            {
                if(x[kk] < s->bg)
                {
                    x[kk] = s->bg;
                }
            }
        }
        here();
        putdot(s);
        dw_iterator_show(it, s);
        benchmark_write(s, it->iter, it->error, x, M, N, P, wM, wN, wP);

    } /* End of main loop */
    dw_iterator_free(it);

    {
        /* Swap back so that x is the final iteration */
        float * t = x;
        x = xp;
        xp = t;
    }

    if(xp != NULL)
    {
        fim_free(xp);
    }

    if(s->verbosity > 0) {
        printf("\n");
    }

    if(W != NULL)
    {
        fim_free(W); /* Allocated as P1 */
    }
    fim_free(cK);

    if(s->fulldump)
    {
        printf("Dumping to fulldump.tif\n");
        fim_tiff_write(ftif, "fulldump.tif", x, NULL, wM, wN, wP);
    }

    float * out = fim_subregion(x, wM, wN, wP, M, N, P);

    if(x != NULL)
    {
        fim_free(x);
    }

    dw_fft_destroy(ff);
    ff = NULL;
    fim_tiff_destroy(ftif);
    return out;
}
