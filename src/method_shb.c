#include "method_shb.h"

float * deconvolve_shb(afloat * restrict im,
                       const int64_t M, const int64_t N, const int64_t P,
                       afloat * restrict psf,
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

    const int nIter = s->nIter;


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

    /* Total number of pixels */
    size_t wMNP = wM*wN*wP;

    if(s->verbosity > 0)
    { printf("image: [%" PRId64 "x%" PRId64 "x%" PRId64 "], psf: [%" PRId64 "x%" PRId64 "x%" PRId64 "], job: [%" PRId64 "x%" PRId64 "x%" PRId64 "]\n",
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

    fft_train(wM, wN, wP,
              s->verbosity, s->nThreads,
              s->log);

    if(s->verbosity > 0)
    {
        printf("Iterating "); fflush(stdout);
    }

    // cK : "full size" fft of the PSF
   afloat * Z = fftwf_malloc(wMNP*sizeof(float));
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

    fftwf_complex * cK = fft(Z, wM, wN, wP);
    //fim_tiff_write("Z.tif", Z, wM, wN, wP);
    fftwf_free(Z);

    putdot(s);

    /* <-- This isn't needed ...
       fftwf_complex * cKr = NULL;
       float * Zr = fftwf_malloc(wMNP*sizeof(float));
       memset(Zr, 0, wMNP*sizeof(float));
       float * psf_flipped = malloc(wMNP*sizeof(float));
       memset(psf_flipped, 0, sizeof(float)*wMNP);
       fim_flipall(psf_flipped, psf, pM, pN, pP);
       fim_insert(Zr, wM, wN, wP, psf_flipped, pM, pN, pP);
       free(psf_flipped);
       fim_circshift(Zr, wM, wN, wP, -(pM-1)/2, -(pN-1)/2, -(pP-1)/2);
       cKr = fft(Zr, wM, wN, wP);
       // Not needed since f(-x) = ifft(conf(fft(f)))
       fftwf_free(Zr);
       --> */



    //printf("P1\n");
    //fim_stats(P1, pM*pN*pP);
    //  writetif("P1.tif", P1, wM, wN, wP);

    putdot(s);

    afloat * W = NULL;
    /* Sigma in Bertero's paper, introduced for Eq. 17 */
    if(s->borderQuality > 0)
    {
        fftwf_complex * F_one = initial_guess(M, N, P, wM, wN, wP);
        afloat * P1 = fft_convolve_cc_conj_f2(cK, F_one, wM, wN, wP); // can't replace this one with cK!
        float sigma = 0.01; // Until 2021.11.25 used 0.001
#pragma omp parallel for shared(P1)
        for(size_t kk = 0; kk<wMNP; kk++)
        {
            if(P1[kk] > sigma)
            {
                P1[kk] = 1/P1[kk];
            } else {
                P1[kk] = 0;
            }
        }
        W = P1;
    }



    float sumg = fim_sum(im, M*N*P);

    float alpha = 0;

    /* x is the initial guess, initially the previous iteration,
     *  xp is
     *  set to be the same */

    afloat * x = fim_constant(wMNP, sumg/wMNP);

    if(0)
    {
        printf("EXPERIMENTAL\n");
        /* If only few iterations can be afforded it makes sense to
         * start with the observed images as the starting guess */
        fim_insert(x, wM, wN, wP, im, M, N, P);
    }
    afloat * xp = fim_copy(x, wMNP);

    int it = 0;
    while(it<nIter)
    {

        if(s->iterdump > 0){
            if(it % s->iterdump == 0)
            {
                afloat * temp = fim_subregion(x, wM, wN, wP, M, N, P);
                char * outname = gen_iterdump_name(s, it);
                //printf(" Writing current guess to %s\n", outname);
                if(s->outFormat == 32)
                {
                    fim_tiff_write_float(outname, temp, NULL, M, N, P);
                } else {
                    fim_tiff_write(outname, temp, NULL, M, N, P);
                }
                free(outname);
                free(temp);
            }
        }

        //afloat * p = fim_copy(x, wMNP);
        afloat * p = xp; /* We don't need xp more */

            /* Eq. 10 in SHB paper */
            alpha = ((float) it-1.0)/((float) it+2.0);
            //alpha /= 1.5;
            alpha < 0 ? alpha = 0: 0;



#pragma omp parallel for shared(p, x, xp)
            /* To be interpreted as p^k in Eq. 7 of SHB */
            for(size_t kk = 0; kk<wMNP; kk++)
            {
                p[kk] = x[kk] + alpha*(x[kk]-xp[kk]);
            }



        if(s->psigma > 0)
        {
            printf("gsmoothing()\n");
            fim_gsmooth(x, wM, wN, wP, s->psigma);
        }

        putdot(s);

        double err = iter_shb(
            &xp, // xp is updated to the next guess
            im,
            cK, // FFT of PSF
            p, // Current guess
            //p,
            W, // Weights (to handle boundaries)
            wM, wN, wP, // Expanded size
            M, N, P, // Original size
            s);
        fftwf_free(p); // i.e. p

        if(1){
            /* Swap so that the current is named x */
            afloat * t = x;
            x = xp;
            xp = t;
        }
        //fftwf_free(p);
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


        putdot(s);
        dw_show_iter(s, it, nIter, err);

        benchmark_write(s, it, err, x, M, N, P, wM, wN, wP);
        it++;
    } // End of main loop

    {
        /* Swap back so that x is the final iteration */
        afloat * t = x;
        x = xp;
        xp = t;
    }

    if(xp != NULL)
    {
        fftwf_free(xp);
    }

    if(s->verbosity > 0) {
        printf("\n");
    }

    if(W != NULL)
    {
        fftwf_free(W); /* Allocated as P1 */
    }
    fftwf_free(cK);

    if(s->fulldump)
    {
        printf("Dumping to fulldump.tif\n");
        fim_tiff_write("fulldump.tif", x, NULL, wM, wN, wP);
    }

    afloat * out = fim_subregion(x, wM, wN, wP, M, N, P);

    if(x != NULL)
    {
        fftwf_free(x);
    }



    return out;
}


float iter_shb(
    afloat ** xp, // Output, f_(t+1)
    const float * restrict im, // Input image
    fftwf_complex * restrict cK, // fft(psf)
    afloat * restrict pk, // p_k, estimation of the gradient
    afloat * restrict W, // Bertero Weights
    const int64_t wM, const int64_t wN, const int64_t wP, // expanded size
    const int64_t M, const int64_t N, const int64_t P, // input image size
    __attribute__((unused)) const dw_opts * s)
{
    // We could reduce memory even further by using
    // the allocation for xp
    const size_t wMNP = wM*wN*wP;

    fftwf_complex * Pk = fft(pk, wM, wN, wP);
    putdot(s);
    afloat * y = fft_convolve_cc_f2(cK, Pk, wM, wN, wP); // Pk is freed
    float error = getError(y, im, M, N, P, wM, wN, wP, s->metric);
    putdot(s);

    const double mindiv = 1e-6; /* Smallest allowed divisor */
#pragma omp parallel for shared(y, im)
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

    fftwf_complex * Y = fft(y, wM, wN, wP);
    fftwf_free(y);

    afloat * x = fft_convolve_cc_conj_f2(cK, Y, wM, wN, wP); // Y is freed

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

    xp[0] = x;
    return error;
}
