#include "method_ave.h"


float alpha_ave(const afloat * restrict g,
                  const afloat * restrict gm,
                  const size_t wMNP, int mode)
{
    /* Biggs, eq. 18 */
    double ggm = 0;
    double gmgm = 0;

#pragma omp parallel for reduction(+: ggm, gmgm) shared(g, gm)
    for(size_t kk = 0; kk<wMNP; kk++)
    {
        ggm += g[kk]*gm[kk];
        gmgm += pow(gm[kk], 2);
    }

    float alpha = 0;
    switch(mode)
    {
    case 1:
        /* Ok for everything seen so far ... */
        alpha = pow(ggm/(gmgm+1e-7), 2);
        break;
    case 2:
        /* Ok for most images, problematic for high contrast regions */
        alpha = ggm/(gmgm+1e-7);
        break;
    case 3:
        /* Used in the paper. Clearly too aggressive. */
        alpha = sqrt(ggm/(gmgm+1e-7));
        break;
    }

    alpha <= 0 ? alpha = 1e-5 : 0;
    alpha >= 1 ? alpha = 1.0-1e-5 : 0;

    return alpha;
}


static float iter_ave(
            afloat ** xp, // Output, f_(t+1)
            const float * restrict im, // Input image
            fftwf_complex * restrict cK, // fft(psf)
            afloat * restrict f, // Current guess
            afloat * restrict W, // Bertero Weights
            const int64_t wM, const int64_t wN, const int64_t wP, // expanded size
            const int64_t M, const int64_t N, const int64_t P, // input image size
            __attribute__((unused)) const dw_opts * s)
 {
     // We could reduce memory even further by using
     // the allocation for xp
     const size_t wMNP = wM*wN*wP;

     fftwf_complex * F = fft(f, wM, wN, wP);
     putdot(s);
     afloat * y = fft_convolve_cc_f2(cK, F, wM, wN, wP); // F is freed
     float error = getError(y, im, M, N, P, wM, wN, wP, s->metric);
     putdot(s);

     const double mindiv = 1e-6; // Before 2021.11.24: 1e-6
#pragma omp parallel for shared(y)
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

     fftwf_complex * F_sn = fft(y, wM, wN, wP);
     fftwf_free(y);

     afloat * x = fft_convolve_cc_conj_f2(cK, F_sn, wM, wN, wP);

     /* Eq. 18 in Bertero */
     if(W != NULL)
     {
#pragma omp parallel for shared(x, f, W)
         for(size_t cc = 0; cc<wMNP; cc++)
         {
             x[cc] *= f[cc]*W[cc];
         }
     } else {
#pragma omp parallel for shared(x, f)
         for(size_t cc = 0; cc<wMNP; cc++)
         {
             x[cc] *= f[cc];
         }
     }

     xp[0] = x;
     return error;
 }


float * deconvolve_ave(afloat * restrict im,
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
         free(psf);
         return fim_copy(im, M*N*P);
     }

     if(fim_maxAtOrigo(psf, pM, pN, pP) == 0)
     {
         printf(" ! AVE PSF is not centered!\n");
         fprintf(s->log, " ! AVE PSF is not centered!\n");
         //return NULL;
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

     fftwf_free(psf);

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

     //printf("initial guess\n"); fflush(stdout);
     afloat * W = NULL;
     if(s->borderQuality > 0)
     {
         fftwf_complex * F_one = initial_guess(M, N, P, wM, wN, wP);
         afloat * P1 = fft_convolve_cc_conj_f2(cK, F_one, wM, wN, wP); // can't replace this one with cK!
         //printf("P1\n");
         //fim_stats(P1, pM*pN*pP);
         //  writetif("P1.tif", P1, wM, wN, wP);

         putdot(s);

         /* Sigma in Bertero's paper, introduced for Eq. 17 */
         if(s->borderQuality > 0)
         {
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
         }
         W = P1;
     }
     //writetif("W.tif", W, wM, wN, wP);

     // Original image -- expanded
     //  float * G = fftwf_malloc(wMNP*sizeof(float));
     // memset(G, 0, wMNP*sizeof(float));
     // fim_insert(G, wM, wN, wP, im, M, N, P);
     //  writetif("G.tif", G, wM, wN, wP);


     float sumg = fim_sum(im, M*N*P);

     float alpha = 0;
     afloat * f = fim_constant(wMNP, sumg/(float) wMNP);
     afloat * y = fim_copy(f, wMNP);

     afloat * x1 = fim_copy(f, wMNP);
     afloat * x2 = fim_copy(f, wMNP);
     afloat * x = x1;
     afloat * xp = x2;
     afloat * xm = x2; // fim_copy(f, wMNP);

     afloat * gm = fim_zeros(wMNP);
     afloat * g = fim_zeros(wMNP);

     dw_iterator_t * it = dw_iterator_new(s);
     while(dw_iterator_next(it) >= 0)
     {

         if(s->iterdump > 0){
             if(it->iter % s->iterdump == 0)
             {
                 afloat * temp = fim_subregion(x, wM, wN, wP, M, N, P);
                 char * outname = gen_iterdump_name(s, it->iter);
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


         if(it == 0)
         {
             memcpy(y, x, wMNP*sizeof(float));
         } else {
#pragma omp parallel for shared(y,x,xm)
         for(size_t kk = 0; kk<wMNP; kk++)
         {
             y[kk] = x[kk] + alpha*(x[kk]-xm[kk]);
         }

         /* Enforce a priori information about the lowest possible value */
         if(s->positivity)
         {
#pragma omp parallel for shared(y)
             for(size_t kk = 0; kk<wMNP; kk++)
             {
                 if(y[kk] < s->bg)
                 {
                     y[kk] = s->bg;
                 }
             }
         }
         }
         if(s->psigma > 0)
         {
             printf("gsmoothing()\n");
             fim_gsmooth(y, wM, wN, wP, s->psigma);
         }

         putdot(s);

         xp = xm;
         fftwf_free(xp);
         double err = iter_ave(
                           &xp, // xp is updated to the next guess
                           im,
                           cK, // FFT of PSF
                           y, // Current guess
                           W, // Weights (to handle boundaries)
                           wM, wN, wP, // Expanded size
                           M, N, P, // Original size
                           s);

         putdot(s);
         dw_iterator_set_error(it, err);
         dw_iterator_show(it, s);


         afloat * swap = g;
         g = gm; gm = swap;

         fim_minus(g, xp, y, wMNP); // g = xp - y

         if(it->iter > 0) { /* 20211127 'it > 0' */
             alpha = alpha_ave(g, gm, wMNP, s->biggs);
         }

         xm = x;
         x = xp;
         xp = NULL;


         benchmark_write(s, it->iter, err, x, M, N, P, wM, wN, wP);

     } // End of main loop


     if(s->verbosity > 0) {
         printf("\n");
     }

     if(W != NULL)
     {
         fftwf_free(W); // is P1
     }

     if(s->fulldump)
     {
         printf("Dumping to fulldump.tif\n");
         fim_tiff_write("fulldump.tif", x, NULL, wM, wN, wP);
     }

     afloat * out = fim_subregion(x, wM, wN, wP, M, N, P);


     //  printf("DEBUG: writing final_full_tif\n");
     //  fim_tiff_write("final_full.tif", x, wM, wN, wP);
     fftwf_free(f);
     if(x != NULL)
     {
         fftwf_free(x); x = NULL;
     }
     if(xm != NULL)
     {
         fftwf_free(xm);
         xm = NULL;
     }
     fftwf_free(g);
     fftwf_free(gm);
     fftwf_free(cK);
     //  if(cKr != NULL)
     //  { fftwf_free(cKr); }
     fftwf_free(y);
     return out;
 }
