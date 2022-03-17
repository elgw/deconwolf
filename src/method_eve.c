#include "method_eve.h"

/* Exponential vector extrapolation (eve) alpha */
 float biggs_alpha_eve(const afloat * restrict Xk,
                       const afloat * restrict Xkm1,
                       const afloat * restrict Ukm1,
                       const afloat * restrict Ukm2,
                       const size_t wMNP)
 {
     /* These calculations take considerable time due to the logf.
      * Full single precision accuracy is probably not needed,
      * so we could think of using approximate calculations
      * or sampling.
      * 0 < \alpha_k < 1
      */

     float numerator = 0;
     float denominator = 0;

#pragma omp parallel for reduction(+:numerator, denominator) shared(Xk, Xkm1, Ukm1, Ukm2)
     for(size_t kk = 0; kk < wMNP; kk++)
     {
         if( (Ukm1[kk] > 0) && (Ukm2[kk] > 0) ) /* Why is U 0 ? */
         {
             float logf_Ukm2_kk = logf(Ukm2[kk]);
             float logf_Ukm1_kk = logf(Ukm1[kk]);
             numerator   += Xk[kk]*logf_Ukm1_kk*Xkm1[kk]*logf_Ukm2_kk;
             denominator += Xkm1[kk]*logf_Ukm2_kk*Xkm1[kk]*logf_Ukm2_kk;
         }
     }

     if(denominator == 0)
     {
         printf("Unexpected denominator=0 in biggs_alpha_eve\n");
         return 0;
     }
     float alpha = numerator / denominator;
     float minalpha = 0.01;
     float maxalpha = 0.98;
     alpha < minalpha ? alpha = minalpha : 0;
     alpha > maxalpha ? alpha = maxalpha : 0;
     return alpha;
 }


/* One RL iteration */
 float iter_eve(
                afloat ** xp, // Output, f_(t+1) xkp1
                const float * restrict im, // Input image
                fftwf_complex * restrict fftPSF,
                afloat * restrict f, // Current guess, xk
                afloat * restrict W, // Bertero Weights
                const int64_t wM, const int64_t wN, const int64_t wP, // expanded size
                const int64_t M, const int64_t N, const int64_t P, // input image size
                __attribute__((unused)) const dw_opts * s)
 {
     const size_t wMNP = wM*wN*wP;

     fftwf_complex * F = fft(f, wM, wN, wP); /* FFT#1 */
     putdot(s);
     afloat * y = fft_convolve_cc_f2(fftPSF, F, wM, wN, wP); /* FFT#2 */
     putdot(s);
     float error = getError(y, im, M, N, P, wM, wN, wP);


     int y_has_zero = 0;
#pragma omp parallel for shared(y, im)
     for(size_t cc = 0; cc < (size_t) wP; cc++){
         for(size_t bb = 0; bb < (size_t) wN; bb++){
             for(size_t aa = 0; aa < (size_t) wM; aa++){
                 size_t yidx = aa + bb*wM + cc*wM*wN;
                 size_t imidx = aa + bb*M + cc*M*N;
                 if(aa < (size_t) M && bb < (size_t) N && cc < (size_t) P)
                 {
                     if(y[yidx] > 0)
                     {
                         y[yidx] = im[imidx]/y[yidx];
                     } else {
                         y_has_zero = 1;
                         y[yidx] = s->bg;
                     }
                 } else {
                     y[yidx]=1e-6;
                 }
             }
         }
     }
     if(y_has_zero == 1)
     {
         if(s->verbosity > 1)
         {
             printf("Zero(s) found in y\n");
         }
     }


     fftwf_complex * F_sn = fft(y, wM, wN, wP); /* FFT#3 */
     fftwf_free(y);
     putdot(s);
     afloat * x = fft_convolve_cc_conj_f2(fftPSF, F_sn, wM, wN, wP); /* FFT#4 */
     putdot(s);

     /* Eq. 18 in Bertero */
     if(W != NULL)
     {
#pragma omp parallel for shared(x,f,W)
         for(size_t cc = 0; cc<wMNP; cc++)
         {
             x[cc] *= f[cc]*W[cc];
         }
     } else {
#pragma omp parallel for shared(x,f)
         for(size_t cc = 0; cc<wMNP; cc++)
         {
             x[cc] *= f[cc];
         }
     }

     xp[0] = x;
     return error;
 }


 float * deconvolve_eve(afloat * restrict im, const int64_t M, const int64_t N, const int64_t P,
                        afloat * restrict psf, const int64_t pM, const int64_t pN, const int64_t pP,
                        dw_opts * s)
 {

     if(s->verbosity > 1)
     {
         printf("Deconvolving\n");
     }

     if(s->nIter == 0)
     {
         free(psf);
         return fim_copy(im, M*N*P);
     }

     const int nIter = s->nIter;

     if(fim_maxAtOrigo(psf, pM, pN, pP) == 0)
     {
         printf("deconv_w: PSF is not centered!\n");
         //return NULL;
     }

     /* This is the 'work dimensions', i.e., dimensions
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

     if(wP % 2 == 1) /* Todo: check if prime instead */
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
     { printf("image: [%" PRId64 "x%" PRId64 "x%" PRId64
              "], psf: [%" PRId64 "x%" PRId64 "x%" PRId64
              "], job: [%" PRId64 "x%" PRId64 "x%" PRId64 "]\n",
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

     /* 1. Expand the PSF to the job size,
      *    transform it and free the original allocation
      */
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
         fim_tiff_write_float("fulldump_PSF.tif", Z, NULL, wM, wN, wP, stdout);
     }

     fftwf_complex * fftPSF = fft(Z, wM, wN, wP);
     fftwf_free(Z);

     putdot(s);


     /* Step 2. Create the Weight map for Bertero boundary handling  */
     afloat * W = NULL;
     if(s->borderQuality > 0)
     {
         /* F_one is 1 over the image domain */
         fftwf_complex * F_one = initial_guess(M, N, P, wM, wN, wP);
         /* Bertero, Eq. 15 */
         W = fft_convolve_cc_conj_f2(fftPSF, F_one, wM, wN, wP);
         F_one = NULL; /* Freed by the call above */

         if(s->borderQuality > 0)
         {
             /* Sigma in Bertero's paper, introduced for Eq. 17 */
             float sigma = 0.01; // Until 2021.11.25 used 0.001
#pragma omp parallel for shared(W)
             for(size_t kk = 0; kk<wMNP; kk++)
             {
                 if(W[kk] > sigma)
                 {
                     W[kk] = 1.0/W[kk];
                 } else {
                     W[kk] = 0;
                 }
             }
         }
     }

     fulldump(s, W, wM, wN, wP, "W");

     putdot(s);

     /* Step 3. */
     float sumg = fim_sum(im, M*N*P);

     float alpha = 0;

     afloat * f = fim_constant(wMNP, sumg/(float) wMNP); /* Initial guess */
     afloat * y = fim_copy(f, wMNP);
     afloat * x1 = fim_copy(f, wMNP);
     afloat * x2 = fim_copy(f, wMNP);

     afloat * x = x1;
     afloat * xp = x2;
     afloat * xm = x2;

     afloat * gm = fim_zeros(wMNP); /* Previous and current quotient */
     afloat * g = fim_zeros(wMNP);
     double alpha_last = 0;

     int it = 0;
     while(it<nIter)
     {

         if(s->iterdump > 0){
             if(it % s->iterdump == 0)
             {
                 afloat * temp = fim_subregion(x, wM, wN, wP, M, N, P);
                 char * outname = gen_iterdump_name(s, it);
                 //fulldump(s, temp, M, N, P, outname);
                 if(s->outFormat == 32)
                 {
                     fim_tiff_write_float(outname, temp, NULL, M, N, P, s->log);
                 } else {
                     fim_tiff_write(outname, temp, NULL, M, N, P, s->log);
                 }
                 free(outname);
                 free(temp);
             }
         }

         struct timespec tstart, tend;
         clock_gettime(CLOCK_REALTIME, &tstart);
         if(it >= 2 && s->biggs > 0)
         {

             clock_gettime(CLOCK_REALTIME, &tstart);


             alpha = biggs_alpha_eve(x, xm, g, gm, wMNP);

             clock_gettime(CLOCK_REALTIME, &tend);
             if(s->showTime){
                 printf("\n--timing alpha: %f\n", clockdiff(&tend, &tstart));
             }
             //printf("alpha = %f\n", alpha);
             alpha = powf(alpha, 1.1);
             //printf("alpha = %f\n", alpha);

             if(alpha > alpha_last)
             {
                 alpha = 0.5*alpha + 0.5*alpha_last;
             }

             alpha_last = alpha;

             clock_gettime(CLOCK_REALTIME, &tstart);
             // printf("Alpha = %f\n", alpha);
             float minvalue = 1e-6;
#pragma omp parallel for shared(xm, x, y)
             for(size_t kk = 0; kk<wMNP; kk++)
             {
                 if(xm[kk] > minvalue)
                 {
                     y[kk] = x[kk]*powf(x[kk]/xm[kk], alpha);
                 }
                 else
                 {
                     y[kk] = minvalue;
                 }
             }
             /* Possibly: Set back to pre-extrapolated value (i.e. x) if
              *  the likelihood didn't increase
              */
         } else {
             /* Biggs accelaration disabled */
#pragma omp parallel for shared(xm, y, x)
             for(size_t kk = 0; kk<wMNP; kk++)
             {
                 if(xm[kk] > 0)
                 {
                     y[kk] = x[kk];
                 }
             }
         }
         clock_gettime(CLOCK_REALTIME, &tend);
         if(s->showTime){
             printf("\n--timing powf (alpha = %f): %f\n", alpha, clockdiff(&tend, &tstart));
         }
         putdot(s);

         xp = xm;
         fftwf_free(xp);
         double err = iter_eve(
                               &xp, // xp is updated to the next guess
                               im,
                               fftPSF,
                               y, // Current guess
                               W, // Weights (to handle boundaries)
                               wM, wN, wP, // Expanded size
                               M, N, P, // Original size
                               s);

         putdot(s);

         if(s->verbosity > 0){
             printf("\r                                             ");
             printf("\rIteration %3d/%3d, fMSE=%e ", it+1, nIter, err);
             fflush(stdout);
         }

         if(s->log != NULL)
         {
             fprintf(s->log, "Iteration %d/%d, fMSE=%e\n", it+1, nIter, err);
             fflush(s->log);
         }

         {
             afloat * swap = g;
             g = gm; gm = swap;
         }

         /* The "correction factor" u_k in Biggs */
         fim_div(g, xp, y, wMNP); // g = xp/y

         /* Any zero pixels have their correction factor set to unity */
#pragma omp parallel for shared(y,g)
             for(size_t kk = 0; kk<wMNP; kk++)
             {
                 if(y[kk] <= s->bg)
                 {
                     g[kk] = 1;
                 }
             }

         if(W != NULL)
         {
#pragma omp parallel for shared(W, g)
             for(size_t kk = 0; kk<wMNP; kk++)
             {

                 if(W[kk] == 0)
                 {
                     g[kk] = 0;
                 }
             }
         }

         if(s->bg > 0)
         {
#pragma omp parallel for shared(xp)
             for(size_t kk = 0; kk<wMNP; kk++)
             {
                 if(xp[kk] < s->bg)
                 {
                     xp[kk] = s->bg;
                 }
             }
         }

         xm = x;
         x = xp;
         xp = NULL;
         benchmark_write(s, it, err, x, M, N, P, wM, wN, wP);
         it++;
     } /* End of main loop */

     if(s->verbosity > 0) {
         printf("\n");
     }

     if(W != NULL)
     {
         fftwf_free(W);
     }

     fulldump(s, x, wM, wN, wP, "fulldump_x.tif");


     fftwf_free(f);

     if(xm != NULL)
     {
         fftwf_free(xm);
         xm = NULL;
     }

     fftwf_free(g);
     fftwf_free(gm);
     fftwf_free(fftPSF);
     fftwf_free(y);

     /* Extract the observed region from the last iteration */
     afloat * out = fim_subregion(x, wM, wN, wP, M, N, P);

     if(x != NULL)
     {
         fftwf_free(x); x = NULL;
     }
     return out;
 }
