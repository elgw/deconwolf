#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <fftw3.h>

#include "tiffio.h"
#include "fft.h"


float iter(float * xp, float * g, 
    fftwf_complex * K,
    fftwf_complex * Kr,
    float * f, 
    fftwf_complex * W,
    fftwf_complex * G)
{
  fftwf_complex * F = fft(f, wM, wN, wP);
  float * y = conv_FF(K, F, wMNP);
  fftwf_free(F);
  float error = getError(y, g, M, N, P, wM, wN, wP);
  sn = float_div(G, y, wMNP); // G isn't fft(g) ? 
  fftwf_free(y);
  fftwf_complex * F_sn = fft(sn, wM, wN, wP);
  float * x = conv_FF(Kr, F_sn);
  for(size_t kk = 0; kk<wMNP; kk++)
  {
    xp[kk] = f[kk]*x[kk]*W[kk];
  }
  return error;
}

fftwf_complex * initial_guess(int M, int N, int P, 
    int wM, int wN, int wP)
{
/* Create initial guess: the fft of an image that is 1 in MNP and 0 outside
 * M, N, P is the dimension of the microscopic image
 */
  assert(wM > M); assert(wN>N); assert(wP>P);

  float * one = fftwf_malloc(wM*wN*wP*sizeof(float));
  for(size_t kk = 0; kk<wM*wN*wP; kk++)
    one[kk] = 0;
  for(int aa = 0; aa < M; aa++) {
    for(int bb = 0; bb < N; bb++) {
      for(int cc = 0; cc < P; cc++) {
        one[aa + wM*bb + wM*wN*cc] = 1;
      }
    }
  }

  fftwf_complex * Fone = fftwf_malloc(wM*wN*wP*sizeof(fftwf_complex));
   fftwf_plan p = fftwf_plan_dft_r2c_3d(wP, wM, wN, 
      one, // Float
      Fone, // fftwf_complex 
      FFTW_ESTIMATE);

  fftwf_execute(p); 
  fftwf_destroy_plan(p);

  fftwf_free(one);
  return Fone;
}

float * expandIm_a(float * in, 
    int pM, int pN, int pP, 
    int M, int N, int P)
{
  assert(pM<=M);
  assert(pN<=N);
  assert(pP<=P);

  float * out = fftwf_malloc(M*N*P*sizeof(float));
  assert(in != NULL);
  assert(out != NULL);
  for(size_t kk = 0; kk<M*N*P; kk++)
    out[kk] = 0;

  for(int mm = 0; mm<pM; mm++)
  {
    for(int nn = 0; nn<pN; nn++)
    {
      for(int pp = 0; pp<pP; pp++)
      {
        float x = in[mm + nn*pM + pp*pM*pN];
        //        printf("%f\n", x);
        out[mm + nn*M + pp*M*N] = x;
      }
    }
  }
  return out;
}

void usage(int argc, char ** argv)
{
  printf("Usage:\n");
  printf("$%s image.tif psf.tif\n", argv[0]);
}

float * deconvolve(float * im, int M, int N, int P,
    float * psf, int pM, int pN, int pP)
{
  int niter = 5;
  /* Deconvolve im using psf */

  // Determine optimal size to work with ...
  // Optimally pP_opt = 2*P+1
  // If nP > pP_opt it should be cropped.
  int pP_opt = 2*P+1;

  float * psf1 = NULL;
  if(pP > 2*P+1)
  {
    printf("Cropping PSF in Z, from %d slices to %d\n", pP, pP_opt);

    psf1 = malloc(pM*pN*pP_opt*sizeof(float));
    assert(psf1 != NULL);
    size_t skip = pM*pN*(pP-pP_opt)/2;
    memcpy(psf1+skip, psf, pM*pN*pP_opt*sizeof(float));
    pP = P;
    free(psf);
  } else {
    psf1 = psf;
  }

  /* PAD etc ... and then continue in new function with less cluttered namespace */

  /*
  assert(psf1 != NULL);
  printf("Expanding\n"); fflush(stdout);
  float * psf2 = expandIm_a(psf1, pM, pN, pP, M, N, P);
  free(psf1);*/

  // This is the dimensions that we have to work with
  int wM = M + pM -1;
  int wN = N + pN -1;
  int wP = P + pP -1;
  size_t wMNP = wM*wN*wP;

  fftwf_complex * F_one = initial_guess(M, N, P, wM, wN, wP);
  // float * P1 = convolve_FF(Kr, Fone, wM, wN, wP);
  // W = invert(P1, wM*wN*wP);
  
  float sigma = 0.001;
  // fftwf_free(P1);
  // float * G = ... the original image expanded
  
  float alpha = 0;

  //float * x = aligned_copy(f, wMNP);
  //float * y = aligned_copy(f, wMNP);
  //float * xp = aligned_copy(f, wMNP);
  //float * gm = aligned_zeros(wMNP);
  //float * g = aligned_zeros(wMNP);
  
  while(it<niter)
  {
    // for(size_t idx = 0; idx<wMNP; idx++)
    // { y[kk] = x[kk] + alpha*(x[kk]-xm[kk]);
    // y[kk] < 0 ? y[kk] = 0 : 0;
    // }
    // double err = iter(xp, g1, K, Kr, y, W, G);
    // printf("Iteration %d, error=%e\n", err);
    // memcpy(gm, g, wMNP*sizeof(float));
    // for(size_t idx = 0; idx<wMNP; idx++)
    // { g[kk] = xp[kk]-y[kk]; }
    // if(it > 0) {
    // alpha = update_alpha(g, gm, wMNP);
    // memcpy(xm, x, wMNP*sizeof(float));
    it = it + 1;
  }

  return(im);
}

int main(int argc, char ** argv)
{

  if(argc < 3)
  {
    usage(argc, argv);
    exit(1);
  }
  char * imFile = argv[1];
  char * psfFile = argv[2];
  char * outFile = malloc(100*sizeof(char));
  sprintf(outFile, "dummy.tif");
  //  sprintf(outFile, "dw_%s", imFile);

  printf("Settings:\n");
  printf("\t image:  %s\n", imFile);
  printf("\t psf:    %s\n", psfFile);
  printf("\t output: %s\n", outFile);
  printf("Press enter to continue\n");
  getchar();

  int M = 0, N = 0, P = 0;
  float * im = readtif_asFloat(imFile, &M, &N, &P);

  int pM = 0, pN = 0, pP = 0;
  float * psf = readtif_asFloat(psfFile, &pM, &pN, &pP);
  assert(psf != NULL);

  float * out = deconvolve(im, M, N, P, psf, pM, pN, pP);

  int exitstatus = 0;
  if(out == NULL)
  {
    printf("Nothing to write to disk :(\n");
  } else 
  {
    printf("Writing\n"); fflush(stdout);
    floatimage_normalize(out, M*N*P);
    writetif(outFile, out, M, N, P);
    exitstatus = 1;
  }

  printf("Closing down the office\n"); fflush(stdout);
  fftwf_free(im);
  fftwf_free(psf);
  free(outFile);
  return(exitstatus);
}


// Graveyard

/*  float * im = fftwf_malloc(M*N*P*sizeof(float));
    memcpy(im0, im, M*N*P*sizeof(float));
    free(im0); */

