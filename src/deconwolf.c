#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <fftw3.h>

#include "tiffio.h"
#include "fft.h"


int psfIsCentered(float * V, int M, int N, int P)
{
  // Check that the PSF looks reasonable
  if( M % 2 == 0)
  { return 0; }

  if( N % 2 == 0)
  { return 0; }

  if (P % 2 == 0)
  { return 0; }

  int mM = (M-1)/2;
  int mN = (N-1)/2;
  int mP = (P-1)/2;
  float maxV = 0;
  for(size_t kk = 0; kk< M*N*P; kk++)
  {
    V[kk] > maxV ? maxV = V[kk] : 0;
  }

  if(maxV < V[mM + mN*M + mP*M*N])
  { return 0; }

  return 1;
}



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

void fArray_flipall(float * T, float * A, int a1, int a2, int a3)
  /* Equivalent to T = flip(flip(flip(A,1),2),3) in matlab */
{
  for(int aa = 0; aa<a1; aa++){
    for(int bb = 0; bb<a2; bb++){
      for(int cc = 0; cc<a3; cc++){
        int idxTo = aa + bb*a1 + cc*a1*a2;
        int idxFrom = (a1-aa) + (a2-bb)*a1 + (a3-cc)*a1*a2;
        T[idxTo] = A[idxFrom];
      }
    }
  }
}

void fArray_insert(float * T, int t1, int t2, int t3, 
    float * F, int f1, int f2, int f3)
  /* Insert F [f1xf2xf3] into T [t1xt2xt3] in the "upper left" corner */
{
  for(int mm = 0; mm<f1; mm++)
  {
    for(int nn = 0; nn<f2; nn++)
    {
      for(int pp = 0; pp<f3; pp++)
      {
        float x = in[mm + nn*f1 + pp*f1*f2];
        out[mm + nn*t1 + pp*t1*t2] = x;
      }
    }
  }
return;
}

void shift_vector(float * V, 
    int S, 
    int N
    int k)
/* Circular shift of a vector of length N with stride S by step k */
{
  assert(abs(k) < N);

  float * buffer = malloc(N*sizeof(float));
    for(size_t pp = 0; pp<N; pp++)
    {
      buffer[pp] = V[pp*S];
    }
    for(size_t pp = 0; pp<N; pp++)
    {
      V[pp*S] = buffer[pp+k % N];
    }
    free(buffer);
    return;
}

void fArray_circshift(float * A, 
    int M, int N, int P, 
    int sm, int sn, int sp)
  /* Shift the image A [MxNxP] by sm, sn, sp in each dimension */
{
// Dimension 1
  for(int bb = 0;bb<N; bb++)
  {
    for(int cc = 0; cc<P; cc++)
    {
    shift_vector(A + bb*M + cc*M*N, 1, M, sm);
  }
}

// Dimension 2
for(aa = 0; aa<M; aa++)
{
  for(int cc = 0;cc<P; cc++)
  {
    shift_vector(A + aa+cc*M*N, M, N, sn);
  }
}

// Dimension 3
for(aa = 0; aa<M; aa++)
{
  for(int bb = 0;bb<N; bb++)
  {
    shift_vector(A + aa+bb*M, M*N, P, sp);
  }
}
return;
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
  fArray_insert(out, M, N, P, in, pM, pN, pP);
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

  if(psfIsCentered(psf, pM, pN, pP) == 0)
  {
    printf("PSF is not centered!\n");
    return NULL;
  }

  // This is the dimensions that we will work with
  // called M1 M2 M3 in the MATLAB code
  int wM = M + pM -1;
  int wN = N + pN -1;
  int wP = P + pP -1;
  size_t wMNP = wM*wN*wP;

  // Prepare the full size cK and cKr
  // (K and Kr in MATLAB)
  float * Z = fftwf_malloc(wMNP*sizeof(float));
  memset(Z, 0, wMNP*sizeof(float));
  fArray_insert(Z, wM, wN, wP, psf, pM, pN, pP);
  float * Zr = fftwf_malloc(wMNP*sizeof(float));
  memset(Zr, 0, wMNP*sizeof(float));
  float * psf_flipped = malloc(pMNP*sizeof(float));
  fArray_flipall(psf_flipped, psf, pM, pN, pP);
  fArray_insert(Zr, wM, wN, wP, psf_flipped, pM, pN, pP);
  free(psf_flipped);
  fArray_circshift(Z, wM, wN, wP, -(pM-1)/2, -(pN-1)/2, -(pP-1)/2);
  fArray_circshift(Zr, wM, wN, wP, -(pM-1)/2, -(pN-1)/2, -(pP-1)/2);

  fftwf_complex * cK = fftwf_malloc(wMNP*sizeof(fftwf_complex));
  fftwf_complex * cKr = fftwf_malloc(wMNP*sizeof(fftwf_complex));


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
      free(im0); 




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



  assert(psf1 != NULL);
  printf("Expanding\n"); fflush(stdout);
  float * psf2 = expandIm_a(psf1, pM, pN, pP, M, N, P);
  free(psf1);


*/

