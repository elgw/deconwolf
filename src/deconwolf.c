#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <fftw3.h>
#include <sys/types.h>
#include <unistd.h>
#include "tiffio.h"
#include "fft.h"
#include "deconwolf.h"

typedef struct{
  int nThreads;
  int nIter;
  char * imFile;
  char * psfFile;
  char * outFile;
  char * logFile;
  FILE * log;

  int verbosity;
  fftwf_plan fft_plan;
  fftwf_plan ifft_plan;
} opts;

opts * opts_new(void)
{
  opts * s = malloc(sizeof(opts));
  s->nThreads = 8;
  s->nIter = 50;
  s->imFile = NULL;
  s->psfFile = NULL;
  s->outFile = NULL;
  s->logFile = NULL;
  s->log = NULL;
  s->verbosity = 1;
  return s;
}

void nullfree(void * p)
{
  if(p != NULL)
  {
    free(p);
  }
}

void opts_free(opts ** sp)
{
  opts * s = sp[0];
  nullfree(s->imFile);
  nullfree(s->psfFile);
  nullfree(s->outFile);
  nullfree(s->logFile);
  fclose(s->log);
  free(s);
}

void opts_print(opts * s, FILE *f)
{
  f == NULL ? f = stdout : 0;
  fprintf(f, "Settings:\n");
  fprintf(f, "\t image:  %s\n", s->imFile);
  fprintf(f, "\t psf:    %s\n", s->psfFile);
  fprintf(f, "\t output: %s\n", s->outFile);
  fprintf(f, "\t nIter:  %d\n", s->nIter);
  fprintf(f, "\t nThreads: %d\n", s->nThreads);
  fprintf(f, "\t verbosity: %d\n", s->verbosity);
  fprintf(f, "\t log file: %s\n", s->logFile);
}

void show_info(FILE * f)
{
  f == NULL ? f = stdout : 0;
  fprintf(f, "Info:\n");
  fprintf(f, "\t PID: %d\n", (int) getpid());
  fprintf(f, "\t FFTW: %s\n", fftwf_version);
  fprintf(f, "\t TIFF: %s\n", TIFFGetVersion());
  fprintf(f, "\n");
  return;
}

void argparsing(int argc, char ** argv, opts * s)
{
  if(argc < 3)
  {
    usage(argc, argv);
    unittests();
    exit(1);
  }

  s->imFile = malloc(strlen(argv[1])+1);
  strcpy(s->imFile, argv[1]);

  s->psfFile = malloc(strlen(argv[2])+1);
  strcpy(s->psfFile, argv[2]);
  s->outFile = malloc(100*sizeof(char));
  sprintf(s->outFile, "dummy.tif");

  s->logFile = malloc(100*sizeof(char));
  sprintf(s->logFile, "%s.log.txt", s->outFile);
  s->log = fopen(s->logFile, "w");
  assert(s->log != NULL);
}

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

float getError(float * y, float * g, int M, int N, int P, int wM, int wN, int wP)
{
  double e = 0;
  for(int a = 0; a<M; a++)
  {
    for(int b = 0; b<N; b++)
    {
      for(int c = 0; c<P; c++)
      {
        double yval = y[a + b*wM + c*wM*wN];
        double gval = g[a + b*M + c * M*N];
        e+=pow(yval-gval, 2);
      }
    }
  }
  e/=(M*N*P);
  return (float) e;
}


float iter(float * xp, float * g, 
    fftwf_complex * cK,
    fftwf_complex * cKr,
    float * f, 
    float * W,
    float * G, 
    int wM, int wN, int wP,
    int M, int N, int P,
    opts * s)
{
  size_t wMNP = wM*wN*wP;

  fftwf_complex * F = fft(f, wM, wN, wP);
  float * y = fft_convolve_cc(cK, F, wM, wN, wP);
  fftwf_free(F);
  float error = getError(y, g, M, N, P, wM, wN, wP);
  float * sn  = y; // alias
  for(size_t kk = 0; kk<wMNP; kk++)
  {
    y[kk] < 1e-6 ? y[kk] = 1e-6 : 0;
    sn[kk] = G[kk]/y[kk];
  }
  fftwf_complex * F_sn = fft(sn, wM, wN, wP);
  fftwf_free(y); // = sn as well
  float * x = fft_convolve_cc(cKr, F_sn, wM, wN, wP);
  fftwf_free(F_sn);
  for(size_t kk = 0; kk<wMNP; kk++)
  {
    xp[kk] = f[kk]*x[kk]*W[kk];
  }
  fftwf_free(x);
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
  float sum = 0;
  float max = 0;
  for(size_t kk = 0; kk<wM*wN*wP; kk++)
  {
    sum+=one[kk];
    one[kk] > max ? max = one[kk] : 0;
  }

  //  printf("sum(one) = %f max(one) = %f\n", sum, max);

  fftwf_complex * Fone = fft(one, wM, wN, wP);

  fftwf_free(one);
  return Fone;
}

void fArray_stats(float * A, size_t N)
{
  float amax = 0;
  for(size_t kk = 0; kk<N; kk++)
  {
    if(A[kk] > amax)
      amax = A[kk];
  }
  printf("max: %f\n", amax);
}

void fArray_flipall(float * T, float * A, int a1, int a2, int a3)
  /* Equivalent to T = flip(flip(flip(A,1),2),3) in matlab */
{
  for(int aa = 0; aa<a1; aa++){
    for(int bb = 0; bb<a2; bb++){
      for(int cc = 0; cc<a3; cc++){
        int idxTo = aa + bb*a1 + cc*a1*a2;
        int idxFrom = (a1-aa-1) + (a2-bb-1)*a1 + (a3-cc-1)*a1*a2;

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
        float x = F[mm + nn*f1 + pp*f1*f2];
        T[mm + nn*t1 + pp*t1*t2] = x;
      }
    }
  }
  return;
}

float * fArray_subregion(float * A, int M, int N, int P, int m, int n, int p)
{
  float * S = fftwf_malloc(m*n*p*sizeof(float));
  assert(S != NULL);
  for(int mm = 0; mm<m; mm++)
  {
    for(int nn = 0; nn<n; nn++)
    {
      for(int pp = 0; pp<p; pp++)
      {
        size_t Aidx = mm + nn*M + pp*M*N;
        size_t Sidx = mm + nn*m + pp*m*n;
        assert(Aidx < M*N*P);
        assert(Sidx < m*n*p);
        S[Sidx] = A[Aidx];
      }
    }
  }
  return S;
}

int mod(int a, int b)
{
  int r = a % b;
  return r < 0 ? r + b : r;
}



void shift_vector_buf(float * V, 
    int S, 
    int N,
    int k, float * buffer)
  /* Circular shift of a vector of length N with stride S by step k */
{

  k = -k;
  for(size_t pp = 0; pp<N; pp++)
  {
    buffer[pp] = V[pp*S];
  }
  for(size_t pp = 0; pp<N; pp++)
  {
    V[pp*S] = buffer[mod(pp+k, N)];
  }
  return;
}

void shift_vector(float * V, 
    int S, 
    int N,
    int k)
  /* Circular shift of a vector of length N with stride S by step k */
{

  float * buffer = malloc(N*sizeof(float));
  shift_vector_buf(V, S, N, k, buffer);
  free(buffer);
  return;
}

float * fArray_copy(float * V, size_t N)
  // Return a newly allocated copy of V
{
  float * C = fftwf_malloc(N*sizeof(float));
  memcpy(C, V, N*sizeof(float));
  return C;
}

float * fArray_constant(size_t N, float value)
  // Allocate and return an array of N floats
{
  float * A = fftwf_malloc(N*sizeof(float));
  for(size_t kk = 0; kk<N; kk++)
  {
    A[kk] = value;
  }
  return A;
}
float * fArray_zeros(size_t N)
  // Allocate and return an array of N floats
{
  float * A = fftwf_malloc(N*sizeof(float));
  memset(A, 0, N*sizeof(float));
  return A;
}

void fArray_circshift(float * A, 
    int M, int N, int P, 
    int sm, int sn, int sp)
  /* Shift the image A [MxNxP] by sm, sn, sp in each dimension */
{

  size_t bsize = fmax(fmax(M, N), P);
  float * buf = malloc(bsize*sizeof(float));

  // Dimension 1
  for(int bb = 0; bb<N; bb++)
  {
    for(int cc = 0; cc<P; cc++)
    {
      //shift_vector(A + bb*M + cc*M*N, 1, M, sm);
      shift_vector_buf(A + bb*M + cc*M*N, 1, M, sm, buf);
    }
  }

  // Dimension 2
  for(int aa = 0; aa<M; aa++)
  {
    for(int cc = 0; cc<P; cc++)
    {
      //shift_vector(A + aa+cc*M*N, M, N, sn);
      shift_vector_buf(A + aa+cc*M*N, M, N, sn, buf);
    }
  }

  // Dimension 3
  for(int aa = 0; aa<M; aa++)
  {
    for(int bb = 0; bb<N; bb++)
    {
      //shift_vector(A + aa+bb*M, M*N, P, sp);
      shift_vector_buf(A + aa+bb*M, M*N, P, sp, buf);
    }
  }

  free(buf);
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


float update_alpha(float * g, float * gm, size_t wMNP)
{
  double ggm = 0;
  double gmgm = 0;
  for(size_t kk = 0; kk<wMNP; kk++)
  {
    ggm += g[kk]*gm[kk];
    gmgm += pow(gm[kk], 2);
  }
  float alpha = ggm/(gmgm+1e-7);
  alpha < 0 ? alpha = 0 : 0;
  alpha > 1 ? alpha = 1 : 0;
  return alpha;
}

float * deconvolve(float * im, int M, int N, int P,
    float * psf, int pM, int pN, int pP,
    opts * s)
{
  if(s->verbosity > 0)
  {
    printf("Deconvolving\n");
  }

  /*Deconvolve im using psf */
  const int nIter = s->nIter;

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

  if(s->verbosity > 0)
  { printf("image: [%dx%dx%d], psf: [%dx%dx%d], job: [%dx%dx%d] (%zu voxels)\n",
      M, N, P, pM, pN, pP, wM, wN, wP, wMNP);
  printf("Estimated peak memory use: %.1f GB\n", wMNP*65.0/1e9);
  }
  fprintf(s->log, "image: [%dx%dx%d], psf: [%dx%dx%d], job: [%dx%dx%d] (%zu voxels)\n",
      M, N, P, pM, pN, pP, wM, wN, wP, wMNP);
  fflush(s->log);

  fft_train(wM, wN, wP, s->verbosity);


  // Prepare the full size cK and cKr
  // (K and Kr in MATLAB)
  float * Z = fftwf_malloc(wMNP*sizeof(float));
  memset(Z, 0, wMNP*sizeof(float));

  //printf("insert\n"); fflush(stdout);
  fArray_insert(Z, wM, wN, wP, psf, pM, pN, pP);
  float * Zr = fftwf_malloc(wMNP*sizeof(float));
  memset(Zr, 0, wMNP*sizeof(float));
  float * psf_flipped = malloc(wMNP*sizeof(float));
  memset(psf_flipped, 0, sizeof(float)*wMNP);
  //fArray_stats(psf, pM*pN*pP);
  fArray_flipall(psf_flipped, psf, pM, pN, pP);
  //fArray_stats(psf_flipped, pM*pN*pP);
  //writetif("psf_flipped.tif", psf_flipped, pM, pN, pP);
  fArray_insert(Zr, wM, wN, wP, psf_flipped, pM, pN, pP);
  //writetif("Zr0.tif", Zr, wM, wN, wP);
  free(psf_flipped);
  fArray_circshift(Z, wM, wN, wP, -(pM-1)/2, -(pN-1)/2, -(pP-1)/2);
  fArray_circshift(Zr, wM, wN, wP, -(pM-1)/2, -(pN-1)/2, -(pP-1)/2);
  //printf("...\n"); fflush(stdout);
  //fArray_stats(Z, pM*pN*pP);
  //fArray_stats(Zr, pM*pN*pP);
  //writetif("Z.tif", Z, wM, wN, wP);
  //writetif("Zr.tif", Zr, wM, wN, wP);

  fftwf_complex * cK = fft(Z, wM, wN, wP);
  fftwf_free(Z);
  fftwf_complex * cKr = fft(Zr, wM, wN, wP);
  fftwf_free(Zr);
  //printf("initial guess\n"); fflush(stdout);
  fftwf_complex * F_one = initial_guess(M, N, P, wM, wN, wP);
  float * P1 = fft_convolve_cc(cKr, F_one, wM, wN, wP);
  fftwf_free(F_one);
  //printf("P1\n");
  //fArray_stats(P1, pM*pN*pP);

  //writetif("P1.tif", P1, wM, wN, wP);
  float sigma = 0.001;
  for(size_t kk = 0; kk<wMNP; kk++)
  {
    if(P1[kk] < sigma)
    { P1[kk] = 0; } else { P1[kk] = 1/P1[kk]; }
  }
  float * W = P1;
  //writetif("W.tif", W, wM, wN, wP);

  // Original image -- expanded
  float * G = fftwf_malloc(wMNP*sizeof(float));
  memset(G, 0, wMNP*sizeof(float));
  fArray_insert(G, wM, wN, wP, im, M, N, P);
  //writetif("G.tif", G, wM, wN, wP);

  float sumg = 0;
  for(size_t kk = 0; kk<M*N*P; kk++)
  {
    sumg+=im[kk];
  }

  float alpha = 0;
  float * f = fArray_constant(wMNP, sumg/wMNP);
  float * x = fArray_copy(f, wMNP);
  float * y = fArray_copy(f, wMNP);
  float * xp = fArray_copy(f, wMNP);
  float * xm = fArray_copy(f, wMNP);
  float * gm = fArray_zeros(wMNP);
  float * g = fArray_zeros(wMNP);

  int it = 0;
  while(it<nIter)
  {

    for(size_t kk = 0; kk<wMNP; kk++)
    { 
      y[kk] = x[kk] + alpha*(x[kk]-xm[kk]);
      y[kk] < 0 ? y[kk] = 0 : 0;
    }

    // printf("y[0] = %f\n", y[0]);

    double err = iter(xp, // xp is updated
        im, 
        cK, cKr, 
        y, W, G, 
        wM, wN, wP, 
        M, N, P, s);

    if(s->verbosity >0){
      printf("Iteration %d/%d, error=%e\n", it+1, nIter, err);}
    if(s->log != NULL)
    { fprintf(s->log, "Iteration %d/%d, error=%e\n", it+1, nIter, err);}

    memcpy(gm, g, wMNP*sizeof(float)); // swap points instead
    //    printf("xp[0] = %f, y[0] = %f\n", xp[0], y[0]);
    for(size_t kk = 0; kk<wMNP; kk++)
    { 
      g[kk] = xp[kk]-y[kk]; 
    }
    //  if(it == 0)
    //   writetif("g_iter0.tif",g, wM, wN, wP);
    if(it > 0) {
      alpha = update_alpha(g, gm, wMNP);
      // printf("alpha=%f\n", alpha);
    }
    memcpy(xm, x, wMNP*sizeof(float)); // swap pointers instead
    memcpy(x, xp, wMNP*sizeof(float)); // really need tree of them?
    it++;
    //writetif("xp.tif", xp, wM, wN, wP);
  }
  fftwf_free(W); // is P1
  fftwf_free(G);
  float * out = fArray_subregion(xp, wM, wN, wP, M, N, P);
  fftwf_free(x);
  fftwf_free(f);
  fftwf_free(xp);
  fftwf_free(xm);
  fftwf_free(g);
  fftwf_free(gm);
  fftwf_free(cK);
  fftwf_free(cKr);
  fftwf_free(y);
  return out;
}

void fArray_flipall_ut()
{

  float * a = malloc(3*3*3*sizeof(float));
  float * b = malloc(3*3*3*sizeof(float));
  float * c = malloc(3*3*3*sizeof(float));

  for(int kk = 0; kk<27; kk++)
  {
    a[kk] = 0;
  }

  a[13] = 1;
  fArray_flipall(b, a, 3, 3, 3);
  assert(b[13] == 1);

  for(int kk = 0; kk<27; kk++)
  {
    a[kk] = rand();
  }

  fArray_flipall(b, a, 3, 3, 3);
  fArray_flipall(c, b, 3, 3, 3);
  for(int kk = 0; kk<27; kk++)
  {
    assert(a[kk] == c[kk]);
  }

  fArray_flipall(b, a, 4, 3, 2);
  fArray_flipall(c, b, 4, 3, 2);
  for(int kk = 0; kk<24; kk++)
    assert(a[kk] == c[kk]);

  fArray_flipall(b, a, 2, 3, 4);
  fArray_flipall(c, b, 2, 3, 4);
  for(int kk = 0; kk<24; kk++)
    assert(a[kk] == c[kk]);

  free(a); free(b); free(c);
  return;
}

void shift_vector_ut()
{
  int N = 5;
  int S = 1; // stride
  float * V = malloc(N*sizeof(float));

  for(int k = -7; k<7; k++)
  {
    for(int kk = 0; kk<N; kk++)
    {V[kk] = kk;}
    printf("shift: % d -> ", k);
    shift_vector(V,S,N,k);
    for(int kk =0; kk<N; kk++)
    { printf("%.0f ", V[kk]);}
    printf("\n");
  }
}

void unittests()
{
  shift_vector_ut();
  fArray_flipall_ut();
  //
}


int main(int argc, char ** argv)
{

  opts * s = opts_new(); 
  argparsing(argc, argv, s);

  if(s->verbosity > 0) 
  {
    opts_print(s, NULL); 
  }
  if(s->verbosity > 1)
  {
    show_info(NULL);
  }

  if(s->verbosity > 0)
  {
    printf("Reading %s\n", s->imFile);
  }
  int M = 0, N = 0, P = 0;
  float * im = readtif_asFloat(s->imFile, &M, &N, &P, s->verbosity);
  // writetif("im.tif", im, M, N, P);

  int pM = 0, pN = 0, pP = 0;
  float * psf = NULL;
  if(1){
    if(s->verbosity > 0)
    {
      printf("Reading %s\n", s->psfFile);
    }
    psf = readtif_asFloat(s->psfFile, &pM, &pN, &pP, s->verbosity);
    if(psf == NULL)
    {
      printf("Failed to open %s\n", s->psfFile);
      exit(1);
    }
  } else {
    pM = 3; pN = 3; pP = 3;
    psf = malloc(27*sizeof(float));
    memset(psf, 0, 27*sizeof(float));
    psf[13] = 1;
  }

  myfftw_start(s->nThreads);
  float * out = NULL;
  out = deconvolve(im, M, N, P, // input image and size
      psf, pM, pN, pP, // psf and size
      s);// settings

  int exitstatus = 1;
  if(out == NULL)
  {
    if(s->verbosity > 0)
    {
      printf("Nothing to write to disk :(\n");
    }
  } else 
  {
    if(s->verbosity > 0)
    {
      printf("Writing to %s\n", s->outFile); fflush(stdout);
    }
    //    floatimage_normalize(out, M*N*P);
    writetif(s->outFile, out, M, N, P);
    exitstatus = 0;
  }

  if(s->verbosity > 0)
  {
    printf("Finalizing\n"); fflush(stdout);
  }
  free(im);
  free(psf);
  free(out);
  myfftw_stop();
  opts_free(&s);
  if(0)
  {
  printf("Do a grep VmPeak /proc/%d/status now ...\n", getpid());
  getchar();
  }
  return(exitstatus);
}


// Graveyard

/* 
 *
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

