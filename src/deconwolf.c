#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <string.h>
#include <fftw3.h>
#include <sys/types.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <libgen.h>
#include <time.h>
#include "tiffio.h"
#include "fft.h"
#include "deconwolf.h"
#include "tiling.h"

// use ldd to check where possibly neede
#define INLINED inline __attribute__((always_inline))

opts * opts_new(void)
{
  opts * s = malloc(sizeof(opts));
  s->nThreads = 4;
  s->nIter = 50;
  s->imFile = NULL;
  s->psfFile = NULL;
  s->outFile = NULL;
  s->logFile = NULL;
  s->log = NULL;
  s->verbosity = 1;
  s->tiling_maxSize = -1;
  s->tiling_padding = 20;
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
  free(s);
}

void opts_print(opts * s, FILE *f)
{
  f == NULL ? f = stdout : 0;
  fprintf(f, "Settings:\n");
  fprintf(f, "\t image:  %s\n", s->imFile);
  fprintf(f, "\t psf:    %s\n", s->psfFile);
  fprintf(f, "\t output: %s\n", s->outFile);
  fprintf(f, "\t log file: %s\n", s->logFile);
  fprintf(f, "\t nIter:  %d\n", s->nIter);
  fprintf(f, "\t nThreads: %d\n", s->nThreads);
  fprintf(f, "\t verbosity: %d\n", s->verbosity);
  if(s->tiling_maxSize > 0)
  {
    fprintf(f, "\t tiling, maxSize: %d\n", s->tiling_maxSize);
    fprintf(f, "\t tiling, padding: %d\n", s->tiling_padding);
  } else {
    fprintf(f, "\t tiling: OFF\n");
  }
}

void show_info(FILE * f)
{
  f == NULL ? f = stdout : 0;
  fprintf(f, "Info:\n");
  fprintf(f, "deconwolf: %s\n", deconwolf_version);
  fprintf(f, "FFTW: %s\n", fftwf_version);
  fprintf(f, "TIFF: %s\n", TIFFGetVersion());
  fprintf(f, "running as PID: %d\n", (int) getpid());
  char * host = getenv("HOSTNAME");
  if(host != NULL)
  { fprintf(f, "HOSTNAME: %s\n", host); }
  fprintf(f, "\n");
  fflush(f);
  return;
}

void argparsing(int argc, char ** argv, opts * s)
{

  struct option longopts[] = {
    { "version",     no_argument,       NULL,   'v' }, 
    { "help",         no_argument,       NULL,   'h' },
    // Data
    { "out",        required_argument, NULL,   'o' },
    // Settings
    { "iter",      required_argument, NULL,   'n' },
    { "threads",      required_argument, NULL,   'c' },
    { "verbose",      required_argument, NULL,   'l' },
    { "test",        no_argument,        NULL,   't' },
    { "tilesize",    required_argument,  NULL,   's' }, 
    { "tilepad",     required_argument,  NULL,   'p' },
    { NULL,           0,                 NULL,   0   }
  };

  int ch;
  while((ch = getopt_long(argc, argv, "vho:n:c:p:s:p:", longopts, NULL)) != -1)
  {
    switch(ch) {
      case 'v':
        show_info(NULL);
        exit(0);
        break;
      case 'h':
        usage(argc, argv, s);
        exit(0);
        break;
      case 'o':
        s->outFile = malloc(strlen(optarg)+1);
        strcpy(s->outFile, optarg);
        break;
      case 'n':
        s->nIter = atoi(optarg);
        break;
      case 'c':
        s->nThreads = atoi(optarg);
        break;
      case 'l':
        s->verbosity = atoi(optarg);
        break;
      case 't':
        unittests();
        exit(0);
        break;
      case 's':
        s->tiling_maxSize = atoi(optarg);
        break;
      case 'p':
        s->tiling_padding = atoi(optarg);
        break;
    }
  }

  /* Take care of the positional arguments */
  if(optind + 2 != argc)
  {
    printf("At least image and PSF has to be specified, hint: try '--help'!\n");
    exit(1);
  }

  s->imFile = realpath(argv[optind++], 0);
  s->psfFile = realpath(argv[optind++], 0);

  if(s->psfFile == NULL || s->imFile == NULL)
  {
    printf("realpath can't interpret the file names.\n");
    assert(0);
    exit(1);
  }

  if(s->outFile == NULL)
  {
    char * dirc = strdup(s->imFile);
    char * basec = strdup(s->imFile);
    char * dname = dirname(dirc);
    char * bname = basename(basec);
    s->outFile = malloc(strlen(dname) + strlen(bname) + 10);
    sprintf(s->outFile, "%s/dcw_%s", dname, bname);
    free(dirc);
    free(basec);
  }

  s->logFile = malloc(strlen(s->outFile) + 10);
  sprintf(s->logFile, "%s.log.txt", s->outFile);

  //  printf("Options received\n"); fflush(stdout);
}

void printVmPeak(FILE * fout)
{
  if(fout == NULL) fout = stdout;

  char * statfile = malloc(100*sizeof(char));
  sprintf(statfile, "/proc/%d/status", getpid());
  FILE * sf = fopen(statfile, "r");
  if(sf == NULL)
  {
    fprintf(fout, "Failed to open %s\n", statfile);
    free(statfile);
    return;
  }

  char * line = NULL;
  size_t len = 0;

  while( getline(&line, &len, sf) > 0)
  {
    if(strlen(line) > 6)
    {
      if(strncmp(line, "VmPeak", 6) == 0)
      {
        fprintf(fout, "%s", line);
      }
    }
  }
  free(line);
  fclose(sf);
  free(statfile);
  return;
}

int psfIsCentered(const float * restrict V, const int M, const int N, const int P)
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

float getError(const float * restrict y, const float * restrict g, const int M, const int N, const int P, const int wM, const int wN, const int wP)
{
  double e = 0;
  for(int c = 0; c<P; c++)
  {
    for(int b = 0; b<N; b++)
    {
      for(int a = 0; a<M; a++)
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

float getError_ref(float * y, float * g, int M, int N, int P, int wM, int wN, int wP)
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


float iter(float * xp, const float * restrict g, 
    fftwf_complex * restrict cK,
    fftwf_complex * restrict cKr,
    float * restrict f, 
    float * restrict W,
    float * restrict G, 
    const int wM, const int wN, const int wP,
    const int M, const int N, const int P,
    const opts * s)
{
  const size_t wMNP = wM*wN*wP;

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

fftwf_complex * initial_guess(const int M, const int N, const int P, 
    const int wM, const int wN, const int wP)
{
  /* Create initial guess: the fft of an image that is 1 in MNP and 0 outside
   * M, N, P is the dimension of the microscopic image
   */
  assert(wM > M); assert(wN>N); assert(wP>P);

  float * one = fftwf_malloc(wM*wN*wP*sizeof(float));
  for(size_t kk = 0; kk<wM*wN*wP; kk++)
    one[kk] = 0;

  for(int cc = 0; cc < P; cc++) {
    for(int bb = 0; bb < N; bb++) {
      for(int aa = 0; aa < M; aa++) {
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

void fArray_flipall(float * restrict T, const float * restrict A, const int a1, const int a2, const int a3)
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

void fArray_insert(float * restrict T, const int t1, const int t2, const int t3, 
    const float * restrict F, const int f1, const int f2, const int f3)
  /* Insert F [f1xf2xf3] into T [t1xt2xt3] in the "upper left" corner */
{
  for(int pp = 0; pp<f3; pp++)
  {
    for(int nn = 0; nn<f2; nn++)
    {
      for(int mm = 0; mm<f1; mm++)
      {
        float x = F[mm + nn*f1 + pp*f1*f2];
        T[mm + nn*t1 + pp*t1*t2] = x;
      }
    }
  }
  return;
}

void fArray_insert_ref(float * T, int t1, int t2, int t3, 
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

float * fArray_get_cuboid(float * restrict A, const int M, const int N, const int P,
    const int m0, const int m1, const int n0, const int n1, const int p0, const int p1)
{
  /* Create a new array from V using [m0, m1]x[n0, n1]x[p0, p1] */
  int m = m1-m0+1;
  int n = n1-n0+1;
  int p = p1-p0+1;

  float * C = fftwf_malloc(m*n*p*sizeof(float));
  assert(C != NULL);

  for(int aa = m0; aa <= m1; aa++)
  {
    for(int bb = n0; bb <= n1; bb++)
    {
      for(int cc = p0; cc <= p1; cc++)
      {
        // printf("aa:%d, bb:%d, cc:%d\n", aa, bb, cc);
        size_t Aidx = aa + bb*M + cc*M*N;
        assert(Aidx < M*N*P);
        // New coordinates are offset ...
        size_t Cidx = (aa - m0) + 
          (bb - n0)*m + 
          (cc - p0)*m*n;
        assert(Cidx < m*n*p);
        C[Cidx] = A[Aidx];
      }
    }
  }
  return C;
}

float * fArray_subregion(float * restrict A, const int M, const int N, const int P, const int m, const int n, const int p)
{
  /* Extract sub region starting at (0,0,0) */
  float * S = fftwf_malloc(m*n*p*sizeof(float));
  assert(S != NULL);
  for(int pp = 0; pp<p; pp++)
  {
    for(int nn = 0; nn<n; nn++)
    {
      for(int mm = 0; mm<m; mm++)
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

float * fArray_subregion_ref(float * A, int M, int N, int P, int m, int n, int p)
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


INLINED int mod(const int a, const int b)
{
  int r = a % b;
  return r < 0 ? r + b : r;
}

void shift_vector_buf(float * restrict V, 
    const int S, 
    const int N,
    int k, float * restrict buffer)
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

void shift_vector(float * restrict V, 
    const int S, 
    const int N,
    const int k)
  /* Circular shift of a vector of length N with stride S by step k */
{

  float * buffer = malloc(N*sizeof(float));
  shift_vector_buf(V, S, N, k, buffer);
  free(buffer);
  return;
}

float * fArray_copy(const float * restrict V, const size_t N)
  // Return a newly allocated copy of V
{
  float * C = fftwf_malloc(N*sizeof(float));
  memcpy(C, V, N*sizeof(float));
  return C;
}

float * fArray_constant(const size_t N, const float value)
  // Allocate and return an array of N floats sets to a constant value
{
  float * A = fftwf_malloc(N*sizeof(float));
  for(size_t kk = 0; kk<N; kk++)
  {
    A[kk] = value;
  }
  return A;
}

float * fArray_zeros(const size_t N)
  // Allocate and return an array of N floats
{
  float * A = fftwf_malloc(N*sizeof(float));
  memset(A, 0, N*sizeof(float));
  return A;
}

void fArray_circshift(float * restrict A, 
    const int M, const int N, const int P, 
    const int sm, const int sn, const int sp)
  /* Shift the image A [MxNxP] by sm, sn, sp in each dimension */
{

  const size_t bsize = fmax(fmax(M, N), P);
  float * restrict buf = malloc(bsize*sizeof(float));

  // Dimension 1
    for(int cc = 0; cc<P; cc++)
    {
  for(int bb = 0; bb<N; bb++)
  {    
      //shift_vector(A + bb*M + cc*M*N, 1, M, sm);
      shift_vector_buf(A + bb*M + cc*M*N, 1, M, sm, buf);
    }
  }

  // Dimension 2
    for(int cc = 0; cc<P; cc++)
    {
  for(int aa = 0; aa<M; aa++)
  {    
      //shift_vector(A + aa+cc*M*N, M, N, sn);
      shift_vector_buf(A + aa+cc*M*N, M, N, sn, buf);
    }
  }

  // Dimension 3
    for(int bb = 0; bb<N; bb++)
    {
  for(int aa = 0; aa<M; aa++)
  {  
      //shift_vector(A + aa+bb*M, M*N, P, sp);
      shift_vector_buf(A + aa+bb*M, M*N, P, sp, buf);
    }
  }

  free(buf);
  return;
}

float * expandIm_a(const float * restrict in, 
    const int pM, const int pN, const int pP, 
    const int M, const int N, const int P)
  /* "expand an image" by making it larger 
   * pM, ... current size
   * M, Nm ... new size
   * */
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

void usage(const int argc, char ** argv, const opts * s)
{
  printf(" Usage:\n");
  printf("\t$ %s <options> image.tif psf.tif\n", argv[0]);
  printf("\n");
  printf(" Options:\n");
  printf(" --version\n\t Show version info\n");
  printf(" --help\n\t Show this measage\n");
  printf(" --out <file>\n\t Specify output image name\n");
  printf(" --iter N\n\t Specify the number of iterations to use (default: %d)\n", s->nIter);
  printf(" --threads N\n\t Specify the number of threads to use\n");
  printf(" --verbose N\n\t Set verbosity level (default: %d)\n", s->verbosity);
  printf(" --test\n\t Run unit tests\n");
  printf(" --tilesize N\n\t Enables tiling mode and sets the largest tile size to N voxels in x and y.\n");
  printf(" --tilepad N\n\t Sets the tiles to overlap by N voxels in tile mode (default: %d)\n", s->tiling_padding);
  printf("\n");
}


float update_alpha(const float * restrict g, const float * restrict gm, const size_t wMNP)
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

float * deconvolve(const float * restrict im, const int M, const int N, const int P,
    const float * restrict psf, const int pM, const int pN, const int pP,
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
  }
  if(s->verbosity > 1)
  {
    printf("Estimated peak memory usage: %.1f GB\n", wMNP*65.0/1e9);
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
  float * y = fArray_copy(f, wMNP);

  float * x1 = fArray_copy(f, wMNP);
  float * x2 = fArray_copy(f, wMNP);
  float * x = x1;
  float * xp = x2;
  float * xm = x2; // fArray_copy(f, wMNP);

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

    xp = xm;
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

    float * swap = g;
    g = gm; gm = swap;

    for(size_t kk = 0; kk<wMNP; kk++)
    { 
      g[kk] = xp[kk]-y[kk]; 
    }

    if(it > 0) {
      alpha = update_alpha(g, gm, wMNP);
    }

    xm = x;
    x = xp;
    xp = NULL;

    it++;
  }

  fftwf_free(W); // is P1
  fftwf_free(G);
  float * out = fArray_subregion(x, wM, wN, wP, M, N, P);
  fftwf_free(f);
  fftwf_free(x1);
  fftwf_free(x2);
  fftwf_free(g);
  fftwf_free(gm);
  fftwf_free(cK);
  fftwf_free(cKr);
  fftwf_free(y);
  return out;
}

float * deconvolve_tiles(const float * restrict im, int M, int N, int P,
    const float * restrict psf, const int pM, const int pN, const int pP,
    opts * s)
{

  tiling * T = tiling_create(M,N,P, s->tiling_maxSize, s->tiling_padding);
  if( T == NULL)
  {
    printf("Tiling failed, please check your settings\n");
    exit(1);
  }

  if( T->nTiles == 1)
  {
    printf("\n"
        "ERROR: Only one tile! Please omit the `--tilesize` parameter if "
        "that is what you intended to to or decrease the value if you "
        "want to process the image in tiles."
        "\n\n");
    exit(1);
  }

  if(s->verbosity > 0)
  {
    printf("-> Divided the [%d x %d x %d] image into %d tiles\n", M, N, P, T->nTiles);
  }

  // Output image
  float * V = malloc(M*N*P*sizeof(float));
  memset(V, 0, M*N*P*sizeof(float));


  for(int tt = 0; tt<T->nTiles; tt++)
  {
    // Temporal copy of the PSF that might be cropped to fit the tile
    float * tpsf = fArray_copy(psf, pM*pN*pP);
    int tpM = pM, tpN = pN, tpP = pP;

    if(s->verbosity > 0)
    {
      printf("-> Processing tile %d / %d\n", tt+1, T->nTiles);
      fprintf(s->log, "-> Processing tile %d / %d\n", tt+1, T->nTiles);
    }

    float * im_tile = tiling_get_tile(T, tt, im);

    int tileM = T->tiles[tt]->xsize[0];
    int tileN = T->tiles[tt]->xsize[1];
    int tileP = T->tiles[tt]->xsize[2];

   tpsf = autocrop_psf(tpsf, &tpM, &tpN, &tpP, 
      tileM, tileN, tileP, s);

    float * dw_im_tile = deconvolve(im_tile, tileM, tileN, tileP, // input image and size
        tpsf, tpM, tpN, tpP, // psf and size
        s);

    tiling_put_tile(T, tt, V, dw_im_tile);
    free(im_tile);
    free(dw_im_tile);
    free(tpsf);
  }
  tiling_free(T);
  free(T);
  return V;
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

static double timespec_diff(struct timespec* end, struct timespec * start)
{
  double elapsed = (end->tv_sec - start->tv_sec);
  elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
  return elapsed;
}

#define tictoc struct timespec tictoc_start, tictoc_end;
#define tic clock_gettime(CLOCK_REALTIME, &tictoc_start);
#define toc(X) clock_gettime(CLOCK_REALTIME, &tictoc_end); printf(#X); printf(" %f s\n", timespec_diff(&tictoc_end, &tictoc_start));

void timings()
{
  printf("-> Timings\n");
  tictoc
    int M = 1024, N = 1024, P = 50;
  float temp = 0;

  tic
    usleep(1000);
  toc(usleep_1000)

    tic
    float * V = malloc(M*N*P*sizeof(float));
  toc(malloc)

    float * A = malloc(M*N*P*sizeof(float));

  tic
    memset(V, 0, M*N*P*sizeof(float));
  toc(memset)
    memset(A, 0, M*N*P*sizeof(float));
  for(int kk = 0; kk<M*N*P; kk++)
  {
    A[kk] = (float) rand()/(float) RAND_MAX;
  }

  // ---
  tic
    fftwf_plan p = fftwf_plan_dft_r2c_3d(P, N, M, 
        V, NULL, 
        FFTW_WISDOM_ONLY | FFTW_MEASURE);
  fftwf_destroy_plan(p);
  toc(fftwf_plan_create_and_destroy)


    // ---
    tic
    fArray_flipall(V, A, M, N, P);
  toc(fArray_flipall)

    // ---
    tic 
    temp = update_alpha(V, A, M*N*P);
  toc(update_alpha)
    V[0]+= temp;

  // ---
  tic
    float e1 = getError_ref(V, A, M, N, P, M, N, P);
  toc(getError_ref)
    V[0]+= e1;

  tic
    float e2 = getError(V, A, M, N, P, M, N, P);
  toc(getError)
    V[0]+=e2;

  printf("e1: %f, e2: %f, fabs(e1-e2): %f\n", e1, e1, fabs(e1-e2));

  // ---
  tic
    float * S1 = fArray_subregion(V, M, N, P, M-1, N-1, P-1);
  toc(fArray_subregion)

    tic
    float * S2 = fArray_subregion_ref(V, M, N, P, M-1, N-1, P-1);
  toc(fArray_subregion_ref)
    printf("S1 - S2 = %f\n", getError(S1, S1, M-1, N-1, P-1, M-1, N-1, P-1));
  free(S1);
  free(S2);

  // ---
  tic
    fArray_insert(V, M, N, P, A, M-1, N-1, P-1);
  toc(fArray_subregion)      

    tic
    fArray_insert_ref(V, M, N, P, A, M-1, N-1, P-1);
  toc(fArray_subregion_ref)

    // ---


    ((float volatile *)V)[0] = V[0];
  printf("V[0] = %f\n", V[0]);
}


void unittests()
{
  timings();
  shift_vector_ut();
  fArray_flipall_ut();
  //
}

void show_time(FILE * f)
{
  f == NULL ? f = stdout : 0;
  time_t now = time(NULL);
  char * tstring = ctime(&now);
  fprintf(f, "%s\n", tstring);
}

float * autocrop_psf(float * psf, int * pM, int * pN, int * pP,  // psf and size
    int M, int N, int P, // image size
    opts * s)
{
  int m = pM[0];
  int n = pN[0];
  int p = pP[0];

  if((p % 2) == 0)
  {
    printf("Error: The PSF should have odd number of slices\n");
    exit(1);
  }

  // Optimal size
  int mopt = (M-1)*2 + 1;
  int nopt = (N-1)*2 + 1;
  int popt = (P-1)*2 + 1;

  if(p < popt)
  {
    fprintf(s->log, "The PSF does seem to have too few slices\n");
    return psf;
  }

  if(m > mopt || n > nopt || p > popt)
  { 
    int m0 = 0, m1 = m-1;
    int n0 = 0, n1 = n-1;
    int p0 = 0, p1 = p-1;
    if(m > mopt)
    {
      m0 = (m-mopt)/2-1;
      m1 = m1-(m-mopt)/2-1;
    }
    if(n > nopt)
    {
      n0 = (n-nopt)/2-1;
      n1 = n1-(n-nopt)/2-1;
    }
    if(p > popt)
    {
      p0 = (p-popt)/2-1;
      p1 = p1-(p-popt)/2-1;
    }

//    printf("! %d %d : %d %d : %d %d\n", m0, m1, n0, n1, p0, p1);
    float * psf_cropped = fArray_get_cuboid(psf, m, n, p,
        m0, m1, n0, n1, p0, p1);
    free(psf);

    pM[0] = m1-m0+1;
    pN[0] = n1-n0+1;
    pP[0] = p1-p0+1;

    if(s->verbosity > 0)
    {
      fprintf(stdout, "The PSF was cropped to [%d x %d x %d]\n", pM[0], pN[0], pP[0]);
    }
    fprintf(s->log, "The PSF was cropped to [%d x %d x %d]\n", pM[0], pN[0], pP[0]);
    return psf_cropped;

  } else {
    return psf;
  }

}

void dcw_init_log(opts * s)
{
  s->log = fopen(s->logFile, "w");
  assert(s->log != NULL);
  show_time(s->log);
  opts_print(s, s->log); 
  show_info(s->log);
}

void dcw_close_log(opts * s)
{
  printVmPeak(s->log);
  show_time(s->log);
  fclose(s->log);
}


int deconwolf(opts * s)
{
  dcw_init_log(s);

  if(s->verbosity > 0) 
  {
    opts_print(s, NULL); 
    printf("\n");
  }

  s->verbosity > 1 ? show_info(NULL) : 0;

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

  // Possibly the PSF will be cropped even more per tile later on
   psf = autocrop_psf(psf, &pM, &pN, &pP, 
      M, N, P, s);

  myfftw_start(s->nThreads);
  float * out = NULL;
  if(s->tiling_maxSize < 0)
  {
   out = deconvolve(im, M, N, P, // input image and size
        psf, pM, pN, pP, // psf and size
        s);// settings
  } else {
    out = deconvolve_tiles(im, M, N, P, // input image and size
        psf, pM, pN, pP, // psf and size
        s);// settings
  }
  // floatimage_show_stats(out, M, N, P);
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
  if(s->verbosity > 1) printVmPeak(NULL);
  dcw_close_log(s);
  opts_free(&s);

  return(exitstatus);
}

int main(int argc, char ** argv)
{
  opts * s = opts_new(); // Load default settings and initialize
  argparsing(argc, argv, s); // Parse command line
  return deconwolf(s); // And go!
}

