/*    Copyright (C) 2020 Erik L. G. Wernersson
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

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
#include "fft.h"
#include "dw.h"
#include "tiling.h"
#include "fim.h"
#include "fim_tiff.h"

dw_opts * dw_opts_new(void)
{
  dw_opts * s = malloc(sizeof(dw_opts));
  s->nThreads = 4;
  s->nIter = 25;
  s->imFile = NULL;
  s->psfFile = NULL;
  s->outFile = NULL;
  s->logFile = NULL;
  s->prefix = malloc(10*sizeof(char));
  sprintf(s->prefix, "dw");
  s->log = NULL;
  s->verbosity = 1;
  s->overwrite = 0;
  s->tiling_maxSize = -1;
  s->tiling_padding = 20;
  s->method = DW_METHOD_W;
  s->iterdump = 0;
  s->relax = 0;
  return s;
}

void nullfree(void * p)
{
  if(p != NULL)
  {
    free(p);
  }
}

void dw_opts_free(dw_opts ** sp)
{
  dw_opts * s = sp[0];
  nullfree(s->imFile);
  nullfree(s->psfFile);
  nullfree(s->outFile);
  nullfree(s->logFile);
  nullfree(s->prefix);
  free(s);
}

void dw_opts_fprint(FILE *f, dw_opts * s)
{
  f == NULL ? f = stdout : 0;
  fprintf(f, "> Settings:\n");
  fprintf(f, "image:  %s\n", s->imFile);
  fprintf(f, "psf:    %s\n", s->psfFile);
  fprintf(f, "output: %s\n", s->outFile);
  fprintf(f, "log file: %s\n", s->logFile);
  fprintf(f, "nIter:  %d\n", s->nIter);
  fprintf(f, "nThreads: %d\n", s->nThreads);
  fprintf(f, "verbosity: %d\n", s->verbosity);
  switch(s->method)
  { 
    case DW_METHOD_RL:
      fprintf(f, "method: Richardson-Lucy\n");
      break;
    case DW_METHOD_ID:
      fprintf(f, "method: Identity\n");
  }
  if(s->overwrite == 1)
  { fprintf(f, "overwrite: YES\n"); } else
  { fprintf(f, "overwrite: NO\n"); }

  if(s->tiling_maxSize > 0)
  {
    fprintf(f, "tiling, maxSize: %d\n", s->tiling_maxSize);
    fprintf(f, "tiling, padding: %d\n", s->tiling_padding);
  } else {
    fprintf(f, "tiling: OFF\n");
  }
  if(s->relax > 0)
  {
    fprintf(f, "PSF relaxation: %f\n", s->relax);
  }
  fprintf(f, "\n");
}

int file_exist(char * fname)
{
  if( access( fname, F_OK ) != -1 ) {
    return 1; // File exist
  } else {
    return 0;
  }
}

void dw_fprint_info(FILE * f)
{
  f == NULL ? f = stdout : 0;
  fprintf(f, "deconwolf: '%s' PID: %d\n", deconwolf_version, (int) getpid());
#ifdef GIT_VERSION
  fprintf(f, "GIT_VERSION: '%s'\n", GIT_VERSION);
#endif
#ifdef CC_VERSION
  fprintf(f, "COMPILER: '%s'\n", CC_VERSION);
#endif
  fprintf(f, "BUILD_DATE: '%s'\n'", __DATE__);
  fprintf(f, "FFTW: '%s'\n", fftwf_version);
  fprintf(f, "TIFF: '%s'\n", TIFFGetVersion());
  char * user = getenv("USER");
  if(user != NULL)
  { fprintf(f, "USER: '%s'\n", user); }

  char * hname = malloc(1024*sizeof(char));
  if(gethostname(hname, 1023) == 0)
  { 
    fprintf(f, "HOSTNAME: '%s'\n", hname); 
    free(hname);
  }

  fprintf(f, "\n");
  fflush(f);
  return;
}

void deconwolf_batch(dw_opts * s)
{
  // Better to do in Lua/Python/bash/...
  printf("Not implemented !\n");
}

void dw_argparsing(int argc, char ** argv, dw_opts * s)
{

  int generate_batch = 0;

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
    { "overwrite",   no_argument,        NULL,   'w' },
    { "prefix",       required_argument, NULL,   'f' },
    { "batch",        no_argument, NULL, 'b' },
    { "method",       required_argument, NULL,   'm' },
    { "relax",        required_argument, NULL,   'r' },
    { NULL,           0,                 NULL,   0   }
  };

  int ch;
  while((ch = getopt_long(argc, argv, "vho:n:c:p:s:p:", longopts, NULL)) != -1)
  {
    switch(ch) {
      case 'v':
        dw_fprint_info(NULL);
        exit(0);
        break;
      case 'h':
        dw_usage(argc, argv, s);
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
        dw_unittests();
        exit(0);
        break;
      case 's':
        s->tiling_maxSize = atoi(optarg);
        break;
      case 'p':
        s->tiling_padding = atoi(optarg);
        break;
      case 'w':
        s->overwrite = 1;
        break;
      case 'f':
        free(s->prefix);
        s->prefix = malloc(strlen(optarg) + 1);
        strcpy(s->prefix, optarg);
        break;
      case 'b':
        generate_batch = 1;
        break;
      case 'm':
        if(strcmp(optarg, "rl") == 0)
        { 
          s->method = DW_METHOD_RL; 
          sprintf(s->prefix, "drl");
        }
        if(strcmp(optarg, "id") == 0)
        { 
          s->method = DW_METHOD_ID; 
          sprintf(s->prefix, "id");
        }
        break;
      case 'r':
        s->relax = atof(optarg);
        break;
    }
  }

  /* Take care of the positional arguments */
  if(optind + 2 != argc)
  {
    printf("At least image and PSF has to be specified, hint: try '--help'!\n");
    exit(1);
  }

  if(generate_batch)
  {
    deconwolf_batch(s);
    exit(0);
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
    sprintf(s->outFile, "%s/%s_%s", dname, s->prefix, bname);
    free(dirc);
    free(basec);
  }

  if( s->overwrite == 0 && file_exist(s->outFile))
  {
    printf("%s already exist. Doing nothing\n", s->outFile);
    exit(0);
  }

  s->logFile = malloc(strlen(s->outFile) + 10);
  sprintf(s->logFile, "%s.log.txt", s->outFile);

  //  printf("Options received\n"); fflush(stdout);
}

#ifdef __APPLE__
size_t get_peakMemoryKB(void)
{
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);
  return (size_t) round((double) r_usage.ru_maxrss/1024.0);
}
#endif

#ifndef __APPLE__
size_t get_peakMemoryKB(void)
{
  char * statfile = malloc(100*sizeof(char));
  sprintf(statfile, "/proc/%d/status", getpid());
  FILE * sf = fopen(statfile, "r");
  if(sf == NULL)
  {
    fprintf(stderr, "Failed to open %s\n", statfile);
    free(statfile);
    return 0;
  }

  char * peakline = NULL;

  char * line = NULL;
  size_t len = 0;

  while( getline(&line, &len, sf) > 0)
  {
    if(strlen(line) > 6)
    {
      if(strncmp(line, "VmPeak", 6) == 0)
      {
        peakline = strdup(line);
      }
    }
  }
  free(line);
  fclose(sf);
  free(statfile);

  // Parse the line starting with "VmPeak"
  // Seems like it is always in kB
  // (reference: fs/proc/task_mmu.c)
  // actually in kiB i.e., 1024 bytes
  // since the last three characters are ' kb' we can skip them and parse in between
  size_t peakMemoryKB = 0;
  //  printf("peakline: '%s'\n", peakline);
  if(strlen(peakline) > 11)
  {
    peakline[strlen(peakline) -4] = '\0';

    //    printf("peakline: '%s'\n", peakline+7);
    peakMemoryKB = (size_t) atol(peakline+7);
  }

  free(peakline);
  return peakMemoryKB;
}
#endif

void fprint_peakMemory(FILE * fout)
{
  size_t pm = get_peakMemoryKB();

  if(fout == NULL) fout = stdout;
  fprintf(fout, "peakMemory: %zu kiB\n", pm);

  return;
}


float getErrorX(const float * restrict y, const float * restrict g, const int M, const int N, const int P, const int wM, const int wN, const int wP)
{
  /* Same as getError with the difference that G is expanded to MxNxP */
  assert(wM>=M);
  double e = 0;
  for(int c = 0; c<P; c++)
  {
    for(int b = 0; b<N; b++)
    {
      for(int a = 0; a<M; a++)
      {
        double yval = y[a + b*wM + c*wM*wN];
        double gval = g[a + b*wM + c*wM*wN];
        e+=pow(yval-gval, 2);
      }
    }
  }
  e/=(M*N*P);
  return (float) e;
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


float iter(float * xp, // Output, f_(t+1)
    fftwf_complex * restrict cK, // fft(psf)
    fftwf_complex * restrict cKr, // = NULL
    float * restrict f, // Current guess
    float * restrict W, // Weights
    float * restrict G, // expanded input image
    const int wM, const int wN, const int wP, // expanded size
    const int M, const int N, const int P, // input image size
    const dw_opts * s)
{
  const size_t wMNP = wM*wN*wP;

  fftwf_complex * F = fft(f, wM, wN, wP);

  float * y = fft_convolve_cc_f2(cK, F, wM, wN, wP);

  float error = getErrorX(y, G, M, N, P, wM, wN, wP);

  for(size_t kk = 0; kk<wMNP; kk++)
  {
    y[kk] < 1e-6 ? y[kk] = 1e-6 : 0;
    y[kk] = G[kk]/y[kk];
  }
  
  fftwf_complex * F_sn = fft(y, wM, wN, wP); 
  fftwf_free(y);

  //  float * x = fft_convolve_cc(cKr, F_sn, wM, wN, wP); // xp could be used for x
float * x = fft_convolve_cc_conj_f2(cK, F_sn, wM, wN, wP); 

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
   *
   * Possibly more stable to use the mean of the input image rather than 1
   */

  assert(wM > M); assert(wN > N); assert(wP > P);

  float * one = fim_zeros(wM*wN*wP);

  for(int cc = 0; cc < P; cc++) {
    for(int bb = 0; bb < N; bb++) {
      for(int aa = 0; aa < M; aa++) {
        one[aa + wM*bb + wM*wN*cc] = 1;
      }
    }
  }
  //  writetif("one.tif", one, wM, wN, wP);

  fftwf_complex * Fone = fft(one, wM, wN, wP);

  fftwf_free(one);
  return Fone;
}

void dw_usage(const int argc, char ** argv, const dw_opts * s)
{
  printf(" Usage:\n");
  printf("\t$ %s <options> image.tif psf.tif\n", argv[0]);
  //  printf("or\n");
  //  printf("\t$ %s --batch <options> image_dir psf_dir\n", argv[0]);
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
  printf(" --prefix str\n\t Set the prefix of the output files (default: '%s')\n", s->prefix);
  printf(" --overwrite\n\t Allows deconwolf to overwrite already existing output files\n");
  printf(" --relax F\n\t Multiply the central pixel of the PSF by F. (F>1 relaxation)\n");
  //  printf(" --batch\n\t Generate a batch file to deconvolve all images in the `image_dir`\n");
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

float * deconvolve_w(float * restrict im, const int M, const int N, const int P,
    const float * restrict psf, const int pM, const int pN, const int pP,
    dw_opts * s)
{
  if(s->verbosity > 0)
  {
    printf("Deconvolving\n");
  }

  /*Deconvolve im using psf */
  const int nIter = s->nIter;

  if(fim_maxAtOrigo(psf, pM, pN, pP) == 0)
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
  fprintf(s->log, "\timage: [%dx%dx%d]\npsf: [%dx%dx%d]\njob: [%dx%dx%d] (%zu voxels)\n",
      M, N, P, pM, pN, pP, wM, wN, wP, wMNP);
  fflush(s->log);

  fft_train(wM, wN, wP, 
      s->verbosity, s->nThreads);

  // cK : "full size" fft of the PSF
  float * Z = fftwf_malloc(wMNP*sizeof(float));
  memset(Z, 0, wMNP*sizeof(float));
  fim_insert(Z, wM, wN, wP, psf, pM, pN, pP);
  fim_circshift(Z, wM, wN, wP, -(pM-1)/2, -(pN-1)/2, -(pP-1)/2);
  fftwf_complex * cK = fft(Z, wM, wN, wP);
  //fim_tiff_write("Z.tif", Z, wM, wN, wP);

  fftwf_free(Z);

  fftwf_complex * cKr = NULL;

  /* <-- This isn't needed ...
     float * Zr = fftwf_malloc(wMNP*sizeof(float));
     memset(Zr, 0, wMNP*sizeof(float));
     float * psf_flipped = malloc(wMNP*sizeof(float));
     memset(psf_flipped, 0, sizeof(float)*wMNP);
     fim_flipall(psf_flipped, psf, pM, pN, pP);
     fim_insert(Zr, wM, wN, wP, psf_flipped, pM, pN, pP);
     free(psf_flipped);
     fim_circshift(Zr, wM, wN, wP, -(pM-1)/2, -(pN-1)/2, -(pP-1)/2);
     cKr = fft(Zr, wM, wN, wP); 
  // Possibly not needed due to f(-x) = ifft(conf(fft(f)))
  fftwf_free(Zr);
  --> */

  //printf("initial guess\n"); fflush(stdout);
  fftwf_complex * F_one = initial_guess(M, N, P, wM, wN, wP);
  float * P1 = fft_convolve_cc_conj_f2(cK, F_one, wM, wN, wP); // can't replace this one with cK!
  //printf("P1\n");
  //fim_stats(P1, pM*pN*pP);
  //  writetif("P1.tif", P1, wM, wN, wP);

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
  fim_insert(G, wM, wN, wP, im, M, N, P);
  //  writetif("G.tif", G, wM, wN, wP);

  if(s->method == DW_METHOD_RL)
  {
    //    float * x = deconvolve_rl(G, );
    //   free(G);
    //   return x;
  }

  float sumg = fim_sum(im, M*N*P);
  fftwf_free(im);

  float alpha = 0;
  float * f = fim_constant(wMNP, sumg/wMNP);
  float * y = fim_copy(f, wMNP);

  float * x1 = fim_copy(f, wMNP);
  float * x2 = fim_copy(f, wMNP);
  float * x = x1;
  float * xp = x2;
  float * xm = x2; // fim_copy(f, wMNP);

  float * gm = fim_zeros(wMNP);
  float * g = fim_zeros(wMNP);

  int it = 0;
  while(it<nIter)
  {

    if(s->iterdump){
      float * temp = fim_subregion(x, wM, wN, wP, M, N, P);
      char * tempname = malloc(100*sizeof(char));
      sprintf(tempname, "x_%03d.tif", it);
      printf("Writing current guess to %s\n", tempname);
      fim_tiff_write(tempname, temp, M, N, P);
      free(temp);
    }

    for(size_t kk = 0; kk<wMNP; kk++)
    { 
      y[kk] = x[kk] + alpha*(x[kk]-xm[kk]);
      y[kk] < 0 ? y[kk] = 0 : 0;
    }

    xp = xm;
    double err = iter(xp, // xp is updated to the next Guess
        cK, cKr, // FFT of PSF
        y, // Current guess
        W, // Weights (to handle boundaries)
        G, // expanded original image
        wM, wN, wP, // Expanded size
        M, N, P, // Original size
        s);
    if(s->verbosity > 0){
      if(s->verbosity> 1 || it+1 == nIter) {
        {printf("Iteration %d/%d, error=%e\n", it+1, nIter, err);}
      } else 
        if(it % 5 == 0)  {printf("Iteration %d/%d, error=%e\n", it+1, nIter, err);}
    }

    if(s->log != NULL)
    { fprintf(s->log, "Iteration %d/%d, error=%e\n", it+1, nIter, err);}

    float * swap = g;
    g = gm; gm = swap;

    fim_minus(g, xp, y, wMNP); // g = xp - y

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
  float * out = fim_subregion(x, wM, wN, wP, M, N, P);
  fftwf_free(f);
  fftwf_free(x1);
  fftwf_free(x2);
  fftwf_free(g);
  fftwf_free(gm);
  fftwf_free(cK);
  if(cKr != NULL)
  { fftwf_free(cKr); }
  fftwf_free(y);
  return out;
}

float * psf_autocrop(float * psf, int * pM, int * pN, int * pP,  // psf and size
    int M, int N, int P, // image size
    dw_opts * s)
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
    fprintf(s->log, "WARNING: The PSF has only %d slices, %d would be better.\n", p, popt);
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
    if(s->verbosity > 2)
    {
      printf("! %d %d : %d %d : %d %d\n", m0, m1, n0, n1, p0, p1);
    }
    float * psf_cropped = fim_get_cuboid(psf, m, n, p,
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


float * deconvolve_tiles(float * restrict im, int M, int N, int P,
    const float * restrict psf, const int pM, const int pN, const int pP,
    dw_opts * s)
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
    float * tpsf = fim_copy(psf, pM*pN*pP);
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

    tpsf = psf_autocrop(tpsf, &tpM, &tpN, &tpP, 
        tileM, tileN, tileP, s);

    fim_normalize_max1(tpsf, tpM, tpN, tpP);
    // Note: deconvolve_w will free it's first argument
    float * dw_im_tile = deconvolve_w(im_tile, tileM, tileN, tileP, // input image and size
        tpsf, tpM, tpN, tpP, // psf and size
        s);

    tiling_put_tile(T, tt, V, dw_im_tile);
    free(dw_im_tile);
    free(tpsf);
  }
  fftwf_free(im);
  tiling_free(T);
  free(T);
  return V;
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
    fim_flipall(V, A, M, N, P);
  toc(fim_flipall)

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
    float * S1 = fim_subregion(V, M, N, P, M-1, N-1, P-1);
  toc(fim_subregion)

    tic
    float * S2 = fim_subregion_ref(V, M, N, P, M-1, N-1, P-1);
  toc(fim_subregion_ref)
    printf("S1 - S2 = %f\n", getError(S1, S1, M-1, N-1, P-1, M-1, N-1, P-1));
  free(S1);
  free(S2);

  // ---
  tic
    fim_insert(V, M, N, P, A, M-1, N-1, P-1);
  toc(fim_subregion)      

    tic
    fim_insert_ref(V, M, N, P, A, M-1, N-1, P-1);
  toc(fim_subregion_ref)

    // ---


    ((float volatile *)V)[0] = V[0];
  printf("V[0] = %f\n", V[0]);
  free(A);
  free(V);
}

void dw_unittests()
{
  fprint_peakMemory(NULL);
  timings();

  fim_ut();
  fim_tiff_ut();
  fft_ut();
  printf("done\n");
}

void show_time(FILE * f)
{
  f == NULL ? f = stdout : 0;
  time_t now = time(NULL);
  char * tstring = ctime(&now);
  fprintf(f, "%s\n", tstring);
}

float * psf_makeOdd(float * psf, int * pM, int * pN, int *pP)
{
  // Expand the psf so that it had odd dimensions it if doesn't already have that
  int m = pM[0];
  int n = pN[0];
  int p = pP[0];
  int reshape = 0;
  if(m % 2 == 0)
  { m++; reshape = 1;}
  if(n % 2 == 0)
  { n++; reshape = 1;}
  if(p % 2 == 0)
  { p++; reshape = 1;}

  if(reshape == 0)
  {  return psf; }
  // printf("%d %d %d -> %d %d %d\n", pM[0], pN[0], pP[0], m, n, p);
  float * psf2 = fim_zeros(m*n*p);
  fim_insert(psf2, m, n, p, psf, pM[0], pN[0], pP[0]);
  free(psf);
  pM[0] = m; 
  pN[0] = n; 
  pP[0] = p; 
  return psf2;
}


void dcw_init_log(dw_opts * s)
{
  s->log = fopen(s->logFile, "w");
  assert(s->log != NULL);
  show_time(s->log);
  dw_opts_fprint(s->log, s); 
  dw_fprint_info(s->log);
}

void dcw_close_log(dw_opts * s)
{
  fprint_peakMemory(s->log);
  show_time(s->log);
  fclose(s->log);
}


int dw_run(dw_opts * s)
{
  struct timespec tstart, tend;
  clock_gettime(CLOCK_REALTIME, &tstart);

  dcw_init_log(s);

  if(s->verbosity > 0) 
  {
    dw_opts_fprint(NULL, s); 
    printf("\n");
  }

  s->verbosity > 1 ? dw_fprint_info(NULL) : 0;

  if(s->verbosity > 0)
  {
    printf("Reading %s\n", s->imFile);
  }
  int M = 0, N = 0, P = 0;
  float * im = fim_tiff_read(s->imFile, &M, &N, &P, s->verbosity);
  if(s->verbosity > 2)
  {
    printf("Image size: [%d x %d x %d]\n", M, N, P);
  }
  if(im == NULL)
  {
    printf("Failed to open %s\n", s->imFile);
    exit(1);
  }

  // writetif("im.tif", im, M, N, P);

  int pM = 0, pN = 0, pP = 0;
  float * psf = NULL;
  if(1){
    if(s->verbosity > 0)
    {
      printf("Reading %s\n", s->psfFile);
    }
    psf = fim_tiff_read(s->psfFile, &pM, &pN, &pP, s->verbosity);
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

  psf = psf_makeOdd(psf, &pM, &pN, &pP);

  // Possibly the PSF will be cropped even more per tile later on
  psf = psf_autocrop(psf, &pM, &pN, &pP, 
      M, N, P, s);

  if(s->relax > 0)
  {
    // Note: only works with odd sized PSF
    fprintf(s->log, "Relaxing the PSF by %f\n", s->relax);
    if(s->verbosity > 0)
    {    
    printf("Relaxing the PSF\n");
    }
    size_t mid = (pM-1)/2 + (pN-1)/2*pM + (pP-1)/2*pM*pN;
    psf[mid] *= s->relax;
    fim_normalize_max1(psf, pM, pN, pP);
  }

  myfftw_start(s->nThreads);
  float * out = NULL;
  if(s->tiling_maxSize < 0)
  {

    fim_normalize_max1(psf, pM, pN, pP);
    // Will free im
    out = deconvolve_w(im, M, N, P, // input image and size
        psf, pM, pN, pP, // psf and size
        s);// settings
  } else {
    // Will free im
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
    fim_tiff_write(s->outFile, out, M, N, P);
    exitstatus = 0;
  }


  if(s->verbosity > 0)
  {
    printf("Finalizing\n"); fflush(stdout);
  }
  free(psf);
  free(out);
  myfftw_stop();
  if(s->verbosity > 1) fprint_peakMemory(NULL);
  clock_gettime(CLOCK_REALTIME, &tend);
  fprintf(s->log, "Took %f s\n", timespec_diff(&tend, &tstart));
  dcw_close_log(s);
  dw_opts_free(&s);

  return(exitstatus);
}


