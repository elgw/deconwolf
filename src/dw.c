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
#ifdef _OPENMP // turned on with -fopenmp
#include <omp.h>
#endif 
#include "fft.h"
#include "dw.h"
#include "tiling.h"
#include "fim.h"
#include "fim_tiff.h"


#define tictoc struct timespec tictoc_start, tictoc_end;
#define tic clock_gettime(CLOCK_REALTIME, &tictoc_start);
#define toc(X) clock_gettime(CLOCK_REALTIME, &tictoc_end); printf(#X); printf(" %f s\n", timespec_diff(&tictoc_end, &tictoc_start)); fflush(stdout);

/* This is what fftw_malloc returns
 * http://www.fftw.org/fftw3_doc/SIMD-alignment-and-fftw_005fmalloc.html#SIMD-alignment-and-fftw_005fmalloc
 */
typedef float afloat __attribute__ ((__aligned__(16)));

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
  s->xycropfactor = 0.001;
  s->commandline = NULL;
  s->onetile = 0;
  s->borderQuality = 2;
  return s;
}

void nullfree(void * p)
{
  if(p != NULL)
  {
    free(p);
  }
}


static int64_t int64_t_max(int64_t a, int64_t b)
{
  if( a > b)
    return a;
  return b;
}

void dw_opts_free(dw_opts ** sp)
{
  dw_opts * s = sp[0];
  nullfree(s->imFile);
  nullfree(s->psfFile);
  nullfree(s->outFile);
  nullfree(s->logFile);
  nullfree(s->prefix);
  nullfree(s->commandline);
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
  fprintf(f, "XY crop factor: %f\n", s->xycropfactor);
  if(s->relax > 0)
  {
    fprintf(f, "PSF relaxation: %f\n", s->relax);
  }
  fprintf(f, "Border Quality: ");
  switch(s->borderQuality){
    case 0:
      fprintf(f, "0 = worst\n");
      break;
    case 1:
      fprintf(f, "1 = low\n");
      break;
    case 2:
      fprintf(f, "2 = normal\n");
        break;
    default:
        ;
  }
  if(s->onetile == 1)
  {
    fprintf(f, "DEBUG OPTION: ONETILE = TRUE (only first tile will be deconvolved)\n");
  }
  fprintf(f, "\n");
}

static double timespec_diff(struct timespec* end, struct timespec * start)
{
  double elapsed = (end->tv_sec - start->tv_sec);
  elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
  return elapsed;
}

int file_exist(char * fname)
{
  if( access( fname, F_OK ) != -1 ) {
    return 1; // File exist
  } else {
    return 0;
  }
}

void dw_fprint_info(FILE * f, dw_opts * s)
{
  f == NULL ? f = stdout : 0;
  fprintf(f, "deconwolf: '%s' PID: %d\n", deconwolf_version, (int) getpid());
  char cwd[1024];
  if (getcwd(cwd, sizeof(cwd)) != NULL) {
    fprintf(f, "PWD: %s\n", cwd);
  } 

  if(s->commandline != NULL)
  {
    fprintf(f, "CMD: %s\n", s->commandline);
  }

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
#ifdef _OPENMP
  fprintf(f, "OpenMP: YES\n");
#endif

  fprintf(f, "\n");
  fflush(f);
  return;
}

void deconwolf_batch(dw_opts * s)
{
  // Better to do in Lua/Python/bash/...
  printf("Not implemented !\n");
}

static void getCmdLine(int argc, char ** argv, dw_opts * s)
{
  // Copy the command line to s->commandline
  int lcmd=0;
  for(int kk = 0; kk<argc; kk++)
  {
    lcmd += strlen(argv[kk]);
  }
  lcmd += argc+2;
  s->commandline = malloc(lcmd);
  int pos = 0;
  for(int kk = 0; kk<argc; kk++)
  {
    sprintf(s->commandline+pos, "%s ", argv[kk]);
    pos += strlen(argv[kk])+1;
  }
  s->commandline[pos-1] = '\0';
  //  printf("argc: %d cmd: '%s'\n", argc, s->commandline);
}

void dw_argparsing(int argc, char ** argv, dw_opts * s)
{

  getCmdLine(argc, argv, s);

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
    { "xyfactor",     required_argument, NULL,   'x' },
    { "onetile",      no_argument,       NULL,   'T' },
    { "bq",           required_argument, NULL,   'B' },
    { NULL,           0,                 NULL,   0   }
  };

  int ch;
  while((ch = getopt_long(argc, argv, "Bvho:n:c:p:s:p:T", longopts, NULL)) != -1)
  {
    switch(ch) {
      case 'B':
        s->borderQuality = atoi(optarg);
        break;
      case 'v':
        dw_fprint_info(NULL, s);
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
      case 'x':
        s->xycropfactor = atof(optarg);
        if(s->xycropfactor > 1 || s->xycropfactor < 0)
        {
          printf("The crop factor in x and y has to be => 0 and < 1\n");
          exit(1);
        }
        break;
      case 'T':
        s->onetile = 1;
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

  s->imFile = realpath(argv[optind], 0);
  if(s->imFile == NULL)
  {
    printf("ERROR: Can't read %s\n", argv[optind]);
    exit(1);
  }
  s->psfFile = realpath(argv[++optind], 0);
  if(s->psfFile == NULL)
  {
    printf("ERROR: Can't read %s\n", argv[optind]);
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


void fsetzeros(const char * fname, size_t N)
  /* Create a new file and fill it with N bytes of zeros
  */
{
  size_t bsize = 1024*1024;
  char * buffer = malloc(bsize);
  memset(buffer, 0, bsize);
  FILE * fid = fopen(fname, "wb");
  size_t written = 0;
  while(written + bsize < N)
  {
    fwrite(buffer, bsize, 1, fid);
    written += bsize;
  }
  //  printf("\r %zu / %zu \n", written, N); fflush(stdout);
  
  fwrite(buffer, N-written, 1, fid);
  fclose(fid);
  free(buffer);
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


float getErrorX(const float * restrict y, const float * restrict g, const int64_t M, const int64_t N, const int64_t P, const int64_t wM, const int64_t wN, const int64_t wP)
{
  /* Same as getError with the difference that G is expanded to MxNxP */
  assert(wM>=M);
  double e = 0;
  for(size_t c = 0; c<P; c++)
  {
    for(size_t b = 0; b<N; b++)
    {
      for(size_t a = 0; a<M; a++)
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


float getError(const afloat * restrict y, const afloat * restrict g, 
    const int64_t M, const int64_t N, const int64_t P, 
    const int64_t wM, const int64_t wN, const int64_t wP)
{
  double e = 0;
#pragma omp parallel for reduction(+: e)
  for(int64_t c = 0; c<P; c++)
  {
    for(int64_t b = 0; b<N; b++)
    {
      for(int64_t a = 0; a<M; a++)
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

float getError_ref(float * y, float * g, int64_t M, int64_t N, int64_t P, int64_t wM, int64_t wN, int64_t wP)
{
  double e = 0;
  for(int64_t a = 0; a<M; a++)
  {
    for(int64_t b = 0; b<N; b++)
    {
      for(int64_t c = 0; c<P; c++)
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


float iter(
    afloat ** xp, // Output, f_(t+1)
    const float * restrict im, // Input image
    fftwf_complex * restrict cK, // fft(psf)
    fftwf_complex * restrict cKr, // = NULL
    afloat * restrict f, // Current guess
    afloat * restrict W, // Weights
    const int64_t wM, const int64_t wN, const int64_t wP, // expanded size
    const int64_t M, const int64_t N, const int64_t P, // input image size
    const dw_opts * s)
{
  // We could reduce memory even further by using 
  // the allocation for xp
  const size_t wMNP = wM*wN*wP;

  fftwf_complex * F = fft(f, wM, wN, wP);
  afloat * y = fft_convolve_cc_f2(cK, F, wM, wN, wP); // F is freed
  float error = getError(y, im, M, N, P, wM, wN, wP);


#pragma omp parallel for
  for(size_t cc =0; cc<wP; cc++){
    for(size_t bb = 0; bb<wN; bb++){
      for(size_t aa = 0; aa<wM; aa++){
        size_t yidx = aa + bb*wM + cc*wM*wN;
        size_t imidx = aa + bb*M + cc*M*N;
        if(aa<M && bb<N && cc<P)
        {
          y[yidx] < 1e-6 ? y[yidx] = 1e-6 : 0;
          y[yidx]=im[imidx]/y[yidx];
        } else {
          y[yidx]=0;
        }
      }
    }
  }

  fftwf_complex * F_sn = fft(y, wM, wN, wP); 
  fftwf_free(y); 

  afloat * x = fft_convolve_cc_conj_f2(cK, F_sn, wM, wN, wP); 

#pragma omp parallel for
  for(size_t cc = 0; cc<wMNP; cc++)
  { x[cc] *= f[cc]*W[cc]; }

  xp[0] = x;
  return error;
}

fftwf_complex * initial_guess(const int64_t M, const int64_t N, const int64_t P, 
    const int64_t wM, const int64_t wN, const int64_t wP)
{
  /* Create initial guess: the fft of an image that is 1 in MNP and 0 outside
   * M, N, P is the dimension of the microscopic image
   *
   * Possibly more stable to use the mean of the input image rather than 1
   */

  assert(wM >= M); assert(wN >= N); assert(wP >= P);

  afloat * one = fim_zeros(wM*wN*wP);

#pragma omp parallel for
  for(int64_t cc = 0; cc < P; cc++) {
    for(int64_t bb = 0; bb < N; bb++) {
      for(int64_t aa = 0; aa < M; aa++) {
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
  printf(" --xyfactor F\n\t Discard outer planes of the PSF with sum < F of the central\n");
  printf(" --bq Q\n\t Set border quality to 0 'worst', 1 'bad', or 2 'normal' which is default\n");
  //  printf(" --batch\n\t Generate a batch file to deconvolve all images in the `image_dir`\n");
  printf("\n");
}


float update_alpha(const afloat * restrict g, const afloat * restrict gm, const size_t wMNP)
{
  double ggm = 0;
  double gmgm = 0;

#pragma omp parallel for reduction(+: ggm, gmgm)
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

float * deconvolve_w(afloat * restrict im, const int64_t M, const int64_t N, const int64_t P,
    const afloat * restrict psf, const int64_t pM, const int64_t pN, const int64_t pP,
    dw_opts * s)
{
  if(s->verbosity > 1)
  {
    printf("Deconvolving\n");
  }

  if(s->nIter == 0)
  {  
    return fim_copy(im, M*N*P);
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
  // Default border quality
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



  size_t wMNP = wM*wN*wP;




  if(s->verbosity > 0)
  { printf("image: [%" PRId64 "x%" PRId64 "x%" PRId64 "], psf: [%" PRId64 "x%" PRId64 "x%" PRId64 "], job: [%" PRId64 "x%" PRId64 "x%" PRId64 "] (%zu voxels)\n",
      M, N, P, pM, pN, pP, wM, wN, wP, wMNP);
  }
  if(s->verbosity > 1)
  {
    printf("Estimated peak memory usage: %.1f GB\n", wMNP*35.0/1e9);
  }
  fprintf(s->log, "image: [%" PRId64 "x%" PRId64 "x%" PRId64 "]\npsf: [%" PRId64 "x%" PRId64 "x%" PRId64 "]\njob: [%" PRId64 "x%" PRId64 "x%" PRId64 "] (%zu voxels)\n",
      M, N, P, pM, pN, pP, wM, wN, wP, wMNP);
  fflush(s->log);

  fft_train(wM, wN, wP, 
      s->verbosity, s->nThreads);

  // cK : "full size" fft of the PSF
  afloat * Z = fftwf_malloc(wMNP*sizeof(float));
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
  afloat * P1 = fft_convolve_cc_conj_f2(cK, F_one, wM, wN, wP); // can't replace this one with cK!
  //printf("P1\n");
  //fim_stats(P1, pM*pN*pP);
  //  writetif("P1.tif", P1, wM, wN, wP);

  float sigma = 0.001;
#pragma omp parallel for 
  for(size_t kk = 0; kk<wMNP; kk++)
  {
    if(P1[kk] < sigma)
    { P1[kk] = 0; } else { P1[kk] = 1/P1[kk]; }
  }
  afloat * W = P1;
  //writetif("W.tif", W, wM, wN, wP);

  // Original image -- expanded
  //  float * G = fftwf_malloc(wMNP*sizeof(float));
  // memset(G, 0, wMNP*sizeof(float));
  // fim_insert(G, wM, wN, wP, im, M, N, P);
  //  writetif("G.tif", G, wM, wN, wP);

  if(s->method == DW_METHOD_RL)
  {
    //    float * x = deconvolve_rl(G, );
    //   free(G);
    //   return x;
  }
  float sumg = fim_sum(im, M*N*P);

  float alpha = 0;
  afloat * f = fim_constant(wMNP, sumg/wMNP);
  afloat * y = fim_copy(f, wMNP);

  afloat * x1 = fim_copy(f, wMNP);
  afloat * x2 = fim_copy(f, wMNP);
  afloat * x = x1;
  afloat * xp = x2;
  afloat * xm = x2; // fim_copy(f, wMNP);

  afloat * gm = fim_zeros(wMNP);
  afloat * g = fim_zeros(wMNP);

  if(s->verbosity > 0)
  {
    printf("Iterating ..."); fflush(stdout);
  }

  int it = 0;
  while(it<nIter)
  {

    if(s->iterdump){
      afloat * temp = fim_subregion(x, wM, wN, wP, M, N, P);
      char * tempname = malloc(100*sizeof(char));
      sprintf(tempname, "x_%03d.tif", it);
      printf("Writing current guess to %s\n", tempname);
      fim_tiff_write(tempname, temp, M, N, P);
      free(temp);
    }

#pragma omp parallel for
    for(size_t kk = 0; kk<wMNP; kk++)
    { 
      y[kk] = x[kk] + alpha*(x[kk]-xm[kk]);
      y[kk] < 0 ? y[kk] = 0 : 0;
    }

    xp = xm;
    free(xp);
    double err = iter(
        &xp, // xp is updated to the next guess
        im,
        cK, NULL, // FFT of PSF
        y, // Current guess
        W, // Weights (to handle boundaries)
        wM, wN, wP, // Expanded size
        M, N, P, // Original size
        s);

    if(s->verbosity > 0){
      printf("\rIteration %3d/%3d, error=%e", it+1, nIter, err);
      fflush(stdout);
    }

    if(s->log != NULL)
    { fprintf(s->log, "Iteration %d/%d, error=%e\n", it+1, nIter, err); fflush(s->log); }

    afloat * swap = g;
    g = gm; gm = swap;

    fim_minus(g, xp, y, wMNP); // g = xp - y

    if(it > 0) {
      alpha = update_alpha(g, gm, wMNP);
    }

    xm = x;
    x = xp;
    xp = NULL;

    it++;
  } // End of main loop

  if(s->verbosity > 0) {
    printf("\n");
  }

  fftwf_free(W); // is P1
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
  if(cKr != NULL)
  { fftwf_free(cKr); }
  fftwf_free(y);
  return out;
}

float * psf_autocrop_centerZ(float * psf, int64_t * pM, int64_t * pN, int64_t * pP,  // psf and size
    dw_opts * s)
{

  const int64_t m = pM[0];
  const int64_t n = pN[0];
  const int64_t p = pP[0];

  const int64_t midm = (m-1)/2;
  const int64_t midn = (n-1)/2;
  const int64_t midp = (p-1)/2;

  //  printf("m: %d, n:%d, p:%d\n", m, n, p);
  //  printf("midm: %d, midn: %d, midp: %d\n", midm, midn, midp);

  float maxvalue = -INFINITY;
  int64_t maxp = -1;

  for(int64_t pp = 0; pp<p; pp++)
  {
    size_t idx = midm + midn*m + pp*m*n;
    if(psf[idx] > maxvalue)
    {
      maxp = pp;
      maxvalue = psf[idx];
    }
  }

  if(maxp == midp)
  {
    if(s->verbosity > 2)
    {
      printf("PSF is Z-centered :)\n");
    }
    return psf;
  }


  int64_t m0 = 0, m1 = m-1;
  int64_t n0 = 0, n1 = n-1;
  int64_t p0 = maxp, p1 = maxp;

  while(p0 > 1 && p1+2 < p)
  {
    p0--; p1++;
  }
  if(s->verbosity > 2)
  {
    printf("PSF has %" PRId64 " slices\n", p);
    printf("brighest at plane %" PRId64 "\n", maxp);
    printf("Selecting Z-planes: %" PRId64 " -- %" PRId64 "\n", p0, p1);
  }

  fprintf(s->log, "Selecting Z-planes %" PRId64 " -- %" PRId64 "\n", p0, p1);

  float * psf_cropped = fim_get_cuboid(psf, m, n, p,
      m0, m1, n0, n1, p0, p1);
  free(psf);
  pP[0] = p1-p0+1;
  return psf_cropped;

}

float * psf_autocrop_byImage(float * psf, int64_t * pM, int64_t * pN, int64_t * pP,  // psf and size
    int64_t M, int64_t N, int64_t P, // image size
    dw_opts * s)
{

  const int64_t m = pM[0];
  const int64_t n = pN[0];
  const int64_t p = pP[0];

  if((p % 2) == 0)
  {
    printf("Error: The PSF should have odd number of slices\n");
    exit(1);
  }

  // Optimal size
  int64_t mopt = (M-1)*2 + 1;
  int64_t nopt = (N-1)*2 + 1;
  int64_t popt = (P-1)*2 + 1;

  if(p < popt)
  {
    fprintf(s->log, "WARNING: The PSF has only %" PRId64 " slices, %" PRId64 " would be better.\n", p, popt);
    return psf;
  }

  if(m > mopt || n > nopt || p > popt)
  { 
    int64_t m0 = 0, m1 = m-1;
    int64_t n0 = 0, n1 = n-1;
    int64_t p0 = 0, p1 = p-1;
    if(m > mopt)
    {
      m0 = (m-mopt)/2;
      m1 = m1-(m-mopt)/2;
    }
    if(n > nopt)
    {
      n0 = (n-nopt)/2;
      n1 = n1-(n-nopt)/2;
    }
    if(p > popt)
    {
      p0 = (p-popt)/2;
      p1 = p1-(p-popt)/2;
    }
    if(s->verbosity > 2)
    {
      printf("! %" PRId64 " %" PRId64 " : %" PRId64 " %" PRId64 " : %" PRId64 " %" PRId64 "\n", m0, m1, n0, n1, p0, p1);
    }
    float * psf_cropped = fim_get_cuboid(psf, m, n, p,
        m0, m1, n0, n1, p0, p1);
    free(psf);

    pM[0] = m1-m0+1;
    pN[0] = n1-n0+1;
    pP[0] = p1-p0+1;

    if(s->verbosity > 0)
    {
      fprintf(stdout, "PSF crop [%" PRId64 " x %" PRId64 " x %" PRId64 "] -> [%" PRId64 " x %" PRId64 " x %" PRId64 "]\n", 
          m, n, p, pM[0], pN[0], pP[0]);
    }
    fprintf(s->log, "PSF crop [%" PRId64 " x %" PRId64 " x %" PRId64 "] -> [%" PRId64 " x %" PRId64 " x %" PRId64 "]\n", 
        m, n, p, pM[0], pN[0], pP[0]);

    return psf_cropped;
  } else {
    return psf;
  }
}

float * psf_autocrop_XY(float * psf, int64_t * pM, int64_t * pN, int64_t * pP,  // psf and size
    int64_t M, int64_t N, int64_t P, // image size
    dw_opts * s)
{
  // Find the y-z plane with the largest sum
  float maxsum = 0;
  for(int64_t xx = 0; xx<pM[0]; xx++)
  {
    float sum = 0;
    for(int64_t yy = 0; yy<pN[0]; yy++)
    {
      for(int64_t zz = 0; zz<pP[0]; zz++)
      {
        sum += psf[xx + yy*pM[0] + zz*pM[0]*pN[0]];
      }
    }
    sum > maxsum ? maxsum = sum : 0;
  }

  //  printf("X maxsum %f\n", maxsum);

  int64_t first=-1;
  float sum = 0;

  while(sum < s->xycropfactor * maxsum)
  {
    first++;
    sum = 0;
    int64_t xx = first;
    for(int64_t yy = 0; yy<pN[0]; yy++)
    {
      for(int64_t zz = 0; zz<pP[0]; zz++)
      {
        sum += psf[xx + yy*pM[0] + zz*pM[0]*pN[0]];
      }
    }
  }

  if(first == 0)
  {
    if(s->verbosity > 1)
    {
      printf("PSF X-crop: Not cropping\n");
    }
    return psf;
  }

  float * p = fim_get_cuboid(psf, pM[0], pN[0], pP[0],
      first, pM[0] - first -1, 
      first, pN[0] - first -1, 
      0, pP[0]-1);
  pM[0] -= 2*first;
  pN[0] -= 2*first;

  if(s->verbosity > 0)
  {
    printf("PSF X-crop: Removing %" PRId64 " planes in XY\n", first);
  }

  free(psf);
  return p;
}

float * psf_autocrop(float * psf, int64_t * pM, int64_t * pN, int64_t * pP,  // psf and size
    int64_t M, int64_t N, int64_t P, // image size
    dw_opts * s)
{
  float * p = psf;
  p = psf_autocrop_centerZ(p, pM, pN, pP, s);
  assert(pM[0] > 0);
  // Crop the PSF if it is larger than necessary
  p = psf_autocrop_byImage(p, pM, pN, pP, M, N, P, s);
  assert(pM[0] > 0);
  // Crop the PSF by removing outer planes that has very little information
  p = psf_autocrop_XY(p, pM, pN, pP, M, N, P, s);
  assert(pM[0] > 0);
  assert(p != NULL);
  return p;
}


int deconvolve_tiles(const int64_t M, const int64_t N, const int64_t P,
    const float * restrict psf, const int64_t pM, const int64_t pN, const int64_t pP,
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
    printf("-> Divided the [%" PRId64 " x %" PRId64 " x %" PRId64 "] image into %d tiles\n", M, N, P, T->nTiles);
  }

  /* Output image initialize as zeros
   * will be updated block by block
   */
  char * tfile = malloc(strlen(s->outFile)+10);
  sprintf(tfile, "%s.raw", s->outFile);

  if(s->verbosity > 0)
  {
    printf("Initializing %s to 0\n", tfile); fflush(stdout);
  }
    fsetzeros(tfile, (size_t) M* (size_t) N* (size_t) P*sizeof(float));

      char * imFileRaw = malloc(strlen(s->imFile) + 10);
    sprintf(imFileRaw, "%s.raw", s->imFile);

      if(s->verbosity > 0)
      {
        printf("Dumping %s to %s (for quicker io)\n", s->imFile, imFileRaw);
      }
    
    fim_tiff_to_raw(s->imFile, imFileRaw);
    if(0){
    printf("Writing to imdump.tif\n");
    fim_tiff_from_raw("imdump.tif", M, N, P, imFileRaw);
    }

    //fim_tiff_write_zeros(s->outFile, M, N, P);
    if(s->verbosity > 0)
    {
      printf("\n"); fflush(stdout);
    }

  int nTiles = T->nTiles;
  if(s->onetile == 1)
  {
    nTiles = 1;
    fprintf(s->log, "DEBUG: only the first tile to be deconvolved\n");
    fprintf(stdout, "DEBUG: only the first tile to be deconvolved\n");
  }

  for(int tt = 0; tt < nTiles; tt++)
  {
    // Temporal copy of the PSF that might be cropped to fit the tile
    float * tpsf = fim_copy(psf, pM*pN*pP);
    int64_t tpM = pM, tpN = pN, tpP = pP;

    if(s->verbosity > 0)
    {
      printf("-> Processing tile %d / %d\n", tt+1, T->nTiles);
      fprintf(s->log, "-> Processing tile %d / %d\n", tt+1, T->nTiles);
    }

//    tictoc
 //   tic
    //float * im_tile = tiling_get_tile_tiff(T, tt, s->imFile);
    float * im_tile = tiling_get_tile_raw(T, tt, imFileRaw);
//    toc(tiling_get_tile_tiff)

    int64_t tileM = T->tiles[tt]->xsize[0];
    int64_t tileN = T->tiles[tt]->xsize[1];
    int64_t tileP = T->tiles[tt]->xsize[2];
    
if(0)
{
  printf("writing to tiledump.tif\n");
    fim_tiff_write("tiledump.tif", im_tile, tileM, tileN, tileP);
    getchar();
}

    fim_normalize_sum1(tpsf, tpM, tpN, tpP);

    tpsf = psf_autocrop(tpsf, &tpM, &tpN, &tpP, 
        tileM, tileN, tileP, s);

    fim_normalize_sum1(tpsf, tpM, tpN, tpP);

    float * dw_im_tile = deconvolve_w(im_tile, tileM, tileN, tileP, // input image and size
        tpsf, tpM, tpN, tpP, // psf and size
        s);
    free(im_tile);
    tiling_put_tile_raw(T, tt, tfile, dw_im_tile);
    free(dw_im_tile);
    free(tpsf);
  }
  tiling_free(T);
  free(T);

  if(s->verbosity > 2)
{
  printf("converting %s to %s\n", tfile, s->outFile);
}
  fim_tiff_from_raw(s->outFile, M, N, P, tfile);

  if(s->verbosity < 5)
  {
    remove(tfile);
  } else {
    printf("Keeping %s for inspection, remove manually\n", tfile);
  }
  free(tfile);

  remove(imFileRaw);
  free(imFileRaw);
  return 0;
}



void timings()
{
  printf("-> Timings\n");
  tictoc
    int64_t M = 1024, N = 1024, P = 50;
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


  tic
    for(size_t kk = 0; kk<M*N*P; kk++)
    {
      A[kk] = (float) rand()/(float) RAND_MAX;
    }
  toc(rand)

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

float * psf_makeOdd(float * psf, int64_t * pM, int64_t * pN, int64_t *pP)
{
  // Expand the psf so that it had odd dimensions it if doesn't already have that
  int64_t m = pM[0];
  int64_t n = pN[0];
  int64_t p = pP[0];
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
  dw_fprint_info(s->log, s);
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


#ifdef _OPENMP
  omp_set_num_threads(s->nThreads);
  omp_set_dynamic(1);
  //  printf("omp_get_max_threads: %d\n", omp_get_max_threads());
#endif

  if(s->verbosity > 1) 
  {
    dw_opts_fprint(NULL, s); 
    printf("\n");
  }

  s->verbosity > 1 ? dw_fprint_info(NULL, s) : 0;



  int64_t M = 0, N = 0, P = 0;
  if(fim_tiff_get_size(s->imFile, &M, &N, &P))
  {
    printf("Failed to open %s\n", s->imFile);
    return -1;
  }

  if(s->verbosity > 0)
  {
    printf("Image dimensions: %" PRId64 " x %" PRId64 " x %" PRId64 "\n", M, N, P);
  }

  int tiling = 0;
  if(s->tiling_maxSize > 0 && (M > s->tiling_maxSize || N > s->tiling_maxSize))
  {
    tiling = 1;
  }


  float * im = NULL;

  if(tiling == 0)
  {
    if(s->verbosity > 0 )
    {
      printf("Reading %s\n", s->imFile);
    }

    im = fim_tiff_read(s->imFile, &M, &N, &P, s->verbosity);
    if(fim_min(im, M*N*P) < 0)
    {
      printf("min value of the image is %f, shifting to 0\n", fim_min(im, M*N*P));
      fim_set_min_to_zero(im, M*N*P);
      if(fim_max(im, M*N*P) < 1000)
      {
        fim_mult_scalar(im, M*N*P, 1000/fim_max(im, M*N*P));
      }
    }
    if(im == NULL)
    {
      printf("Failed to open %s\n", s->imFile);
      exit(1);
    }

  }

  // fim_tiff_write("identity.tif", im, M, N, P);

  int64_t pM = 0, pN = 0, pP = 0;
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
  psf = psf_autocrop_centerZ(psf, &pM, &pN, &pP, s);

  if(fim_maxAtOrigo(psf, pM, pN, pP) == 0)
  {
    printf("PSF is not centered!\n");
    return -1;
  }

  // Possibly the PSF will be cropped even more per tile later on

  fim_normalize_sum1(psf, pM, pN, pP);
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
    fim_normalize_sum1(psf, pM, pN, pP);
  }

  if(s->verbosity > 0)
  {
    printf("Output: %s(.log.txt)\n", s->outFile);
  }

  myfftw_start(s->nThreads);

  float * out = NULL;

  if(tiling)
  {
    deconvolve_tiles(M, N, P, psf, pM, pN, pP, // psf and size
        s);// settings
  } else {
    fim_normalize_sum1(psf, pM, pN, pP);
    out = deconvolve_w(im, M, N, P, // input image and size
        psf, pM, pN, pP, // psf and size
        s);// settings
  }

  if(tiling == 0)
  {
    fftwf_free(im);


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
    }
  }

  if(s->verbosity > 1)
  {
    printf("Finalizing\n"); fflush(stdout);
  }
  free(psf);
  if(out != NULL) free(out);
  myfftw_stop();
  if(s->verbosity > 1) fprint_peakMemory(NULL);
  clock_gettime(CLOCK_REALTIME, &tend);
  fprintf(s->log, "Took: %f s\n", timespec_diff(&tend, &tstart));
  dcw_close_log(s);
  dw_opts_free(&s);

  return 0;
}


