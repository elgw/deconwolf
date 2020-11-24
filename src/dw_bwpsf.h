#ifndef dw_bwpsf_h
#define dw_bwpsf_h

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <getopt.h>
#include <time.h>
#include "fim_tiff.h"
#include "dw_version.h"

// Inspiration:
//https://github.com/Biomedical-Imaging-Group/PSFGenerator/tree/master/src/psf/bornwolf
// gcc -Wall dw_bwpsf.c fim.c fim_tiff.c -lm -ltiff -lfftw3f -lpthread

typedef struct {
  int verbose;
  int overwrite;

  char * cmd; // Command line
  char * outFile;
  char * logFile;
  FILE * log;

  float lambda; // Emission maxima
  float NA; // Numerical aperture of lens
  float ni;
  float TOL;
  size_t K; // Number of iterations
  float resLateral;
  float resAxial;
  int Simpson; // Use this number of points for Simpson integration
  int fast_li; // Use Li's method for the integral
  int oversampling_R; // Oversampling in radial direction

  float * V;
  // shape of output image
  int M;
  int N;
  int P;
  // For parallel processing
  int thread;
  int nThreads;
} bw_conf;

void bw_conf_printf(FILE * out, bw_conf * conf);
void BW_slice(float * , float z, bw_conf * conf);
void bw_argparsing(int , char ** , bw_conf * conf);
void * BW_thread(void * data);

/* Generate the Born-Wolf PSF
 * allocates conf->V is already allocated
 */
void BW(bw_conf * conf);


void getCmdLine(int argc, char ** argv, bw_conf * s);


#endif
