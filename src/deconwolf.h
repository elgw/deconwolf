#ifndef deconwolf_h
#define deconwolf_h

#define deconwolf_version "alpha-0.002"

/* deconwolf
 *
 * fftw3 wisdom data is stored and loaded from
 * $home/.config/ 
 *
 * Internally column major indexing is used, i.e., the distance 
 * between elements is 1 for the first dimension and increases with 
 * each new dimension.
 *
 * Erik Wernersson, 2020
 */

#define DW_METHOD_W 0
#define DW_METHOD_RL 1
#define DW_METHOD_ID 2

typedef struct{
  int nThreads;
  int nIter;
  char * imFile;
  char * psfFile;
  char * outFile;
  char * logFile;
  char * prefix;
  FILE * log;
  int tiling_maxSize;
  int tiling_padding;
  int overwrite; // overwrite output tif file?
  int method;
  int verbosity;
  fftwf_plan fft_plan;
  fftwf_plan ifft_plan;
  int iterdump; // Dump each iteration to file ... 
} dw_opts;

dw_opts * dw_opts_new(void);
void dw_argparsing(int argc, char ** argv, dw_opts * s);
void dw_opts_fprint(FILE *f, dw_opts * s);
void dw_opts_free(dw_opts ** sp);
void dw_usage(const int argc, char ** argv, const dw_opts * );

void dw_fprint_info(FILE * f);
void dw_unittests();

int  dw_run(dw_opts *);

#endif
