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

#ifndef deconwolf_h
#define deconwolf_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <inttypes.h>
#include <string.h>
#include <fftw3.h>
#include <sys/types.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <libgen.h>
#include <time.h>
#include <wchar.h>
#include <locale.h>
#include "fft.h"
#include "tiling.h"
#include "fim.h"
#include "fim_tiff.h"
#include "dw_version.h"

#include "deconwolf_tif_max.h"

/* fftw3 wisdom data is stored and loaded from
 * $home/.config/
 *
 * Internally column major indexing is used, i.e., the distance
 * between elements is 1 for the first dimension and increases with
 * each new dimension.
 */

#define DW_METHOD_W 0
#define DW_METHOD_RL 1
#define DW_METHOD_ID 2

typedef float afloat;

typedef struct{
  int nThreads;
  int nIter;
  char * imFile;
  char * psfFile;
  char * outFile;
  char * logFile;
  char * prefix;
  char * flatfieldFile;
  FILE * log;
  int tiling_maxSize;
  int tiling_padding;
  int overwrite; /* overwrite output if exist */
  int method;
  int verbosity;
  int showTime; /* For dev: show detailed timings */
  fftwf_plan fft_plan;
  fftwf_plan ifft_plan;
  int iterdump; // Dump each iteration to file ...

  float relax;
  float bg; /* Background level, 0 by default */

  /* How aggressive should the Biggs acceleration be.
  *  0 = off,
  *  1 = low/default, safe for most images
  *  2 = intermediate
  *  3 = full, according to the paper
  */
  int biggs;

  int positivity; // Positivity constraint
  float xycropfactor; // discard outer slices that are less than this of the central one
  char * commandline;
  int borderQuality;
  int outFormat; // 16 (=16 bit int) or 32 (=32 bit float)

  /* sigma of Gaussian used for pre filtering of the image and the PSF
   * this was found beneficial a paper by Van Kempen
   *  https://doi.org/10.1046/j.1365-2818.1997.d01-629.x
   * Not used if psigma <= 0 */
  double psigma;

  /* Debug/dev options */
  int onetile; /* For debugging -- only process the first tile if set */
  int experimental1;
  int fulldump; /* write also what is outside of the image */
  /* How far should bigger image sizes be considered? */
  int lookahead;
  int eve; /* Use Exponential vector extrapolation */
} dw_opts;

dw_opts * dw_opts_new(void);
void dw_argparsing(int argc, char ** argv, dw_opts * s);
void dw_opts_fprint(FILE *f, dw_opts * s);
void dw_opts_free(dw_opts ** sp);
void dw_usage(const int argc, char ** argv, const dw_opts * );

void dw_fprint_info(FILE * f, dw_opts * s);
void dw_unittests();

int  dw_run(dw_opts *);

/* Additive Vector Extrapolation (AVE) */
float * deconvolve_ave(afloat * restrict im, const int64_t M, const int64_t N, const int64_t P,
                       afloat * restrict psf, const int64_t pM, const int64_t pN, const int64_t pP,
                       dw_opts * s);

/* Determine Biggs' acceleration parameter alpha  */
float biggs_alpha(const afloat * restrict g,
                  const afloat * restrict gm,
                  const size_t wMNP, int mode);


/* Exponential Vector Extrapolation (EVE) */
float biggs_alpha_eve(const afloat * restrict Xk,
                      const afloat * restrict Xkm1,
                      const afloat * restrict Ukm1,
                      const afloat * restrict Ukm2,
                      const size_t wMNP);

/* Exponential vector extrapolation (eve) alpha */
float biggs_alpha_eve(const afloat * restrict Xk,
                      const afloat * restrict Xkm1,
                      const afloat * restrict Ukm1,
                      const afloat * restrict Ukm2,
                  const size_t wMNP);


/* Autocrop the PSF by:
 * 1/ Cropping if the size is larger than needed by the image.
 * 2/ Optionally, trim the sides in x and y where the PSF vanishes.
 */
float * psf_autocrop(float * psf, int64_t * pM, int64_t * pN, int64_t * pP,  // psf and size
    int64_t M, int64_t N, int64_t P, // image size
    dw_opts * s);

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"


/* Write A to disk as fulldump_<name>.tif if s->fulldump = 1 */
void fulldump(dw_opts * s, float * A, size_t M, size_t N, size_t P, char * name);


#endif
