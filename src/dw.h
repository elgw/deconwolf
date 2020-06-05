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
#include <string.h>
#include <fftw3.h>
#include <sys/types.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <libgen.h>
#include <time.h>
#include "fft.h"
#include "tiling.h"
#include "fim.h"
#include "fim_tiff.h"
#include "dw_version.h"


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
  float relax;
  float xycropfactor; // discard outer slices that are less than this of the central one
  char * commandline;

  int onetile; // For debugging -- only process the first tile if set
} dw_opts;

dw_opts * dw_opts_new(void);
void dw_argparsing(int argc, char ** argv, dw_opts * s);
void dw_opts_fprint(FILE *f, dw_opts * s);
void dw_opts_free(dw_opts ** sp);
void dw_usage(const int argc, char ** argv, const dw_opts * );

void dw_fprint_info(FILE * f, dw_opts * s);
void dw_unittests();

int  dw_run(dw_opts *);

/* Autocrop the PSF by:
 * 1/ Cropping if the size is larger than needed by the image.
 * 2/ Optionally, trim the sides in x and y where the PSF vanishes.
 */
float * psf_autocrop(float * psf, int64_t * pM, int64_t * pN, int64_t * pP,  // psf and size
    int64_t M, int64_t N, int64_t P, // image size
    dw_opts * s);


#endif
