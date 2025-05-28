#pragma once

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

#include <assert.h>
#include <inttypes.h>
#include <fftw3.h>
#include <getopt.h>
#include <locale.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <wchar.h>


#ifdef _OPENMP // turned on with -fopenmp
#include <omp.h>
#endif
#ifdef MKL
#include <mkl.h>
#endif

#include "dw_util.h"
#include "dw_maxproj.h"
#include "dw_tiff_merge.h"
#include "dw_imshift.h"
#include "dw_version.h"
#ifndef WINDOWS
#include "dw_nuclei.h"
#endif
#include "dw_background.h"

/* Uncomment to include, requires linking with libpng and libz
 * can be build separately by the makefile in the src folder
 */
// #include "dw_otsu.h"
#ifndef WINDOWS
#include "dw_dots.h"
#include "dw_psf.h"
#include "dw_psf_sted.h"
#include "dw_align_dots.h"
#endif

#include "fim.h"
#include "fim_tiff.h"
#include "fft.h"
#include "tiling.h"
#include "sparse_preprocess_cli.h"

/* fftw3 wisdom data is stored and loaded from
 * $home/.config/
 *
 * Internally column major indexing is used, i.e., the distance
 * between elements is 1 for the first dimension and increases with
 * each new dimension.
 */

typedef enum {
    DW_METHOD_RL, /* Richardson Lucy */
    DW_METHOD_ID, /* Identity/nothing. For checking image loading/saving. */
    DW_METHOD_SHB, /* Wang and Miller, Scaled Heavy Ball */
    #ifdef OPENCL
    DW_METHOD_SHBCL, /* GPU used only for FFT */
    DW_METHOD_SHBCL2, /* GPU used as much as possible */
    #endif
} dw_method;

/* What should be the initial guess
   for the iterations? */
typedef enum {
    /* The average of the input image. Was default up till version
     * 0.3.7 */
    DW_START_FLAT,
    /* A low pass filtered version of the input image */
    DW_START_LP,
    /* The input image itself */
    DW_START_IDENTITY,
} dw_start_condition;

typedef enum {
    DW_METRIC_MSE, /* Mean Square Error */
    DW_METRIC_IDIV /* I-Divergence */
} dw_metric;

typedef enum {
    DW_ITER_ABS,
    DW_ITER_REL,
    DW_ITER_FIXED
} dw_iter_type;


struct _dw_opts; /* Forward declaration */
typedef struct _dw_opts dw_opts;

typedef float * (*dw_function) (float * restrict im, const int64_t M, const int64_t N, const int64_t P,
                              float * restrict psf, const int64_t pM, const int64_t pN, const int64_t pP,
                              dw_opts * s);

struct _dw_opts{
    int nThreads_FFT;
    int nThreads_OMP;

    char * imFile;
    char * psfFile;
    char * refFile; /* Name of reference image */
    char * tsvFile; /* Where to write tsv benchmark data */
    float * ref; /* Reference image */
    char * outFile;
    char * logFile;
    char * outFolder; /* Where iterdump files should go etc */
    char * prefix;
    char * flatfieldFile;
    FILE * log;
    FILE * tsv;
    int tiling_maxSize;
    int tiling_padding;
    int overwrite; /* overwrite output if exist */

    int nIter_auto; /* Automatic stopping? */
    int nIter; /* Fixed number of iterations, used when nIter_auto = 0 */
    int maxiter; /* Max number of iter for rel and abs mode */
    float err_rel;
    float err_abs;
    dw_iter_type iter_type;

    int verbosity;
    int color; /* Show colored things in terminal */
    int showTime; /* For dev: show detailed timings */
    int iterdump; /* Dump each iteration to file */

    float bg; /* Background level, set automatically by default */
    int bg_auto; /* Indicator if background was set on CLI */
    float offset;
    float psf_pass; /* For low-pass filtering by the PSF */

    /* Selection of method */
    dw_method method; /* what algorithm to use */
    dw_function fun; /* Function pointer */
    dw_metric metric;

    /* Select what the initial guess should be */
    dw_start_condition start_condition;
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
    float scaling; // fixed scaling for 16 bit output, automatic scaling is used if this value <= 0

    int auto_zcrop; // Set to > 0 for an attempt to automatically crop out the middle planes of the image
    int zcrop; // How many planes to remove from top and bottom

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

    float alphamax;


    int cl_device; /* OpenCL device number, default 0 */


    fftwf_plan fft_plan;
    fftwf_plan ifft_plan;
    int fftw3_planning;
    int fft_inplace;
    struct timespec tstart;
};


dw_opts * dw_opts_new(void);
void dw_argparsing(int argc, char ** argv, dw_opts * s);
void dw_opts_fprint(FILE *f, dw_opts * s);
void dw_opts_free(dw_opts ** sp);
void dw_usage(const int argc, char ** argv, const dw_opts * );

void dw_fprint_info(FILE * f, dw_opts * s);
void dw_unittests();

int  dw_run(dw_opts *);

/* Additive Vector Extrapolation (AVE) */
float * deconvolve_ave(float * restrict im, const int64_t M, const int64_t N, const int64_t P,
                       float * restrict psf, const int64_t pM, const int64_t pN, const int64_t pP,
                       dw_opts * s);


/* Determine Biggs' acceleration parameter alpha  */
float biggs_alpha(const float * restrict g,
                  const float * restrict gm,
                  const size_t wMNP, int mode);


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
#define ANSI_UNDERSCORE    "\x1b[4m"


/* Write A to disk as fulldump_<name>.tif if s->fulldump = 1 */
void fulldump(dw_opts * s, float * A, size_t M, size_t N, size_t P, char * name);

/* Passes on to got_fMSE at the moment */
float getError(const float * restrict y, const float * restrict g,
               const int64_t M, const int64_t N, const int64_t P,
               const int64_t wM, const int64_t wN, const int64_t wP, dw_metric);


/* Mean squared error between the input, y, and the forward propagated
 * current guess */
float get_fMSE(const float * restrict y, const float * restrict g,
               const int64_t M, const int64_t N, const int64_t P,
               const int64_t wM, const int64_t wN, const int64_t wP);

/* Idiv between the input, y, and the forward propagated current
 * guess */
float get_fIdiv(const float * restrict y, const float * restrict g,
                const int64_t M, const int64_t N, const int64_t P,
                const int64_t wM, const int64_t wN, const int64_t wP);

/* Show the current iteration, can be called by all method_ */
void dw_show_iter(dw_opts * s, int it, int nIter, float error);

/* Create initial guess: the fft of an image that is 1 in MNP and 0 outside
 * M, N, P is the dimension of the microscopic image
 *
 * Possibly more stable to use the mean of the input image rather than 1
 */
fftwf_complex * initial_guess(const int64_t M, const int64_t N, const int64_t P,
                              const int64_t wM, const int64_t wN, const int64_t wP);

/* Generate a name for the an iterdump file
   at iteration it */
char * gen_iterdump_name(
    __attribute__((unused)) const dw_opts * s,
    int it);


/* UTILS */
double clockdiff(struct timespec* end, struct timespec * start);
/* Show a green dot and flush stdout */
void putdot(const dw_opts *s);
int64_t int64_t_max(int64_t a, int64_t b);

/* Write diagostics to s->tsv if open
 * To do: add timings as well (excluding) this function
*/
void benchmark_write(dw_opts * s, int iter, double fMSE, const float * x,
                     const int64_t M, const int64_t N, const int64_t P,
                     const int64_t wM, const int64_t wN, const int64_t wP);



/* "WARNING" in some formatting */
void warning(FILE * fid);


typedef struct{
    float error; /* Current error */
    float lasterror; /* Last error */
    int iter; /* Current iteration */
    int niter; /* Max number of iterations */
    float relerror; /* Relative error to stop at */
    float abserror; /* Absolute error to stop at */
    dw_iter_type itertype; /* Stop condition class */
} dw_iterator_t;

dw_iterator_t * dw_iterator_new(const dw_opts *);
int dw_iterator_next(dw_iterator_t * );
void dw_iterator_set_error(dw_iterator_t *, float);
void dw_iterator_show(dw_iterator_t *, const dw_opts *);
void dw_iterator_free(dw_iterator_t * );

#ifdef OPENCL
#include "method_shb_cl.h"
#include "method_shb_cl2.h"
#endif

#include "method_identity.h"
#include "method_rl.h"
#include "method_shb.h"
