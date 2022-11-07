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

/* Read and write tiff files to/from single precision floats
 */

/* fim_tiff is not thread safe
 * You need to initialize with a call to
 * fim_tiff_init()
 * and should probably also redirect the output by
 * fim_tiff_set_log(FILE *)
 *
 * TODO:
 * Flags to dw_write_tif for scaling on/off
*/

#ifndef fim_tiff_h
#define fim_tiff_h

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <tiffio.h>
#include <unistd.h>
#include <fftw3.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include "fim.h"
#include "ftab.h"
#include "dw_version.h"

#define INLINED inline __attribute__((always_inline))
#define XTAG_IJIJUNKNOWN 50838
#define XTAG_IJIJINFO 50839


/* Tiff tags -- for simple transfer from one image to another */
typedef struct{
    float xresolution;
    float yresolution;
    float zresolution;
    char * imagedescription;
    char * software;
    uint16_t resolutionunit;
    char * IJIJinfo; // Tag 50839 contains a string, used by Imagej.
    uint32_t nIJIJinfo;
    // Image size
    int M;
    int N;
    int P;
} ttags;

// new with everything set to defaults
ttags * ttags_new();
void ttags_get(TIFF *, ttags *);
void ttags_show(FILE *, ttags *);
void ttags_set(TIFF *, ttags *);
void ttags_set_software(ttags * , char *);
void ttags_set_imagesize(ttags *, int M, int N, int P);
void ttags_set_pixelsize(ttags *, double, double, double);
void ttags_free(ttags **);

/* Initialization, sets the output file to stdout */
void fim_tiff_init(void);

/* Redirect all output here */
void fim_tiff_set_log(FILE * fp);

/* Scale and write data */
int fim_tiff_write(const char * fName, const float * V,
                   ttags * T,
    int64_t M, int64_t N, int64_t P);

/* Don't scale data */
int fim_tiff_write_noscale(const char * fName, const float * V,
                           ttags * T,
                           int64_t N, int64_t M, int64_t P);

int fim_tiff_write_float(const char * fName, const float * V,
                         ttags * T,
    int64_t M, int64_t N, int64_t P);

int fim_tiff_write_zeros(const char * fName, int64_t M, int64_t N, int64_t P);

// Write to the tif file fName as uint16, using the raw data
// in rName
int fim_tiff_from_raw(const char * fName, int64_t M, int64_t N, int64_t P,
    const char * rName);

// Convert tiff image to raw float image
int fim_tiff_to_raw(const char *fName, const char * oName);

// Read a 3D tif stack as a float array
float * fim_tiff_read(const char * fName,
                      ttags * T,
    int64_t * M0, int64_t * N0, int64_t * P0, int verbosity);

// Read a sub region of a 3D stack as float array
// set sub to 1
// reads sM:sM+wM-1, sN:sN+wN-1, sP:sP+wP-1
float * fim_tiff_read_sub(const char * fName,
                          ttags *,
    int64_t * M0, int64_t * N0, int64_t * P0, int verbosity,
    int sub,
   int64_t sM, int64_t sN, int64_t sP, // start
   int64_t wM, int64_t wN, int64_t wP); // width

void fim_tiff_ut();

// Get the size of a tiff file (by name)
// Returns 0 upon success.
int fim_tiff_get_size(char * fname,
    int64_t * M, int64_t * N, int64_t * P);

/* Max projection from input to output file */
int fim_tiff_maxproj(char * in, char * out);

/* Extract a single slice from input to output file */
int fim_tiff_extract_slice(char *in, char *out, int slice);

#endif
