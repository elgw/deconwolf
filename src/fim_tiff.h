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

#ifndef fim_tiff_h
#define fim_tiff_h
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <tiffio.h>
#include <unistd.h>
#include "fim.h"

#define INLINED inline __attribute__((always_inline))

int fim_tiff_write(const char * fName, const float * V, 
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
    int64_t * M0, int64_t * N0, int64_t * P0, int verbosity);

// Read a sub region of a 3D stack as float array
// set sub to 1
// reads sM:sM+wM-1, sN:sN+wN-1, sP:sP+wP-1
float * fim_tiff_read_sub(const char * fName, 
    int64_t * M0, int64_t * N0, int64_t * P0, int verbosity,
    int sub,
   int64_t sM, int64_t sN, int64_t sP, // start
   int64_t wM, int64_t wN, int64_t wP); // width

void fim_tiff_ut();

// Get the size of a tiff file (by name)
// Returns 0 upon success.
int fim_tiff_get_size(char * fname, 
    int64_t * M, int64_t * N, int64_t * P);

#endif
