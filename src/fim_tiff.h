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
    int M, int N, int P);

int fim_tiff_write_zeros(const char * fName, int M, int N, int P);

// Read a 3D tif stack as a float array
float * fim_tiff_read(const char * fName, 
    int * M0, int * N0, int * P0, int verbosity);

// Read a sub region of a 3D stack as float array
// set sub to 1
// reads sM:sM+wM-1, sN:sN+wN-1, sP:sP+wP-1
float * fim_tiff_read_sub(const char * fName, 
    int * M0, int * N0, int * P0, int verbosity,
    int sub,
   int sM, int sN, int sP, // start
   int wM, int wN, int wP); // width

void fim_tiff_ut();

// Get the size of a tiff file (by name)
// Returns 0 upon success.
int fim_tiff_get_size(char * fname, 
    int * M, int * N, int * P);

#endif
