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

int fim_tiff_write(char * fName, float * V, 
    int M, int N, int P);

// Read a 3D tif stack as a float array
float * fim_tiff_read(char * fName, 
    int * M0, int * N0, int * P0, int verbosity);

void fim_tiff_ut();

#endif
