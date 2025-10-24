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


/* This programs calculate Points Spread Functions for deconvolution of
 * wide field microscopy images (stacks).
 * Original implementation:
 * https://github.com/Biomedical-Imaging-Group/PSFGenerator/tree/master/src/psf/bornwolf
 */

#include <stdio.h>

/* Mode for calculating the 1D integral */
typedef enum {
    MODE_BW_GSL,
    MODE_BW_LI } dw_bw_integral_type;

typedef struct {
    int verbose;
    int overwrite;

    char * cmd; // Command line
    char * outFile;
    char * logFile;
    FILE * log;

    /* Physical parameters */
    float lambda; // Emission maxima
    float NA; // Numerical aperture of lens
    float ni;
    float resLateral; // pixel size in x and y
    float resAxial; // pixel size in z

    /* Integration options */
    dw_bw_integral_type mode_bw; /* How to calculate the 1D integral */
    int oversampling_R; // Oversampling in radial direction
    double epsabs;
    double epsrel;
    size_t limit;
    int key;

    // For parallel processing
    int thread;
    int nThreads;

    int testing;

    /* Outputs */
    float * V; // output image
    // shape of output image
    int M;
    int N;
    int P;
} bw_conf;

/* New configuration with default settings */
bw_conf * bw_conf_new(void);

/* Check the configuration, and possibly make some minor adjustments
 * returns EXIT_SUCCESS if the config looks valid. */
int bw_conf_validate(bw_conf * s);

/* Print out the configuration */
void bw_conf_printf(FILE * out, bw_conf * conf);

/* Free conf and all fields that it points to and finally sets it to NULL */
void bw_conf_free(bw_conf ** conf);

/* Generate the Born-Wolf PSF
 */
void BW(bw_conf * conf);

/* CLI */
int dw_bwpsf(int argc, char ** argv);
