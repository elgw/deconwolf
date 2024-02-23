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

/* Calculates the value of a Born-Wolf PSF at a point specified by axial (z) and
 * radial (r) position. This piece of code uses the method described in [1] to
 * speed up the calculations of the integral.
 *
 * The accuracy depends on the number of sample points and the choice of basis
 * functions. More sample points and basis functions needed for large r and z.
 *
 * For usage, see the 'main' function in this file.
 *
 * [1] Jizhou Li, Feng Xue & Thierry Blu, "Fast and accurate three-dimensional
 *     point spread function computation for fluorescence microscopy,"
 *     J. Opt. Soc. Am. A 34, 1029-1034 (2017)
 *     https://doi.org/10.1364/JOSAA.34.001029
 *
 * TODO:
 * - Add an option to check the reconstruction error in debug mode.
 * - Consider j0f and j1f (the single precision counterparts to j0 and j1)
 *   for some extra speed.
 */

#include <stdlib.h>
#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>
#include <string.h>
#include <complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#ifndef dcomplex
#ifdef WINDOWS
typedef _Dcomplex dcomplex;
#else
typedef double complex dcomplex;
#endif
#endif

typedef struct
{
    /* Public, can be set between li_new and the first call to li_calc */
    double lambda; /* wave length of light */
    double NA; /* numerical aperture */
    double ni; /* refractive index */
    int M; /* Number of sample points in the integral */
    int N; /* Number of coefficients in the series expansion */

    /* Internal, don't fiddle */
    double z; /* Axial position */
    /* Pre-calculated data */
    int new; /* Indicates if the arrays below are calculated or not */
    double * Z; /* N scaling values for the Bessel functions */
    double * j0Z; /* Pre computed Bessel values */
    double * j1Z; /* Pre computed bessel values */
    double * Creal; /* Real part of the coefficients */
    double * Cimag; /* complex part of the coefficients */
} li_conf;

li_conf * li_new(double z);

li_conf * li_free(li_conf ** LP);

dcomplex li_calc(li_conf * L, const double r);

void li_show(li_conf * L);
