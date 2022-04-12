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

#ifndef __fwhm_h__
#define __fwhm_h__

/* Only thread safe when logFile == NULL */

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_spline.h>
#include "fim.h"

struct dbg_data {
    double xm; // location of maxima
    double ym; // value at maxima
    double xleft;
    double xright;
};

/* Calculate the fwhm for the y at the locations x
 * Stores the result in fwhm
 *
 * Algorithm:
 * 1. Find the interpolated location of the maxima.
 * 2. Use all pixels in the profile to find the background level as the
 *    smallest value
 * 3. Use approx 50% of the central pixels to find the zero crossing of
 *    the 50% value and the profile.
 * 4. Return the distance between the left and right zero crossing if both
 *    look ok, or use 2 mid-zero crossing of only one found.
*/
int fwhm1d(const double * x, const double * y, size_t N, double * fwhm);

/* Calculate the lateral fwhm for one point */
float fwhm_lateral(fim_t * I,
                   int x, int y, int z,
                   int verbose);

// my_f and my_f_params defines the interpolation function that the
// root-finding functions are run on

struct my_f_params {
    gsl_spline * spline;
    gsl_interp_accel * acc;
    double offset; // shifts the function
};

double my_f(double , void * );
int findmin(double , double , gsl_spline *, gsl_interp_accel *, double * , double * );

#ifdef STANDALONE_FWHM
static void createGaussian(double * , double * , size_t , float );
#endif


#endif
