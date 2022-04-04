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
float fwhm_lateral(fim_image_t * I,
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

#ifdef STANDALONE
static void createGaussian(double * , float * , size_t , float );
#endif


#endif
