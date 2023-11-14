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

typedef struct {
    float fwhm_lateral;
    float fwhm_axial;
} fwhm_res_t;

/** @brief FWHM for N points in I */
fwhm_res_t * fwhm(fim_t * I,
                  const int * X, const int * Y, const int * Z,
                  size_t N);

/** @brief Calculate the lateral fwhm for one point */
float fwhm_lateral(fim_t * I,
                   int x, int y, int z,
                   int verbose);

/** @brief Calculate the lateral fwhm for one point */
float fwhm_axial(fim_t * I,
                   int x, int y, int z,
                   int verbose);


/** @brief run some self tests */
int fwhm_ut(int argc, char ** argv);
