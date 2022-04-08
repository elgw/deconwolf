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

#ifndef __method_identity_h__
#define __method_identity_h__

#include "dw.h"

/* Exponential Vector Extrapolation (EVE) */
float * deconvolve_identity(float * restrict im, const int64_t M, const int64_t N, const int64_t P,
                       float * restrict psf, const int64_t pM, const int64_t pN, const int64_t pP,
                       dw_opts * s);

#endif
