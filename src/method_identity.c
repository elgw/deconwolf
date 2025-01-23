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

#include "method_identity.h"

 float * deconvolve_identity(float * restrict im, const int64_t M, const int64_t N, const int64_t P,
                        float * restrict psf, const int64_t pM, const int64_t pN, const int64_t pP,
                        dw_opts * s)
 {

     if(s->verbosity > 1)
     {
         printf("method_identity: Doing nothing\n");
     }

     if(s->verbosity > 0)
     {
         printf("image: [%" PRId64 "x%" PRId64 "x%" PRId64
              "], psf: [%" PRId64 "x%" PRId64 "x%" PRId64
              "]\n",
              M, N, P, pM, pN, pP);
     }

     fprintf(s->log, "image: [%" PRId64 "x%" PRId64 "x%" PRId64 "]\n"
             "psf: [%" PRId64 "x%" PRId64 "x%" PRId64 "]\n",
             M, N, P, pM, pN, pP);
     fflush(s->log);

     fim_free(psf);
     float * out = fim_malloc(M*N*P*sizeof(float));
     memcpy(out, im, M*N*P*sizeof(float));
     return out;
 }
