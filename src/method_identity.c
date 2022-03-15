#include "method_identity.h"


 float * deconvolve_identity(afloat * restrict im, const int64_t M, const int64_t N, const int64_t P,
                        afloat * restrict psf, const int64_t pM, const int64_t pN, const int64_t pP,
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

     fftwf_free(psf);
     afloat * out = fftwf_malloc(M*N*P*sizeof(afloat));
     memcpy(out, im, M*N*P*sizeof(afloat));
     return out;
 }
