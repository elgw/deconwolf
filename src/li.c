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
 * The method is less accurate than computing the integral using standard
 * numerical integration but should be good enough in most cases.
 *
 * It seems like the method is less accurate for large values of z.
 *
 * For usage, see the main function in this file.
 *
 * 1] Jizhou Li, Feng Xue, and Thierry Blu, "Fast and accurate three-dimensional
 * point spread function computation for fluorescence microscopy,"
 * J. Opt. Soc. Am. A 34, 1029-1034 (2017)
 * https://doi.org/10.1364/JOSAA.34.001029
 *
 * TODO:
 * - Add an option to check the reconstruction error in debug mode.
 * - Consider j0f and j1f (the single precision counterparts to j0 and j1)
 *   for some extra speed.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

typedef struct
{
    // Public, can be set between li_new and the first call to li_calc
    double lambda; // wave length of light
    double NA; // numerical aperture
    double ni; // refractive index
    int M; // Number of sample points in the integral
    int N; // Number of coefficients in the series expansion

    // Internal
    double z; // Axial position
    // Pre-calculated data
    int new; // Indicates if the arrays below are calculated or not
    double * Z; // N scaling values for the Bessel functions
    double * j0Z; // Pre computed Bessel values
    double * j1Z; // Pre computed bessel values
    double * Creal; // Real part of the coefficients
    double * Cimag; // complex part of the coefficients
} li_conf;

li_conf * li_new(double z)
{
    li_conf * L = malloc(sizeof(li_conf));

    // These can be changed by the user
    L->z = z;
    L->lambda = 436;
    L->NA = 1.45;
    L->ni = 1.512;

    // Values below should not be set/changed by the user
    // In the original paper N = 100, and M = 1000;

    L->N = 0; // Dynamic -- this is set at the first call
    L->M = 0; // Dynamic -- this is set at the first call

    L->new = 1;
    L->Z = NULL;
    L->j0Z = NULL;
    L->j1Z = NULL;
    L->Cimag = NULL;
    L->Creal = NULL;
    return L;
}

li_conf * li_free(li_conf ** LP)
{
    li_conf * L = LP[0];
    if(L->Z != NULL)
        free(L->Z);
    if(L->j0Z != NULL)
        free(L->j0Z);
    if(L->j1Z != NULL)
        free(L->j1Z);
    if(L->Cimag != NULL)
        free(L->Cimag);
    if(L->Creal != NULL)
        free(L->Creal);
    free(L);
    return NULL;
}

double li_calc(li_conf * L, const double r)
{

    const double alpha = 2*M_PI/L->lambda*L->NA*r;
    const double beta = 2*M_PI/L->lambda*pow(L->NA,2)/(2*L->ni)*L->z;

    if(L->new == 1)
    {

        L->N = 6+ceil(fabs(beta));
        L->M = 2*L->N;


        gsl_vector * Z = gsl_vector_calloc(L->N);
        L->Z = malloc(L->N*sizeof(double));
        gsl_vector * Freal = gsl_vector_calloc(L->M);
        gsl_vector * Fimag = gsl_vector_calloc(L->M);

        //
        // Set the scaling coefficients (\sigma)
        //


        /* This is the version from the paper
           for(int kk = 0; kk < L->N; kk++)
           {
           double skk = L->NA*(3.0*kk-2.0)*436.0/L->lambda;
           Z->data[kk] = skk;
           }
        */

        double s0 = 0.6;
        Z->data[0] = s0;
        double sd = 0.8;
        double sdd = 0.05;
        for(int kk = 1; kk< L->N; kk++)
        {
            double skk = Z->data[kk-1];
            skk += fmin(M_PI, sd+sdd*(kk-1));
            Z->data[kk] = skk;
        }

        memcpy(L->Z, Z->data, L->N*sizeof(double));

        //
        // Sample the function to expand as a series
        //

        double delta = 1.0/(L->M-1.0);
        for(int kk = 0; kk < L->M; kk++)
        {
            double t = kk*delta;
            Freal->data[kk] = cos(beta*pow(t,2));
            Fimag->data[kk] = sin(beta*pow(t,2));
        }

        gsl_matrix * B = gsl_matrix_calloc(L->M, L->N);

        //
        // Sample the Bessel functions
        //

        for(int cc = 0; cc<L->M; cc++)
        {
            for(int rr = 0; rr<L->N; rr++)
            {
                B->data[cc*L->N + rr] = j0(cc*delta*Z->data[rr]);
            }
        }

        //
        // Find the coefficients
        //

        gsl_vector * Creal = gsl_vector_calloc(L->N);
        gsl_vector * Cimag = gsl_vector_calloc(L->N);
        gsl_matrix * V = gsl_matrix_calloc(L->N, L->N);
        gsl_vector * S = gsl_vector_calloc(L->N);
        gsl_vector * work = gsl_vector_calloc(L->N);
        gsl_linalg_SV_decomp(B, V, S, work);
        gsl_matrix * U = B;
        gsl_linalg_SV_solve(U, V, S, Freal, Creal);
        gsl_linalg_SV_solve(U, V, S, Fimag, Cimag);

        L->Creal = malloc(L->M*sizeof(double));
        L->Cimag = malloc(L->M*sizeof(double));

        for(int kk = 0; kk < L->N; kk++)
        {
            L->Creal[kk] = Creal->data[kk];
            L->Cimag[kk] = Cimag->data[kk];
        }
        L->j0Z = malloc(L->N*sizeof(double));
        L->j1Z = malloc(L->N*sizeof(double));
        for(int kk = 0; kk < L->N; kk++)
        {
            L->j0Z[kk] = j0(L->Z[kk]);
            L->j1Z[kk] = j1(L->Z[kk]);
        }
        L->new = 0;
        gsl_vector_free(Creal);
        gsl_vector_free(Cimag);
        gsl_vector_free(work);
        gsl_matrix_free(V);
        gsl_vector_free(S);
        gsl_matrix_free(B);
        gsl_vector_free(Freal);
        gsl_vector_free(Fimag);
        gsl_vector_free(Z);
    }

    //
    // Calculate the value of the integral
    //
    double vreal = 0;
    double vcomp = 0;

    double u = alpha;
    double j1u = j1(u);
    double j0u = j0(u);
    double * j1Z = L->j1Z;
    double * j0Z = L->j0Z;

    for(int kk = 0; kk<L->N; kk++)
    {
        double v = L->Z[kk];
        double k = (u*j1u*j0Z[kk] - v*j0u*j1Z[kk])/(pow(u,2)-pow(v,2));
        vreal += k*L->Creal[kk];
        vcomp += k*L->Cimag[kk];
    }

    return pow(vreal,2) + pow(vcomp,2);
}

#ifdef LI_TEST
int main(int argc, char ** argv)
{
    // Here is an example of how to use this "class"

    /* Create a new configuration, L, use defaults
       and initialize it for z = 0 nm. */
    double z = 0;
    li_conf * L = li_new(z);
    assert(L != NULL);

    for(double r = 0; r<300; r+=130)
        printf("PSF(r=%f, z=%f) = %f dV\n", r, L->z, li_calc(L, r));

    L = li_free(&L);
    assert(L == NULL);

    return 0;
}
#endif
