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

#include "li.h"

void li_show(li_conf * L)
{
    printf("lambda = %f\n", L->lambda);
    printf("NA = %f\n", L->NA);
    printf("ni = %f\n", L->ni);
    printf("z = %f\n", L->z);
    return;
}

li_conf * li_new(double z)
{
    li_conf * L = malloc(sizeof(li_conf));

    /* These can be changed by the user */
    L->z = z;
    L->lambda = 461;
    L->NA = 1.45;
    L->ni = 1.512;

    /* Values below should not be set/changed by the user
     * In the original paper N = 100, and M = 1000; */

    L->N = 0; /* Dynamic -- this is set at the first call */
    L->M = 0; /* Dynamic -- this is set at the first call */

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

double complex li_calc(li_conf * L, const double r)
{

    /* Updated 2021-11-18 */
    const double alpha = 2*M_PI/L->lambda*L->NA*r;
    const double beta = 2*M_PI/L->lambda*pow(L->NA,2)/(2*L->ni)*L->z;

    // const double alpha = (2.0*M_PI/L->lambda)*L->NA/L->ni*r;
    // const double beta =  (2.0*M_PI/L->lambda)*pow(L->NA/L->ni,2)/2.0*L->z;

    if(L->new == 1)
    {
        L->N = 6+ceil(fabs(beta));
        L->M = 2*L->N;

        gsl_vector * Z = gsl_vector_calloc(L->N);
        L->Z = malloc(L->N*sizeof(double));
        gsl_vector * Freal = gsl_vector_calloc(L->M);
        gsl_vector * Fimag = gsl_vector_calloc(L->M);

        /*
         * Set the scaling coefficients (\sigma)
         */


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

        /*
         * Sample the function to expand as a series
         */

        double delta = 1.0/(L->M-1.0);
        for(int kk = 0; kk < L->M; kk++)
        {
            double t = kk*delta;
            Freal->data[kk] = cos(beta*pow(t,2));
            Fimag->data[kk] = sin(beta*pow(t,2));
        }

        gsl_matrix * B = gsl_matrix_calloc(L->M, L->N);

        /*
         * Sample the Bessel functions
         */

        for(int cc = 0; cc<L->M; cc++)
        {
            for(int rr = 0; rr<L->N; rr++)
            {
                B->data[cc*L->N + rr] = j0(cc*delta*Z->data[rr]);
            }
        }

        /*
         * Find the coefficients
         */

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

    /*
     * Calculate the value of the integral
     */

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
        double k = (u*j1u*j0Z[kk] - v*j0u*j1Z[kk])/(pow(u, 2)-pow(v, 2));
        vreal += k*L->Creal[kk];
        vcomp += k*L->Cimag[kk];
    }

    return vreal + I*vcomp;
}

#ifdef LI_TEST
int main(int argc, char ** argv)
{
    /* Here is an example of how to use this "class" */

    /* Create a new configuration, L, use defaults
       and initialize it for z = 0 nm. */
    double z = 0;
    li_conf * L = li_new(z);

    li_show(L);
    assert(L != NULL);

    for(double r = 0; r<500; r+=130)
    {
        double complex v = li_calc(L, r);
        double vreal = pow(creal(v),2) + pow(cimag(v),2);
        printf("LI(r=%f, z=%f) = (%f %fi) %f dV\n", r, L->z, creal(v), cimag(v), vreal);
    }

    L = li_free(&L);
    assert(L == NULL);

    return 0;
}
#endif /* #ifdef LI_TEST */
