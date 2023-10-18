#pragma once

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <gsl/gsl_integration.h>

typedef struct {
    double NA;
    double ni;

    double lambda;
    size_t ncalls; /* Counting the number of function calls */

    /* Integration quality */
    size_t limit;      /* Max number of sub-intervals */
    double epsabs; /* Absolute error */
    double epsrel; /* Relative error */
    int key; /* Integration method for xyz*/
    int keybw; /* Integration method for the bw integral */
    /* 1-6 corresponds to
     * 15, 21, 31, 41, 51 and 61 point Gauss-Kronrod rules
     * GSL_INTEG_GAUSS15 clearly fastest in comparison for this
     * problem.
     */

    /* One workspace per dimension since the integrals are nested */
    gsl_integration_workspace *wspx, *wspy, *wspz, *wspb;

    /* Integration region specification */
    double x;
    double y;
    double z;
    double r;
    double x0;
    double x1;
    double y0;
    double y1;
    double z0;
    double z1;

} bw_gsl_conf_t;

bw_gsl_conf_t * bw_gsl_new(size_t ninterval);
void bw_gsl_free(bw_gsl_conf_t * );
void bw_gsl_fprint(FILE *, bw_gsl_conf_t *);

/* Integrate bw at a r, z */
double bw_gsl_integrate(bw_gsl_conf_t *,
                        double r, double z);

double bw_gsl_integrate_y(bw_gsl_conf_t * conf,
                          double x, double y0, double y1, double z);

double bw_gsl_integrate_z(bw_gsl_conf_t *,
                          double r, double z0, double z1);

double bw_gsl_integrate_xyz(bw_gsl_conf_t *,
                            double x0, double x1,
                            double y0, double y1,
                            double z0, double z1);

double bw_gsl_integrate_xy(bw_gsl_conf_t *,
                           double x0, double x1,
                           double y0, double y1,
                           double z);

#ifndef j0f
#define j0f j0
#endif
