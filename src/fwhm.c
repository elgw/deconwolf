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

#include "fwhm.h"

/*
 * Forward declarations
 */

static void createGaussian(
    double * restrict x,
    double * restrict y,
    size_t N,
    const double x0);
static double my_f(double , void * );
static int findmin(double , double , gsl_spline *, gsl_interp_accel *, double * , double * );
static int fwhm1d(const double * x,
                  const double * y,
                  size_t N,
                  double * fwhm);


/*
 * End of forward declarations
 */

struct dbg_data {
    double xm; // location of maxima
    double ym; // value at maxima
    double xleft;
    double xright;
};

// my_f and my_f_params defines the interpolation function that the
// root-finding functions are run on

struct my_f_params {
    gsl_spline * spline;
    gsl_interp_accel * acc;
    double offset; // shifts the function
};


static void * ck_calloc(size_t nmemb, size_t size)
{
    void * p = calloc(nmemb, size);
    if(p == NULL)
    {
        fprintf(stderr, "ERROR: fwhm failed to allocate memory\n");
        exit(EXIT_FAILURE);
    }
    return p;
}

static double vec_min(const double * v, size_t N)
{
    double min = v[0];
    for(size_t kk = 1; kk<N; kk++)
    {
        v[kk] < min ? min = v[kk] : 0;
    }
    return min;
}

/* Search for the intersection of the signal and
 * bg + 0.5*(ym-bg)
 */
static int get_intersection(const double x_mid, const double x_edge,
                            const double bg, const double ym,
                            gsl_spline *spline,
                            gsl_interp_accel *acc,
                            gsl_root_fsolver * rsolve,
                            double * pos)
{
    /* TODO: move x_edge close to x_mid until only one zero crossing
     * remain. Or just call with smaller intervals ...
     */

    double x_lo0 = x_mid;
    double x_hi0 = x_edge;
    if(x_lo0 > x_hi0)
    {
        x_lo0 = x_edge;
        x_hi0 = x_mid;
    }

    struct my_f_params f_params;
    f_params.spline = spline;
    f_params.acc = acc;
    f_params.offset = (bg+(-ym-bg)/2.0);

    gsl_function F;
    F.function = &my_f;
    F.params = (void *) &f_params;

    /*
     * Function: int gsl_root_fsolver_set (gsl_root_fsolver * s,
     * gsl_function * f, double x_lower, double x_upper) This function
     * initializes, or reinitializes, an existing solver s to use the
     * function f and the initial search interval [x_lower, x_upper]
     */

    int status =  gsl_root_fsolver_set(rsolve, &F, x_lo0, x_hi0);

    size_t iter = 0;
    size_t max_iter = 10;

#if verbose > 0
    double r_expected = x_lo0 + 0.5*(x_hi0-x_lo0);
    printf ("using %s method\n", gsl_root_fsolver_name (rsolve));
    printf ("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "root", "err", "err(est)");
#endif

    double x_lo = x_lo0;
    double x_hi = x_hi0;
    do {
        iter++;
        status = gsl_root_fsolver_iterate (rsolve);
        // update bounds
        x_lo = gsl_root_fsolver_x_lower (rsolve);
        x_hi = gsl_root_fsolver_x_upper (rsolve);
        // see if converged
        status = gsl_root_test_interval (x_lo, x_hi, 0, 0.001);
#if verbose > 0
        double r = gsl_root_fsolver_root (rsolve);
        if (status == GSL_SUCCESS)
            printf ("Converged:\n");
        printf ("%5lu [%.7f, %.7f] %.7f %+.7f %.7f\n", iter, x_lo, x_hi,
                r, r - r_expected,
                x_hi - x_lo);
#endif
    }
    while (status == GSL_CONTINUE && iter < max_iter);

#if verbose > 0
    printf("Found intersection at %f\n", (x_lo+x_hi)/2);
#endif

    pos[0] = (x_lo+x_hi)/2;

    /* Validate the intersection, it is unvalid if values > ym can be
     * found between x_lo and x_hi
     */

    double elow = pos[0];
    double ehigh = x_mid;
    if(elow > ehigh)
    {
        elow = x_mid;
        ehigh = pos[0];
    }

    const double neval_points = 7;
    int valid = 1;
    double delta = (ehigh-elow)/(neval_points + 2);
    double ymax = -ym;
    //printf("[%f, %f], %f\n", elow, ehigh, bg + 0.5*(ymax-bg));
    for(double x = elow+delta; x < ehigh-delta; x+=delta)
    {
        double y =  gsl_spline_eval(spline, x, acc);
        if(y > ymax)
        {
            //printf("x=%f, y=%f > ym=%f\n", x, y, ymax);
            valid = 0;
        }
        if(y < bg + 0.5*(ymax-bg))
        {
            //printf("x=%f, y=%f < mid=%f\n", x, y, 0.5*(ymax-bg));
            valid = 0;
        }
    }

    if(valid)
    {
        return EXIT_SUCCESS;
    } else {
        pos[0] = -1;
        return EXIT_FAILURE;
    }
}

static double my_f(double x, void * p)
{
    struct my_f_params * params = (struct my_f_params*) p;
#if verbose > 0
    printf("f(%f) + %f=  ", x, params->offset); fflush(stdout);
#endif

    double y =  -gsl_spline_eval(params->spline, x, params->acc) + params->offset;
#if verbose > 0
    printf("%f\n", y); fflush(stdout);
#endif

    return(y);
}

static int findmin(double a, double b,
                   gsl_spline *spline, gsl_interp_accel * acc,
                   double * xm, double * ym)
{

/* domain to test [a, b]*/
    int iter = 0;
    int max_iter = 25;
    double m = (a+b)/2.0; /* expected location of minima */

    int status;
    gsl_min_fminimizer *s;

    struct my_f_params f_params;
    f_params.spline = spline;
    f_params.acc = acc;
    f_params.offset = 0;

    gsl_function F;
    F.function = &my_f;
    F.params = (void *) &f_params;

    const gsl_min_fminimizer_type * T = gsl_min_fminimizer_goldensection;
    //const gsl_min_fminimizer_type * T = gsl_min_fminimizer_brent;
    s = gsl_min_fminimizer_alloc (T);
    if(0){
        fprintf(stderr, "Starting with f(%f) = %f, f(%f) = %f, f(%f) = %f\n",
                a , gsl_spline_eval(spline, a, acc),
                m , gsl_spline_eval(spline, m, acc),
                b , gsl_spline_eval(spline, b, acc));
        fflush(stderr);
    }
    if(gsl_min_fminimizer_set (s, &F, m, a, b) == GSL_EINVAL)
    {
        // If no minima is enclosed ... i.e. the value at (a+b)/2 isn't lower than at a and b
        gsl_min_fminimizer_free(s);
        return EXIT_FAILURE;
    }

#if verbose > 0
    fflush(stdout);
    printf("Interval: [%f, %f]\n", a, b);
    double m_expected = -10;
    printf ("method: '%s'\n", gsl_min_fminimizer_name (s));
    printf ("%5s [%9s, %9s] %9s %10s %9s\n", "iter", "lower", "upper", "min", "err", "err(est)");
    printf ("%5d [%.7f, %.7f] %.7f %+.7f %.7f\n", iter, a, b, m, m - m_expected, b - a);
#endif

    const double epsabs = 1e-4;
    const double epsrel = 0.0;

    do {
        iter++;
        status = gsl_min_fminimizer_iterate(s);
        // GSL_FAILURE, GSL_EBADFUNC
        a = gsl_min_fminimizer_x_lower(s);
        b = gsl_min_fminimizer_x_upper(s);

        status = gsl_min_test_interval(a, b, epsabs, epsrel);

        //double argmax = gsl_min_fminimizer_x_minimum(s);
        //printf ("%5d [%.7f, %.7f] " "%.7f %.7f\n", iter, a, b, argmax, b - a);

    } while (status == GSL_CONTINUE && iter < max_iter);

    xm[0] = gsl_min_fminimizer_x_minimum(s);
    ym[0] = my_f(xm[0], &f_params);
    gsl_min_fminimizer_free(s);

    return EXIT_SUCCESS;
}

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
int fwhm1d(const double * x, const double * y, size_t N, double * w)
{
    int retvalue = EXIT_SUCCESS;
    w[0] = -1; /* To be updated if successful */

    gsl_set_error_handler_off();

    /* Basic sanity check */
    if(N%2 == 0)
    {
        fprintf(stderr, "Can only use arrays of odd size\n");
        retvalue = EXIT_FAILURE;
        goto exit2;
    }
    int N2 = (N-1)/2; /* an even number */
    if(y[N2] < y[N2+1] || y[N2] < y[N2-1])
    {
        //printf("%f %f %f\n", y[N2-1], y[N2], y[N2+1]);
        //fprintf(stderr, "The central pixel is not a local maxima\n");
        retvalue = EXIT_FAILURE;
        goto exit2;
    }

    /* Set up interpolation */
    gsl_interp_accel *acc = gsl_interp_accel_alloc();  // Neccessary?
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, N);
    gsl_spline_init(spline, x, y, N);


#if verbose > 0
    /* Double check the interpolation. See that the correct value is
     * given at x=0 */
    printf("y[%f] = %f\n", x[N2], gsl_spline_eval(spline, x[N2], acc));
    fflush(stdout);
#endif

    /* Determine background */
    double bg = vec_min(y, N);
#if verbose > 0
    printf("Background level: %f\n", bg);
#endif

    if(y[N2] <= bg)
    {
        fprintf(stderr, "Signal looks constant\n");
        retvalue = EXIT_FAILURE;
        goto exit1;
    }

    /* Find the sub pixel position of the max, i.e., centre in [xmin, xmax] */
    double xmin = -1.0;
    double xmax = 1.0;

    double xm = 0;
    double ym = 0;
    if( findmin(xmin, xmax, spline, acc, &xm, &ym) )
    {
        if(0){
            fprintf(stderr, "Unable to refine the maxima position\n");
            fprintf(stderr, "D = [ ");
            for(size_t kk = 0; kk<N; kk++)
            {
                fprintf(stderr, "%f, %f\n", x[kk], y[kk]);
            }
            fprintf(stderr, "]\n");
        }
        retvalue = EXIT_FAILURE;
        goto exit1;
    }
    //printf("bg=%f, xm=%f, ym=%f\n", bg, xm, ym);



    /* Find intersection .5*(max-bg) for each side */

    const gsl_root_fsolver_type * Tsolve = gsl_root_fsolver_bisection;
    gsl_root_fsolver * rsolve = gsl_root_fsolver_alloc(Tsolve);

    /* Search for left and right intersections */
    double intersections[2] = {0,0};
    /* Don't use the full line profile for the search, the edge
     * pixels are only used to find the bg level */
    double edge = (xm-(double) N2)/2.0;
    //printf("b edge = %f -> ", edge);
    for(int kk = 0; kk+1 < N2; kk++)
    {
        y[kk] < bg + 0.5*(bg-ym) ? edge = x[kk] : 0;
    }
    //printf("edge = %f\n", edge);
    int ok_left = get_intersection(xm, edge, //-N2,
                                   bg, ym,
                                   spline, acc, rsolve,
                                   &intersections[0]);
    edge = (xm+(double) N2)/2.0;
    //printf("a edge = %f -> ", edge);
    for(int kk = N+1; kk < (int) N; kk--)
    {
        y[kk] < bg + 0.5*(bg-ym) ? edge = x[kk] : 0;
    }
    //printf("edge = %f\n", edge);
    int ok_right = get_intersection(xm, edge, // N2,
                                    bg, ym,
                                    spline, acc, rsolve,
                                    &intersections[1]);
    //printf("left = %d (%f), right = %d (%f)\n", ok_left, intersections[0], ok_right, intersections[1]);
    gsl_root_fsolver_free(rsolve);

    if(ok_left == 0 && ok_right == 0)
    {
        w[0] = intersections[1] - intersections[0];
    }
    /* One-sided */
    if(ok_left == 0 && ok_right != 0)
    {
        w[0] = 2.0*(xm-intersections[0]);
    }
    if(ok_left != 0 && ok_right == 0)
    {
        w[0] = 2.0*(intersections[1] - xm);
    }
    /* Non good :() */
    if(ok_left != 0 && ok_right != 0)
    {
        w[0] = -1;
        retvalue = EXIT_FAILURE;
    }

    /* Clean up */
exit1:
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);

    /* fwhm = right - left positions */

exit2:
    return retvalue;
}


/* Create a gaussian shaped test signal, y
 * over the domain x
 * x0 is the offset from center
 * */

static void createGaussian(
    double * restrict x,
    double * restrict y,
    size_t N,
    const double x0)
{

    assert(N % 2 == 1);
    int N2 = (N-1)/2; /* middle element index */
    double s = 1.5;

    for(int kk = 0; kk< (int) N; kk++)
    {
        double p = kk-N2;
        p+=x0;
        if(kk == N2)
        {
            //printf(" [mid x = %f]\n", p);
        }
        x[kk] = kk-N2; // Set x to pixel position from center
        y[kk] = exp(-0.5*pow(p/s,2));

    }
}



/**
 */

float fwhm_lateral(fim_t * I,
                   int x, int y, int z,
                   int verbose)
{
    int nPix = 11*2+1;
    assert(nPix % 2 == 1);
    /* Line length */
    double * X = ck_calloc(nPix, sizeof(double));
    assert(X != NULL);
    for(int kk = 0; kk<nPix; kk++)
    {
        X[kk] = kk - (nPix-1)/2;
    }

    double * lx = fim_get_line_double(I, x, y, z, 0, nPix);
    double * ly = fim_get_line_double(I, x, y, z, 1, nPix);

    if(verbose > 4)
    {
        for(int kk = 0; kk<nPix; kk++)
        {
            printf("%d %f %f %f\n", kk, X[kk], lx[kk], ly[kk]);
        }
    }

    double fwhm_x = 0;
    int status_x = fwhm1d(X, lx, nPix, &fwhm_x);
    double fwhm_y = 0;
    int status_y = fwhm1d(X, ly, nPix, &fwhm_y);

    if(verbose > 3)
    {
        printf("Status x: %d, y: %d", status_x, status_y);
        printf(" fwhm x: %f, y: %f\n", fwhm_x, fwhm_y);
    }
    free(X);
    free(lx);
    free(ly);

    /** Return -1 or the average of the fwhm values that were ok */
    float fwhm = -1;

    if(status_x == 0 && status_y == 0)
    {
        fwhm = 0.5*(fwhm_x + fwhm_y);
    }

    if(status_x == 0 && status_y != 0)
    {
        fwhm = fwhm_x;
    }

    if(status_x != 0 && status_y == 0)
    {
        fwhm = fwhm_y;
    }

    return fwhm;
}

float fwhm_axial(fim_t * I,
                 int x, int y, int z,
                 int verbose)
{
    int nPix = 11*2+1;
    assert(nPix % 2 == 1);
    /* Line length */
    double * X = ck_calloc(nPix, sizeof(double));
    assert(X != NULL);
    for(int kk = 0; kk<nPix; kk++)
    {
        X[kk] = kk - (nPix-1)/2;
    }

    double * lz = fim_get_line_double(I, x, y, z,
                                      2, // dimension
                                      nPix);
    //for(int kk = 0; kk<nPix; kk++)
    //{ printf("%d %f %f\n", kk, X[kk], lz[kk]); }

    double fwhm_z = 0;
    int status_z = fwhm1d(X, lz, nPix, &fwhm_z);
    // printf("status_z = %d fwhm_z = %f\n", status_z, fwhm_z);
    free(X);
    free(lz);

    /** Return -1 or the average of the fwhm values that were ok */
    float fwhm = -1;

    if(status_z == 0 )
    {
        fwhm = fwhm_z;
    } else {
        if(verbose > 2)
        {
            printf("fwhm_axial: status=%d for (%d, %d, %d)\n", status_z, x, y, z);
        }
    }

    return fwhm;
}


fwhm_res_t * fwhm(fim_t * I,
                  const int * X, const int * Y, const int * Z,
                  const size_t N)
{
    fwhm_res_t * res = ck_calloc(N, sizeof(fwhm_res_t));
    for(size_t kk = 0; kk<N; kk++)
    {
        res[kk].fwhm_lateral = fwhm_lateral(I, X[kk], Y[kk], Z[kk], 0);
        res[kk].fwhm_axial = fwhm_axial(I, X[kk], Y[kk], Z[kk], 0);
    }
    return res;
}

static void fwhm_ut_3d(void)
{
    printf("-> fwhm_ut_3d\n");

    for(int N = 3; N < 21; N+=2)
    {
        printf("   N=%d\n", N);
        int mid = (N-1)/2;
        fim_t * I = fimt_zeros(N, N, N);
        I->V[mid + N*mid + N*N*mid] = 1;

        size_t n_points = 1;
        int X[1] = {mid};
        int Y[1] = {mid};
        int Z[1] = {mid};

        fwhm_res_t * fwhm_res = fwhm(I, X, Y, Z, n_points);
        for(size_t kk = 0; kk < n_points; kk++)
        {
            printf("(%d, %d, %d) fwhm_lateral = %f, fwhm_axial = %f\n",
                   X[kk], Y[kk], Z[kk],
                   fwhm_res[kk].fwhm_lateral,
                   fwhm_res[kk].fwhm_axial);
        }
        free(fwhm_res);
        fimt_free(I);
    }
    return;
}

int fwhm_ut(int argc, char ** argv)
{
    if(argc > 1)
    {
        printf("%s does not use any command line arguments\n", argv[0]);
    }

    fwhm_ut_3d();

    size_t N = 11; // Size of test signal
    double * y = ck_calloc(N, sizeof(double));
    double * x = ck_calloc(N, sizeof(double));
    createGaussian(x, y, N, 0);

    double w = -1; // output

    printf("-> Testing shifted signals (1D)\n");
    /* Might not be consistent due to background estimation */
    for(double delta = -1; delta <=1; delta +=0.1)
    {

        createGaussian(x, y, N, delta);
        printf("offset: % .2f pixels, ", delta);

        int status = fwhm1d(x, y, N, &w);
        if(status == EXIT_SUCCESS)
        {
            printf("fwhm: %f pixels\n", w);
        }
        if(status == EXIT_FAILURE)
        {
            printf("failed.\n");
            if(y[(N-1)/2] > y[(N-1)/2+1] && y[(N-1)/2] > y[(N-1)/2-1])
            {
                goto fail;
            }
        }
    }

    printf("-> Testing scaling x (constant y)\n");
    /* should be consistent ...  */
    for(double scale = 1; scale < 4; scale*=1.5)
    {
        createGaussian(x, y, N, 0);
        for(size_t kk = 0; kk<N; kk++)
        {
            x[kk] *= scale;
        }
        int status = fwhm1d(x, y, N, &w);
        if(status == EXIT_SUCCESS)
        {
            printf("scale = %f, fwhm = %f, fwhm/scale = %f\n", scale, w, w/scale);
        } else {
            printf("failed\n");
        }
    }


    printf("-> Testing with constant y\n");
    for(size_t kk = 0; kk<N; kk++)
    {
        y[kk] = 0;
    }
    int status = fwhm1d(x, y, N, &w);
    if(status == EXIT_SUCCESS)
    {
        goto fail;
    }

    printf("-> Previous bad case #1\n");
    double x2[] = {-11, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    double y2[] = {868, 1567, 1603, 1169, 1119, 1882, 2316, 1564, 1317, 4389, 13012, 18979, 15861, 17195, 21831, 28497, 28527, 20831, 10808, 3611, 2739, 2371, 1493};
    status = fwhm1d(x2, y2, 23, &w);
    printf("fwhm: %f\n", w);
    printf("-> Previous bad case #2\n");
    double y3[] = {1465, 1755, 1701, 2760, 10001, 14491, 9156, 2146, 4028, 10073, 18612, 18979, 16675, 18196, 19939, 19877, 16368, 12430, 10919, 8418, 4582, 1788, 1171};
    status = fwhm1d(x2, y3, 23, &w);
    printf("fwhm: %f\n", w);
    printf("-> Previous bad case #3\n");
    //double y4[] = {1298.551514, 1062.451294, 927.929016, 977.345337, 1056.960571, 2259.424805, 3914.872070, 4510.613770, 3234.024902, 13202.398438, 27393.123047, 27469.994141, 12494.097656, 5581.300781, 4689.061523, 2668.482178, 1257.371216, 1507.198364, 2094.703613, 1715.845093, 1133.830444, 958.127869, 1468.763428};
    double y4[] = {3308.149414, 3401.491211, 3017.142090, 3901.145264, 4705.533691, 6061.737305, 12203.089844, 27099.371094, 35074.617188, 31480.953125, 25067.810547, 25971.031250, 18761.736328, 13122.783203, 13089.838867, 7903.868652, 2459.835449, 894.984802, 1229.917725, 911.456909, 463.964508, 554.561157, 1196.973511};
    status = fwhm1d(x2, y4, 23, &w);
    printf("fwhm: %f\n", w);

    free(y);
    free(x);
    return(EXIT_SUCCESS);
fail:
    free(y);
    free(x);
    return(EXIT_FAILURE);
}
