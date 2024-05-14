#include "gmlfit.h"

static const int nFeatures = 10;

/* Precision for the patch representation */
typedef float pfloat; // 32-bit

/* Precision for the input image */
typedef float ifloat;

/* Optimization constants */
typedef struct {
    const pfloat * R; /* Region of interest */
    size_t M, N, P; /* Size of ROI */
} optParams;

enum features {
    f_bg = 0, /* Background level  */
    f_G = 1, /* Number of photons in signal */
    f_G0 = 2, /* Peak signal */
    f_x = 3, /* Fitted x */
    f_y = 4,
    f_z = 5,
    f_sxy = 6, /* Lateral sigma */
    f_sz = 7, /* Axial sigma */
    f_status = 8, /* 0=converged */
    f_error = 9 /* \sum -Log Likelihood over the local ROI */
};

/** Multi variate gaussian model evaluated at (0,0,0)
 * for a diagonal covariance matrix. Here specified by
 * sigmaxy, i.e. sigmax and sigmay, in the lateral plane
 * and sigmaz in the axial direction.
 */
static double mvnpdf0(double sxy, double sz)
{
    return pow(2*M_PI, -3.0/2.0)*pow( pow(sxy,4.0)*pow(sz,2.0), -1.0/2.0);
}

/** Extract a 3D region centered at D[0], D[1], D[2]
 * from V and write it in W
 * returns the number of copied pixels.
 * sets the pixel value to NaN for points outside of V
 * @param W : the location for writing.
 * @param wM, wN, wP : the size of M. The size has to be odd in all dimensions.
 * @param V : the location to read from
 * @param Vm, Vn, Vp : the size of V
 * @param D: the coordinate to center around
 * @param onebased: set to 1 for 1-indexed coordinates in D
 */

static size_t
get_roi(pfloat * restrict W,
        const size_t Wm, const size_t Wn, const size_t Wp,
        const ifloat * restrict V,
        const size_t Vm, const size_t Vn, const size_t Vp,
        const double * D, int onebased)
{

    assert(Wm % 2 == 1);
    assert(Wn % 2 == 1);
    assert(Wp % 2 == 1);

    /* Coordinate of the upper left corner to read from */
    int64_t x = nearbyint(D[0]) - (Wm-1)/2;
    int64_t y = nearbyint(D[1]) - (Wn-1)/2;
    int64_t z = nearbyint(D[2]) - (Wp-1)/2;

    if(onebased)
    {
        x--; y--; z--;
    }

    size_t writepos = 0;
    /* Loop over the ROI pixels */
    for(size_t pp = 0; pp < Wp; pp++)
    {
        int p = pp+z; // Coordinate in original volume
        for(size_t nn = 0; nn < Wn; nn++)
        {
            int n = nn+y;
            for(size_t mm = 0; mm < Wm; mm++)
            {
                int m = mm+x;

                double value = NAN;
                if(m >= 0 && m < (int) Vm)
                {
                    if( n >= 0  && n < (int) Vn)
                    {
                        if(p >= 0 && p < (int) Vp)
                        {
                            assert(m + Vm*n + Vm*Vn*p < Vm*Vn*Vp);
                            value = V[m + Vm*n + Vm*Vn*p];
                        }
                    }
                }

                W[writepos++] = value;
            }
        }
    }

    return writepos;
}



double error_fun(const gsl_vector * v, void * _params)
{
    /* Check background value */
    if(gsl_vector_get(v, 0) < 0)
    {
        return 1e99;
    }
    /* Check photon count*/
    if(gsl_vector_get(v, 1) < 0)
    {
        return 1e99;
    }
    /* Check size of signal */
    if(gsl_vector_get(v, 5) < 0)
    {
        return 1e99;
    }
    if(gsl_vector_get(v, 6) < 0)
    {
        return 1e99;
    }

    const optParams * restrict params = (optParams*) _params;
    const pfloat * restrict R = params->R; /* Image data */
    const size_t M = params->M;
    const size_t N = params->N;
    const size_t P = params->P;


    double mx = gsl_vector_get(v, 2);
    double my = gsl_vector_get(v, 3);
    double mz = gsl_vector_get(v, 4);
    double sxy = gsl_vector_get(v, 5);
    double sz = gsl_vector_get(v, 6);

    double G0 = mvnpdf0(sxy, sz);

    double E = 0;
    for(size_t pp = 0; pp < P; pp++)
    {

        double z = (double) pp - (P-1)/2;
        for (size_t nn = 0; nn < N; nn++)
        {
            double y = (double) nn - (N-1)/2;
            for(size_t mm = 0; mm < M; mm++)
            {
                double x = (double) mm - (M-1)/2;

                double Pixel = R[pp*M*N + nn*M + mm];
                if(isnan(Pixel))
                {
                    continue;
                }

                /* Gaussian value */
                double G = G0*exp( - 0.5*(
                                         pow( (x-mx)/sxy, 2)
                                       + pow( (y-my)/sxy, 2)
                                       + pow( (z-mz)/sz, 2) ) );
                /* Full model */
                double Model = gsl_vector_get(v, 0) + gsl_vector_get(v,1)*G;

                E += Pixel*log(Model) - Model;
                assert(isnan(E) == 0);
            }
        }
    }

    return -E;
}

void
error_fun_fdf (const gsl_vector *x, void *params,
               double *f, gsl_vector *df)
{
    double e0 = error_fun(x, params);
    double delta = 1e-6;
    gsl_vector * x_delta = gsl_vector_alloc(x->size);
    for(size_t kk = 0; kk < df->size; kk++)
    {
        gsl_vector_memcpy(x_delta, x);
        gsl_vector_set(x_delta, kk, gsl_vector_get(x, kk) + delta);
        double e = error_fun(x_delta, params);
        gsl_vector_set(df, kk, (e-e0)/delta);
        //printf("%zu : %e \n", kk, gsl_vector_get(df, kk));
    }
    gsl_vector_free(x_delta);
    *f = e0;
    return;
}

void
error_fun_df (const gsl_vector *x, void *params,
              gsl_vector *df)
{
    double e = 0;
    error_fun_fdf(x, params, &e, df);
    return;
}


/** Localize a dot centered in V (or at least approximately)
 * Model: b + c*G(mu, sigma)
 */

static int
localize_dot(const pfloat * restrict V,
             const size_t M, const size_t N, const size_t P,
             const double sigma_xy, const double sigma_z,
             double * F, size_t max_iterations) /* Output */
{

    /* Parameters to the error functional */
    optParams par;
    par.R = V;
    par.M = M;
    par.N = N;
    par.P = P;

    double bg = 0;
    {
        size_t pos = 0;
        while(isnan(V[pos]))
        {
            pos++;
        }
        bg = V[pos];
    }

    //printf("ROI: %zu x %zu x %zu\n", M, N, P);
    double nphot = 0;
    size_t nok = 0;
    for(size_t kk = 0; kk<M*N*P; kk++)
    {
        //printf("%f ", V[kk]);
        if( !isnan(V[kk]) )
        {
            nok++;
            V[kk] < bg ? bg = V[kk] : 0;
            nphot += V[kk];
        }
    }
    //printf("\n");

    //printf("nok = %zu\n", nok);
    //printf("nphot: %f\n", nphot);
    assert(nok > 0);
    assert(nok <= M*N*P);
    nphot -= bg*(double) nok;
    //printf("bg: %f, nphot: %f\n", bg, nphot);

    /* Starting point for parameters to be optimized */
    gsl_vector * x = gsl_vector_alloc(7);
    gsl_vector_set(x, 0, bg);
    gsl_vector_set(x, 1, nphot);
    gsl_vector_set(x, 2, 0); // x
    gsl_vector_set(x, 3, 0); // y
    gsl_vector_set(x, 4, 0); // z
    gsl_vector_set(x, 5, sigma_xy);
    gsl_vector_set(x, 6, sigma_z);


    /* Initialize method and iterate */


    gsl_multimin_function_fdf minex_func_fdf;
    minex_func_fdf.n = 7;
    minex_func_fdf.f = &error_fun;
    minex_func_fdf.df = &error_fun_df;
    minex_func_fdf.fdf = &error_fun_fdf;
    minex_func_fdf.params = &par;

    /* Contestants: conjugate_fr, conjugate_pr, vector_bfgs2
     * nmsimplex2: 0.9 s
     * steepest_descent: 0.7 s
     * conjugate_fr: 1.7 s
     * conjugate_pr: 2.0 s
     * vector_bfgs2: 0.3 s
     */

    /* Note: Setting tol to a very low value results in very poor
    * results, could be because the optimization is more prone to get
    * stuck in local minima. The value 0.1 is the recommended according to the
    * GSL documentation.
    * 0.01 is faster but leads to worse results.
    * 0.5 is slower but results are not better.
    */
    const gsl_multimin_fdfminimizer_type * T = gsl_multimin_fdfminimizer_vector_bfgs2;
    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc (T, 7);
    gsl_multimin_fdfminimizer_set(s, // Optimizer
                                  &minex_func_fdf, // Function to optimize
                                  x,// Start point
                                  0.01, // step size
                                  0.1); // line search tol

    int status = 0;
    size_t iter = 0;
    /* Possibly fewer iterations for really bright signals
     * and a few more for very noisy data */

    do
    {
        iter++;
        status = gsl_multimin_fdfminimizer_iterate(s);

        if (status)
        {
            // printf("Stop beacause of %s\n", gsl_strerror (status));
            break;
        }

        status =  gsl_multimin_test_gradient(gsl_multimin_fdfminimizer_gradient(s), 1e-5);

        if (status == GSL_SUCCESS)
        {
            //printf("Stop beacause of convergence 2\n");
            break;
        }
        assert(status == GSL_CONTINUE);

    } while (status == GSL_CONTINUE && iter < max_iterations);

    //printf("iter = %zu\n", iter);

    F[f_bg] = gsl_vector_get(s->x, 0);
    F[f_G] = gsl_vector_get(s->x, 1);
    F[f_x] = gsl_vector_get(s->x, 2);
    F[f_y] = gsl_vector_get(s->x, 3);
    F[f_z] = gsl_vector_get(s->x, 4);
    F[f_sxy] = gsl_vector_get(s->x, 5);
    F[f_sz] = gsl_vector_get(s->x, 6);


    F[f_G0] = F[f_G]*mvnpdf0(F[f_sxy], F[f_sz]);
    F[f_status] = status;
    F[f_error] = gsl_multimin_fdfminimizer_minimum(s);

    gsl_vector_free(x);
    gsl_multimin_fdfminimizer_free(s);

    return status;
}


gmlfit * gmlfit_new(void)
{
    gmlfit * conf = calloc(1, sizeof(gmlfit));
    assert(conf != NULL);

    conf->max_iter = 1000;
    conf->verbose = 1;
    return conf;
}

/** Fit the dots specified by (x,y,z) in X.
 *  For staring guess use the size from sigma_xy and sigma_z
 *
 */

double *
gauss_mlfit3(const ifloat * restrict V,
             const size_t M, const size_t N, const size_t P,
             const double sigma_xy, const double sigma_z,
             const double * X, const size_t nX, int onebased)
{
    gmlfit * conf = gmlfit_new();
    conf->image = V;
    conf->M = M;
    conf->N = N;
    conf->P = P;
    conf->sigma_xy = sigma_xy;
    conf->sigma_z = sigma_z;
    conf->X = X;
    conf->nX = nX;
    conf->x_base = onebased;
    conf->log = stdout;
    conf->verbose = 1;
    double * F = gmlfit_run(conf);
    free(conf);
    return F;
}

static int gmlfit_validate(gmlfit * conf)
{
    int status = 0;
    if(conf->nX == 0)
    {
        if(conf->log != NULL)
        {
            fprintf(conf->log, "No points to fit (nX = =0)\n");
        }
        status = 1;
    }

    if(conf->image == NULL)
    {
        if(conf->log != NULL)
        {
            fprintf(conf->log, "Image missing (image == NULL)\n");

        }
        status = 1;
    }

    if(conf->X == NULL)
    {
        if(conf->log != NULL)
        {
            fprintf(conf->log, "Start coordinates missing (X == NULL)\n");
        }
        status = 1;
    }

    return status;
}

static void gmlfit_print(gmlfit * conf)
{
    if(conf->log == NULL)
    { return; }

    if(conf->verbose == 0)
    { return; }

    fprintf(conf->log, "gmlfit settings:\n");
    fprintf(conf->log, "   Image size: %zu x %zu x %zu\n",
            conf->M, conf->N, conf->P);
    fprintf(conf->log, "   Points to fit: %zu\n", conf->nX);
    fprintf(conf->log, "   Size guess: %f, %f\n", conf->sigma_xy, conf->sigma_z);
    fprintf(conf->log, "   max_iter: %zu\n", conf->max_iter);
    if(conf->num_threads > 0)
    {
        fprintf(conf->log, "   num_threads = %d\n", conf->num_threads);
    } else {
        fprintf(conf->log, "    num_threads: not set\n");
    }
    fprintf(conf->log, "\n");
    return;
}

double * gmlfit_run(gmlfit * conf)
{
    if(gmlfit_validate(conf))
    {
        return NULL;
    }

    if(conf->num_threads > 0)
    {
        omp_set_num_threads(conf->num_threads);
    }

    gmlfit_print(conf);


/* Allocate memory for output table */
    double * F = calloc(conf->nX*nFeatures, sizeof(double));
    assert(F != NULL);

    /** Determine a good windows size
     * Currently 4*sigma pixels in each dimension, and at least 5^3 pixels.
     * If the window is too small:
     * - More noisy estimates (less data to base the fitting on)
     * - hard to guess the initial background value.
     * Too large:
     * - More biased by nearby dots
     * - Slower
     *
     * The setting 4 was found to be better than 3 and 5 on synthetic
     * Gaussian dots.
     */
    size_t wMN = ceil(4*conf->sigma_xy);
    wMN < 5 ? wMN = 5 : 0;
    wMN % 2 == 0 ? wMN++ : 0;

    size_t wP = ceil(4*conf->sigma_z);
    wP < 5 ? wP = 5: 0;
    wP % 2 == 0 ? wP++ : 0;

    if(conf->verbose > 1)
    {
        if(conf->log != NULL)
        {
            printf("Using a ROI size of %zu x %zu x %zu pixels\n",
                   wMN, wMN, wP);
        }
    }

#pragma omp parallel
    {
        /* Allocate memory for a local patch */
        pfloat * R = calloc(wMN*wMN*wP, sizeof(pfloat));
        assert(R != NULL);

#pragma omp for schedule(dynamic)
        for(size_t kk = 0; kk < conf->nX; kk++)
        {
            /* Copy the neighbourhood around each D into W
               This fails if the dot is completely outside of the image
            */
            if(get_roi(R, wMN, wMN, wP,
                       conf->image, conf->M, conf->N, conf->P,
                       conf->X+kk*3, conf->x_base) > 0)
            {
                //printf("Processing dot %zu / %zu\n", kk+1, nX);
                localize_dot(R, wMN, wMN, wP,
                             conf->sigma_xy, conf->sigma_z,
                             F + nFeatures*kk, conf->max_iter);
            } else {
                for(int ff = 0; ff < nFeatures; ff++)
                {
                    F[ kk*nFeatures + ff ] = 0;
                }
            }
        }
        free(R);
    }

    return F;
}

int gmlfit_ut(int argc, char ** argv)
{
    size_t M = 100, N = 100, P = 100;
    ifloat * I = malloc(M*N*P*sizeof(ifloat));
    assert(I != NULL);
    for(size_t kk =0 ; kk < M*N*P; kk++)
    {
        I[kk] = rand();
    }

    size_t nX = 2;
    if(argc > 1)
    {
        nX = atol(argv[1]);
    }

    double * X = malloc(3*nX*sizeof(double));
    assert(X!=NULL);
    for(size_t kk = 0; kk < nX; kk++)
    {
        double * P = X + 3*kk;
        P[0] = 20; P[1] = 20; P[2] = 20;
    }
    double sigma_xy = 3.51;
    double sigma_z = 3.16;

    gmlfit * conf = gmlfit_new();
    conf->image = I;
    conf->M = M;
    conf->N = N;
    conf->P = P;
    conf->sigma_xy = sigma_xy;
    conf->sigma_z = sigma_z;
    conf->X = X;
    conf->nX = nX;
    conf->log = stdout;
    conf->verbose = 2;


    double * F = gmlfit_run(conf);
    free(conf);

    assert(F!=NULL);
    size_t nX_show = 5;
    nX_show > nX ? nX_show = nX : 0;
    if(nX > nX_show)
    {
        printf("Showing results for %zu / %zu dots\n", nX_show, nX);
    }
    for(size_t kk = 0; kk < nX_show; kk++)
    {
        for(int pp = 0; pp<nFeatures; pp++)
        {
            printf("%f ", F[nFeatures*kk + pp]);
        }
        printf("\n");
    }
    free(F);

    free(I);
    free(X);
    return EXIT_SUCCESS;
}
