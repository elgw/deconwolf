/* (C) 2020 Erik L. G. Wernersson
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

#define _USE_MATH_DEFINES
#include <math.h>
#include <assert.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include <wchar.h>
#include <locale.h>
#include <gsl/gsl_integration.h>
#include <omp.h>

#include "dw_util.h"
#include "fim_tiff.h"
#include "dw_version.h"
#include "lanczos.h"
#include "li.h"
#include "bw_gsl.h"

#include "dw_bwpsf.h"
#include "npio.h"

/* Windows does not agree on the definition of complex numbers,
 * we will call them dcomplex regardless of the underlying library
 */
#ifdef WINDOWS
typedef _Dcomplex dcomplex;
#else
typedef double complex dcomplex;
#endif

int GLOB_N_GSL_EROUND = 0; /* Counter for GSL_EROUND */

static void BW_slice(float * , float z, bw_conf * conf);
static void bw_argparsing(int , char ** , bw_conf * conf);

// Get the command line options
static void getCmdLine(int argc, char ** argv, bw_conf * s);


/* Integration over pixel */
typedef struct {
    bw_conf * conf;
    double * radprofile;
    size_t nr;
    int radsample;
    double y0;
    double y1;
    double x;
    gsl_integration_workspace * wspx;
    gsl_integration_workspace * wspy;
} pixely_t;


static void dw_bw_gsl_err_handler(const char * reason,
                           const char * file,
                           int line,
                           int gsl_errno)
{
    /* The model oscillate heavily when far from the origin.
     * We use high precision when possible and accept that the
     * numerical integration sometimes can't reach that.
     */

    if(gsl_errno == GSL_EROUND)
    {
        GLOB_N_GSL_EROUND++;
        return;
    }

    printf("An unexpected error happened and dw_bw can't continue. Please file a bug report\n");
    printf("Reason: %s\n", reason);
    printf("File: %s\n", file);
    printf("Line: %d\n", line);
    printf("gsl_errno: %d\n", gsl_errno);
    exit(EXIT_FAILURE);
}

void bw_conf_printf(FILE * out, bw_conf * conf)
{
    if(out != stdout)
    {
        fprintf(out, "cmd: %s\n", conf->cmd);
    }
    fprintf(out, "lambda = %.2f nm\n", conf->lambda);
    fprintf(out, "NA = %f\n", conf->NA);
    fprintf(out, "ni = %f\n", conf->ni);
    fprintf(out, "resLateral = %.2f nm\n", conf->resLateral);
    fprintf(out, "resAxial = %.2f nm\n", conf->resAxial);
    fprintf(out, "out size: [%d x %d x %d] pixels\n", conf->M, conf->N, conf->P);
    fprintf(out, "nThreads: %d\n", conf->nThreads);
    fprintf(out, "Verbosity: %d\n", conf->verbose);
    fprintf(out, "Overwrite: %d\n", conf->overwrite);
    fprintf(out, "File: %s\n", conf->outFile);
    fprintf(out, "Log: %s\n", conf->logFile);

    fprintf(out, "Oversampling of radial profile: %d X\n", conf->oversampling_R);

    if(conf->mode_bw == MODE_BW_LI)
    {
        fprintf(out, "BW Integral: Li's fast method enabled\n");
    }
    if(conf->mode_bw == MODE_BW_GSL)
    {
        fprintf(out, "BW Integral: gsl_integration_qag\n");
    }
    fprintf(out, "epsabs: %e\n", conf->epsabs);
    fprintf(out, "epsrel: %e\n", conf->epsrel);

    // Show theoretical FWHM of the PSF
    // And look at the number of pixels per fwhm
    double fwhm_r = 1.616340*conf->lambda/M_PI/conf->NA;
    double fwhm_z = 2*2.783115/M_PI*conf->ni/pow(conf->NA, 2)*conf->lambda;

    fprintf(out, "FWHM_r (lateral plane) %.2f nm\n", fwhm_r);
    fprintf(out, "Rayleigh Resolution (lateral) = %.2f nm\n",
            0.61*conf->lambda / conf->NA);
    double qlateral = fwhm_r/conf->resLateral;
    fprintf(out, "FWHM_r / dr = %.2f\n", qlateral);
    if(qlateral < 2)
    {
        fprintf(out,"! Warning: This is suboptimal, aim for at least 2\n");
    }

    fprintf(out, "FWHM_z (axial direction) %.2f nm\n", fwhm_z);
    double qaxial = fwhm_z/conf->resAxial;
    fprintf(out, "FWHM_z / dz = %.2f\n", qaxial);
    if(qaxial < 2)
    {
        fprintf(out,"! Warning: This is suboptimal, aim for at least 2\n");
    }
    return;
}


bw_conf * bw_conf_new()
{
    bw_conf * conf = calloc(1, sizeof(bw_conf));
    assert(conf != NULL);

    /* Physical settings */
    conf->lambda = 600;
    conf->NA = 1.45;
    conf->ni = 1.515;
    conf->resAxial = 300;
    conf->resLateral = 130;

    /* Output image Settings*/
    conf->M = 181;
    conf->N = 181;
    conf->P = 181; // Does only need to be 2*size(I,3)+1

    conf->nThreads = 8;
    conf->V = NULL;
    conf->verbose = 1;
    conf->V = NULL;
    conf->outFile = NULL;
    conf->logFile = NULL;
    conf->log = NULL;
    conf->overwrite = 0;
    conf->testing = 0;

    /* BW integration */
    conf->oversampling_R = 17; /* Expose from CLI? */
    conf->mode_bw = MODE_BW_GSL;

    /* For XY integration */
    conf->epsabs = 1e-6; /* As good as we can get without any GSL_EROUND */
    conf->epsrel = 1e-6;
    conf->key = 1;
    /* For XY integration and BW integration */
    conf->limit = 10000000;

    return conf;
}

void bw_conf_free(bw_conf ** _conf)
{
    bw_conf * conf = *_conf;
    if(conf->log != NULL)
    {
        fclose(conf->log);
    }
    if(conf->ftif != NULL)
    {
        fim_tiff_destroy(conf->ftif);
        conf->ftif = NULL;
    }
    free(conf->outFile);
    free(conf->logFile);
    free(conf->cmd);
    free(conf);
    _conf = NULL;
    return;
}


static void getCmdLine(int argc, char ** argv, bw_conf * s)
{
    // Copy the command line to s->cmd
    int lcmd=0;
    for(int kk = 0; kk<argc; kk++)
    {
        lcmd += strlen(argv[kk]);
    }
    lcmd += argc+2;
    s->cmd = malloc(lcmd);
    assert(s->cmd != NULL);
    int pos = 0;
    for(int kk = 0; kk<argc; kk++)
    {
        sprintf(s->cmd+pos, "%s ", argv[kk]);
        pos += strlen(argv[kk])+1;
    }
    s->cmd[pos-1] = '\0';
}

static void fprint_time(FILE * f)
{
    f == NULL ? f = stdout : 0;
    time_t now = time(NULL);
    char * tstring = ctime(&now);
    fprintf(f, "%s\n", tstring);
}

static void bw_version(FILE* f)
{
#ifdef WINDOWS
    fprintf(f, "deconwolf: '%s'\n", deconwolf_version);
#else
    fprintf(f, "deconwolf: '%s' PID: %d\n", deconwolf_version, (int) getpid());
#endif
}

static void usage(__attribute__((unused)) int argc, char ** argv, bw_conf * s)
{
    printf("Usage:\n");
    printf("\t%s <options> psf.tif \n", argv[0]);
    printf("\n");
    printf("Required parameters:\n");
    printf(" --NA\n\t Set numerical aperture.\n");
    printf(" --lambda\n\t Set emission wavelength [nm].\n");
    printf(" --ni\n\t Set refractive index for the immersion (also used for coverslip and sample).\n");
    printf(" --resxy\n\t Set lateral pixel size (x-y) [nm].\n");
    printf(" --resz\n\t Set the axial pixel size in (z) [nm].\n");
    printf("\n");
    printf("Optional:\n");
    printf(" --size N \n\t Set output size to N x N x P [pixels].\n"
           "\t This should be set big enough so that the PSF cone\n"
           "\t isn't cropped at the first and last slice of the PSF.\n");
    printf(" --nslice P\n\t Set output size to N x N x P [pixels].\n");
    printf("\t P has to be an odd number. Please not that P should be at\n"
           "\t least (2Z-1) where Z is the number of slices in the image\n"
           "\t to be deconvolved\n");
    printf(" --threads\n\t Set number of threads. Currently set to %d\n", s->nThreads);
    printf("\t Q has to be an odd number.\n");

    printf(" --li\n\t Use Li's method (our implementation) for the BW integral.\n");
    printf("\t About 2x faster for the default PSF size. "
           "Not thoroughly tested.");
    if(s->mode_bw == MODE_BW_LI)
    {
        printf(" Default.");
    }
    printf("\n");
    printf(" --gsl\n\t Use GSL for the BW integral as well as the integration over pixels.");
    if(s->mode_bw == MODE_BW_GSL)
    {
        printf(" Default.");
    }
    printf("\n");
    printf(" --epsabs \n\t Set absolute error tolerance for the integration. Default: %e\n",
           s->epsabs);
    printf(" --epsrel \n\t Set relative error tolerance for the integration. Default: %e\n",
           s->epsrel);

    printf(" --overwrite \n\t Overwrite the target image if it exists.\n");
    printf(" --version\n\t Show version info and quit.\n");
    printf(" --help\n\t Show this message.\n");
    printf(" --verbose L\n\t Set verbosity level to L.\n"
           "\t where 0 = quiet, 2 = full\n");
    //   printf(" --test\n\t Run some tests\n");
    printf("\n");
    printf("Web page: https://www.github.com/elgw/deconwolf/\n");
    return;
}

/* Returns 0 on success */
int bw_conf_validate(bw_conf * s)
{
    if(s->lambda < 50)
    {
        fprintf(stderr, "Error: lambda has be be at least 50 nm\n");
        return EXIT_FAILURE;
    }

    if(s->lambda > 2000)
    {
        fprintf(stderr, "Error lambda can be at most 2000 nm\n");
        return EXIT_FAILURE;
    }

    if(s->M % 2 == 0)
    {
        if(s->verbose > 0)
        {
            printf("Warning: The size has to be odd, 1, 3, ...\n");
            printf("Changing from %d to %d.\n", s->M, s->M+1);
        }
        s->M++;
        s->N++;
    }

    if(s->P % 2 == 0)
    {
        if(s->verbose > 0)
        {
            printf("Warning: The number of slices has to be odd, 1, 3, ...\n");
            printf("          Changing from %d to %d.\n", s->P, s->P+1);
        }
        s->P++;
    }

    /* Not more threads than slices */
    if(s->nThreads > s->P)
    {
        s->nThreads = s->P;
    }
    return EXIT_SUCCESS;
}

static void bw_argparsing(int argc, char ** argv, bw_conf * s)
{
    bw_conf * defaults = calloc(1, sizeof(bw_conf));
    assert(defaults != NULL);
    memcpy(defaults, s, sizeof(bw_conf));

    if(argc < 2)
    {
        printf("Deconwolf PSF generator: Too few input arguments.\n");
        printf("See `%s --help` or 'man dw_bw'\n", argv[0]);
        exit(-1);
    }

    getCmdLine(argc, argv, s);

    int gotna = 0;
    int gotlambda = 0;
    int gotni = 0;
    int gotdx = 0;
    int gotdz = 0;

    struct option longopts[] =
        {
            { "epsabs",     required_argument, NULL,   'a'},
            { "li",         no_argument,       NULL,   'f'},
            { "gsl",        no_argument,       NULL,   'g'},
            { "help",       no_argument,       NULL,   'h'},
            { "ni",         required_argument, NULL,   'i'},
            { "lambda",     required_argument, NULL,   'l'},
            { "NA",         required_argument, NULL,   'n'},
            { "verbose",    required_argument, NULL,   'p'},
            { "epsrel",     required_argument, NULL,   'r'},
            { "threads",    required_argument, NULL,   't'},
            { "version",    no_argument,       NULL,   'v'},
            { "overwrite",  no_argument,       NULL,   'w'},
            { "resxy",      required_argument, NULL,   'x'},
            { "resz",       required_argument, NULL,   'z'},
            { "size",       required_argument, NULL,   'N'},
            { "nslice",     required_argument, NULL,   'P'},
            { "testing",    no_argument,       NULL,   'T'},
            { NULL,         0,                 NULL,    0 }
        };
    int Pset = 0;
    int ch;

    while((ch = getopt_long(argc, argv, "a:fghi:l:n:p:r:t:vwx:z:N:P:T", longopts, NULL)) != -1)
    {
        switch(ch) {
        case 'a':
            s->epsabs = atof(optarg);
            break;
        case 'f':
            s->mode_bw = MODE_BW_LI;
            break;
        case 'g':
            s->mode_bw = MODE_BW_GSL;
            break;
        case 'h':
            usage(argc, argv, defaults);
            exit(0);
        case 'i':
            s->ni = atof(optarg);
            gotni = 1;
            break;
        case 'l':
            s->lambda = atof(optarg);
            gotlambda = 1;
            break;
        case 'n':
            s->NA = atof(optarg);
            gotna = 1;
            break;
        case 'p':
            s->verbose = atoi(optarg);
            break;
        case 'r':
            s->epsrel = atof(optarg);
            break;
        case 't':
            s->nThreads = atoi(optarg);
            break;
        case 'v':
            bw_version(stdout);
            exit(0);
            break;
        case 'w':
            s->overwrite = 1;
            break;
        case 'x':
            s->resLateral = atof(optarg);
            gotdx = 1;
            break;
        case 'z':
            s->resAxial = atof(optarg);
            gotdz = 1;
            break;
        case 'N':
            s->M = atoi(optarg);
            s->N = s->M;
            break;
        case 'P':
            s->P = atoi(optarg);
            Pset = 1;
            break;
        case 'T':
            s->testing = 1;
            break;
        default:
            exit(EXIT_FAILURE);
            break;
        }
    }

    free(defaults);

    int stay = 1;
    if(gotna == 0)
    {
        printf("Numerical Aperture (--NA) not specified\n");
        stay = 0;
    }
    if(gotlambda == 0)
    {
        printf("Wave Length (--lambda) not specified\n");
        stay = 0;
    }
    if(gotni == 0)
    {
        printf("Refractive index (--ni) not specified\n");
        stay = 0;
    }
    if(gotdx == 0)
    {
        printf("Laterial pixel size not specified (--resxy) not specified\n");
        stay = 0;
    }
    if(gotdz == 0)
    {
        printf("Axial distance between planes not specified (--resz) not specified\n");
        stay = 0;
    }
    if(stay == 0)
    {
        exit(EXIT_FAILURE);
    }


    if(Pset == 0)
    {
        // Assume that samples are imaged with about the same thickness
        // regardless of resolution
        int nslice = floor(181.0 * 300.0/(s->resAxial)/2.0);
        s->P = nslice*2+3;

        printf("Setting the number of slices to %d, based on the pixel size.\n", s->P);
        printf("     Please adjust manually with the --nslice P argument.\n");
    }

    if(optind + 1 == argc)
    {
        s->outFile = malloc(strlen(argv[argc-1]) + 1);
        assert(s->outFile != NULL);
        sprintf(s->outFile, "%s", argv[argc-1]);
    }

    if(s->outFile == NULL)
    {
        s->outFile = malloc(100*sizeof(char));
        assert(s->outFile != NULL);
        sprintf(s->outFile, "PSFBW_%.2f_%.2f_%.1f_%.1f_%.1f.tif",
                s->NA, s->ni, s->lambda, s->resLateral, s->resAxial);
    }


    s->logFile = malloc(strlen(s->outFile) + 10);
    assert(s->logFile != NULL);
    sprintf(s->logFile, "%s.log.txt", s->outFile);


    return;
}

static double cabs2(dcomplex v)
{
    return pow(creal(v), 2) + pow(cimag(v), 2);
}



void BW(bw_conf * conf)
{
    assert(conf->V == NULL);
    assert(conf->nThreads > 0);
    conf->V = malloc(conf->M*conf->N*conf->P*sizeof(float));
    assert(conf->V != NULL);

#pragma omp parallel for
    for(int z = 0; z < (conf->P+1)/2; z++)
    {
        float defocus = conf->resAxial * (z - (conf->P - 1.0) / 2.0);
        BW_slice(conf->V + z*conf->M*conf->N, defocus, conf);
    }

    /* symmetry in Z */
    for(int z = 0; z< (conf->P-1)/2; z++)
    {
        memcpy(conf->V+(conf->P-z-1)*conf->M*conf->N,
               conf->V+z*conf->M*conf->N,
               conf->M*conf->N*sizeof(float));
    }
    return;
}

static double pixelpointfun(double y, void * _conf)
{
    /* interpolate radial profile at (x,y) */
    pixely_t * iconf = (pixely_t *) _conf;
    double r2 = pow(iconf->x, 2) + pow(y, 2);
    double r = 0;
    if(fabs(r2) > 1e-10)
    {
        r = sqrt(r2);
    }

    double res = lanczos5(iconf->radprofile,
                          iconf->nr,
                          r*(double) iconf->radsample);
    return res;
}

static double pixelyfun(double x, void *_conf)
{
    pixely_t * iconf = (pixely_t *) _conf;
    //printf("pixelyfun: iconf->nr = %zu, iconf->radsample=%d\n", iconf->nr, iconf->radsample);
    iconf->x = x;

    gsl_function fun;
    fun.function = &pixelpointfun;
    fun.params = iconf;
    double result = 0;
    double abserr = 0;
    //printf("y0=%f, y1=%f\n", iconf->y0, iconf->y1);
    gsl_integration_qag(&fun, iconf->y0, iconf->y1,
                        iconf->conf->epsabs, iconf->conf->epsrel,
                        iconf->conf->limit, iconf->conf->key, iconf->wspy,
                        &result, &abserr);

    return result;
}

static double integrate_pixel(bw_conf * conf,
                       gsl_integration_workspace * wspx,
                       gsl_integration_workspace * wspy,
                       double * radprofile, size_t nr, int radsample,
                       double x0, double x1,
                       double y0, double y1)
{
    assert(x1 > x0);
    assert(y1 > y0);

    gsl_function fun;
    fun.function = &pixelyfun;

    pixely_t iconf;
    iconf.radprofile = radprofile;
    iconf.nr = nr;
    iconf.radsample = radsample;
    iconf.y0 = y0;
    iconf.y1 = y1;
    iconf.conf = conf;
    iconf.wspx = wspx;
    iconf.wspy = wspy;
    //printf("intgrate_pixel: iconf->nr = %zu, iconf->radsample=%d\n", iconf.nr, iconf.radsample);

    fun.params = &iconf;

    double result = 0;
    double abserr = 0;

    //printf("x0=%f, x1=%f, y0=%f, y1=%f, epsabs=%e, epsrel=%e key=%d\n", x0, x1, y0, y1, conf->epsabs, conf->epsrel, conf->key);

    gsl_integration_qag(&fun, x0, x1,
                        conf->epsabs, conf->epsrel,
                        conf->limit, conf->key, wspx,
                        &result, &abserr);

    return result/((x1-x0)*(y1-y0));
}

static void BW_slice(float * V, float z, bw_conf * conf)
{
    /* Calculate the PSF for a single plane/slice
     * V is a pointer to the _slice_ not the whole volume
     * This function is called from different threads with
     * a copy of the initial conf
     */

    gsl_integration_workspace * wspx =
        gsl_integration_workspace_alloc(conf->limit);
    gsl_integration_workspace * wspy =
        gsl_integration_workspace_alloc(conf->limit);


    /* The center of the image in units of [pixels] */
    double x0 = ((double) conf->M - 1.0) / 2.0;
    double y0 = ((double) conf->N - 1.0) / 2.0;

    /* The bw integral will be sampled up to radmax,
     * with radnsamples samples per pixel
     */
    int radmax = (int) round(sqrt(pow(conf->M - x0, 2) + pow(conf->N - y0, 2))) + 1;
    int radnsamples = conf->oversampling_R;

    size_t nr = (radmax+1)*conf->oversampling_R+2;
    double * r = malloc(nr * sizeof(double));
    assert(r != NULL);

    // abs(x)^2
    double * radprofile = malloc(nr * sizeof(double));
    assert(radprofile != NULL);

    /* Calculate the radial profile
       possibly by integrating over Z */

    for(size_t n = 0; n < nr; n++)
    {
        r[n] = ((double) n) / ((double) radnsamples);
    }

    /* 1D integral using Bessel series, li */
    if(conf->mode_bw == MODE_BW_LI)
    {
        li_conf * L = li_new(z);
        L->lambda = conf->lambda;
        L->NA = conf->NA;
        L->ni = conf->ni;
        for(size_t n = 0; n < nr; n++) // over r
        {
            dcomplex v = li_calc(L, r[n]*conf->resLateral);
            radprofile[n] = cabs2(v);
        }
        li_free(&L);
    }


    /* BW integral by GSL */
    if(conf->mode_bw == MODE_BW_GSL)
    {
        bw_gsl_conf_t * bw_gsl_conf = bw_gsl_new(conf->limit);
        bw_gsl_conf->NA = conf->NA;
        bw_gsl_conf->ni = conf->ni;
        bw_gsl_conf->epsabs = conf->epsabs;
        bw_gsl_conf->epsrel = conf->epsrel;
        bw_gsl_conf->lambda = conf->lambda;


        if(z == 0 && (conf->log != NULL) ) /* Only write from the thread that processes z==0 */
        {
            fprintf(conf->log, "Settings for BW integration over Z:\n");
            bw_gsl_fprint(conf->log, bw_gsl_conf);
        }


        /* Integrate bw for a equidistant set of radii */
        for(size_t n = 0; n < nr; n++)
        {
            double radius = r[n]*conf->resLateral;
            //double value  = bw_gsl_integrate_z(bw_gsl_conf, radius, z0, z1);
            double value = bw_gsl_integrate(bw_gsl_conf, radius, z);
            radprofile[n] = value;
        }
        bw_gsl_free(bw_gsl_conf);
    }
    assert(r[0] == 0);

    /*
     * For each pixel in the quadrant integrate over x and y
     */

    //printf("!\n");
    for (int x = 0; 2*x <= conf->M; x++) {
        for (int y = x; 2*y <= conf->N; y++) {

            /* radius of the current pixel in units of [pixels] */
            double rPixel = sqrt( pow((double) x - x0, 2)
                                  + pow((double) y - y0, 2));


            /* Pixel integration */
            double pIntensity = 0;
            if(rPixel < radmax)
            {
                //printf("x=%d, y=%d, [%f, %f]x[%f, %f]x%f\n", x, y,
                //       (double) x - 0.5 - x0, (double) x + 0.5 - x0,
                //       (double) y - 0.5 - y0, (double) y + 0.5 - y0,z );

                pIntensity = integrate_pixel(conf,
                                             wspx, wspy,
                                             radprofile, nr, radnsamples,
                                             (double) x - 0.5 - x0, (double) x + 0.5 - x0,
                                             (double) y - 0.5 - y0, (double) y + 0.5 - y0);
            }

            /* Use symmetry to fill the output image plane */
            int xf = conf->M-x-1;
            int yf = conf->N-y-1;
            V[x + conf->M * y] = pIntensity;
            V[x + conf->M * yf] = pIntensity;

            V[y + conf->M * x] = pIntensity;
            V[y + conf->M * xf] = pIntensity;

            V[xf + conf->M * y] = pIntensity;
            V[xf + conf->M * yf] = pIntensity;

            V[yf + conf->M * x] = pIntensity;
            V[yf + conf->M * xf] = pIntensity;
        }
    }

    gsl_integration_workspace_free(wspx);
    gsl_integration_workspace_free(wspy);
    free(r);
    free(radprofile);
    return;
}

static void unit_tests(bw_conf * conf)
{
    printf("No unit tests to run\n");
    bw_conf_printf(stdout, conf);
    return;
}

int dw_bwpsf(int argc, char ** argv)
{
    struct timespec tstart, tend;
    dw_gettime(&tstart);

    bw_conf * conf = bw_conf_new();
    bw_argparsing(argc, argv, conf);

    if(bw_conf_validate(conf))
    {
        printf("Can't continue with the provided settings\n");
        return EXIT_FAILURE;
    }

    if( conf->overwrite == 0 && dw_isfile(conf->outFile))
    {
        printf("%s already exist. Doing nothing\n", conf->outFile);
        bw_conf_free(&conf);
        return EXIT_SUCCESS;
    }


    if(conf->verbose > 1)
    {
        bw_conf_printf(stdout, conf);
    }
    conf->log = fopen(conf->logFile, "w");
    if(conf->log == NULL)
    {
        fprintf(stderr, "Failed to open %s for writing\n", conf->logFile);
        exit(-1);
    }
    conf->ftif = fim_tiff_new(conf->log, conf->verbose);

    fprint_time(conf->log);
    bw_conf_printf(conf->log, conf);
    fflush(conf->log);

    /* run */
    if(conf->testing)
    {
        unit_tests(conf);
        exit(0);
    }

    /* Turn off error handling in GSL. Otherwise GSL will close the program
     * if the tolerances can't be achieved. We set it back later on.
     */

    gsl_error_handler_t * old_handler =
        gsl_set_error_handler(dw_bw_gsl_err_handler); /* ignore GSL_EROUND */
    /* Do the calculations */
    BW(conf);

    gsl_set_error_handler(old_handler);
    if(GLOB_N_GSL_EROUND > 0)
    {
        printf("Warning: GSL_EROUND was raised %d times\n", GLOB_N_GSL_EROUND);
        printf("         i.e. \"roundoff error prevents tolerance from being achieved\"\n");
        fprintf(conf->log, "GSL_EROUND was raised %d times\n", GLOB_N_GSL_EROUND);
    }

    /* Write to disk  */
    if(conf->verbose > 0) {
        printf("Writing '%s' (32-bit float)\n", conf->outFile);
    }


    if(npyfilename(conf->outFile))
    {
        int shape[3] = {conf->P, conf->N, conf->M};
        npio_write(conf->outFile, 3, shape, conf->V, NPIO_F32, NPIO_F32);
    } else
    {
        ttags * T = ttags_new();
        char * swstring = malloc(1024);
        assert(swstring != NULL);
        sprintf(swstring, "deconwolf %s", deconwolf_version);
        ttags_set_software(conf->ftif, T, swstring);
        ttags_set_imagesize(conf->ftif, T, conf->M, conf->N, conf->P);
        ttags_set_pixelsize(conf->ftif, T, conf->resLateral, conf->resLateral, conf->resAxial);
        free(swstring);

        fim_tiff_write_float(conf->ftif, conf->outFile, conf->V, T, conf->M, conf->N, conf->P);
        ttags_free(&T);
    }

    fprint_time(conf->log);
    dw_gettime(&tend);
    fprintf(conf->log, "Took: %f s\n", timespec_diff(&tend, &tstart));
    fprint_peak_memory(conf->log);
    fprintf(conf->log, "done!\n");

    /* Clean up */
    free(conf->V);
    bw_conf_free(&conf);

    return EXIT_SUCCESS;
}


#ifdef STANDALONE
int main(int argc, char ** argv)
{
    return dw_bwpsf(argc, argv);
}
#endif
