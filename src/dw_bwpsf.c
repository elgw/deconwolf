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

/*
 * TODO:
 * Use something better than linear interpolation for r
*/



#include "dw_bwpsf.h"
#include "li.c"
#include "bw_int1_gsl.c"
#include "lanczos3.c"

// j0f, i.e., the float version of j0 is not available on mac
#ifdef __APPLE__
#define j0f j0
#endif

pthread_mutex_t stdout_mutex;


void dw_bw_gsl_err_handler(const char * reason,
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
    fprintf(out, "TOL = %f\n", conf->TOL);
    fprintf(out, "resLateral = %.2f nm\n", conf->resLateral);
    fprintf(out, "resAxial = %.2f nm\n", conf->resAxial);
    fprintf(out, "out size: [%d x %d x %d] pixels\n", conf->M, conf->N, conf->P);
    fprintf(out, "nThreads: %d\n", conf->nThreads);
    fprintf(out, "Verbosity: %d\n", conf->verbose);
    fprintf(out, "Overwrite: %d\n", conf->overwrite);
    fprintf(out, "File: %s\n", conf->outFile);
    fprintf(out, "Log: %s\n", conf->logFile);

    fprintf(out, "Pixel integration: %dx%dx%d samples per pixel\n",
            conf->samples_xy, conf->samples_xy, conf->samples_z);

    fprintf(out, "Oversampling of radial profile: %d X\n", conf->oversampling_R);

    if(conf->mode_int1 == MODE_INT1_LI)
    {
        fprintf(out, "BW Integral: Li's fast method enabled\n");
    }
    if(conf->mode_int1 == MODE_INT1_PG)
    {
        fprintf(out, "BW Integral: PSFGenerator clone\n");
    }
    if(conf->mode_int1 == MODE_INT1_GSL)
    {
        fprintf(out, "BW Integral: gsl_integration_qag\n");
    }
    /* How Z is handled */
    if((conf->mode_int1 == MODE_INT1_LI) | (conf->mode_int1 == MODE_INT1_PG))
    {
        fprintf(out, "Z Integral: Trapetzoid using %d samples\n", conf->samples_z);
    }
    if(conf->mode_int1 == MODE_INT1_GSL)
    {
        fprintf(out, "Z integral: gsl_integration_qag\n");
    }


    // Show theoretical FWHM of the PSF
    // And look at the number of pixels per fwhm
    double fwhm_r = 1.616340*conf->lambda/M_PI*conf->ni/conf->NA;
    double fwhm_z = 2*2.783115/M_PI*pow(conf->ni/conf->NA, 2)*conf->lambda;
    fprintf(out, "FWHM_r (lateral plane) %.2f nm\n", fwhm_r);
    double qlateral = fwhm_r/conf->resLateral;
    fprintf(out, "FWHM_r / resLateral = %.2f\n", qlateral);
    if(qlateral < 2)
    {
        fprintf(out,"! Warning: This is suboptimal, aim for at least 2\n");
    }

    fprintf(out, "FWHM_z (axial direction) %.2f nm\n", fwhm_z);
    double qaxial = fwhm_z/conf->resAxial;
    fprintf(out, "FWHM_z / resAxial = %.2f\n", qaxial);
    if(qaxial < 2)
    {
        fprintf(out,"! Warning: This is suboptimal, aim for at least 2\n");
    }
    return;
}


bw_conf * bw_conf_new()
{
    bw_conf * conf = malloc(sizeof(bw_conf));
    conf->lambda = 600;
    conf->NA = 1.45;
    conf->ni = 1.515;
    conf->TOL = 1e-1;
    conf->K = 9; // corresponding to "best" in PSFGenerator
    conf->M = 181;
    conf->N = 181;
    conf->P = 181; // Does only need to be 2*size(I,3)+1
    conf->resAxial = 300;
    conf->resLateral = 130;
    conf->nThreads = 8;
    conf->V = NULL;
    conf->verbose = 1;
    conf->V = NULL;
    conf->outFile = NULL;
    conf->logFile = NULL;
    conf->overwrite = 0;
    conf->samples_xy = 7;
    conf->samples_z = conf->samples_xy;
    /* Expose from CLI? */
    conf->oversampling_R = 15;
    conf->mode_int1 = MODE_INT1_GSL;
    conf->testing = 0;
    return conf;
}

static double timespec_diff(struct timespec* end, struct timespec * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}

void getCmdLine(int argc, char ** argv, bw_conf * s)
{
    // Copy the command line to s->cmd
    int lcmd=0;
    for(int kk = 0; kk<argc; kk++)
    {
        lcmd += strlen(argv[kk]);
    }
    lcmd += argc+2;
    s->cmd = malloc(lcmd);
    int pos = 0;
    for(int kk = 0; kk<argc; kk++)
    {
        sprintf(s->cmd+pos, "%s ", argv[kk]);
        pos += strlen(argv[kk])+1;
    }
    s->cmd[pos-1] = '\0';
}


int file_exist(char * fname)
{
    if( access( fname, F_OK ) != -1 ) {
        return 1; // File exist
    } else {
        return 0;
    }
}

void fprint_time(FILE * f)
{
    f == NULL ? f = stdout : 0;
    time_t now = time(NULL);
    char * tstring = ctime(&now);
    fprintf(f, "%s\n", tstring);
}

void bw_version(FILE* f)
{
    fprintf(f, "deconwolf: '%s' PID: %d\n", deconwolf_version, (int) getpid());
}

void usage(__attribute__((unused)) int argc, char ** argv, bw_conf * s)
{
    printf("Usage:\n");
    printf("\t%s <options> psf.tif \n", argv[0]);
    printf("\n");
    printf("Required parameters:\n");
    printf(" --NA\n\t Set numerical aperture.\n");
    printf(" --lambda\n\t Set emission wavelength [nm].\n");
    printf(" --ni\n\t Set refractive index.\n");
    printf(" --resxy\n\t Set lateral pixel size (x-y) [nm].\n");
    printf(" --resz\n\t Set the axial pixel size in (z) [nm].\n");
    printf("Optional:\n");
    printf(" --size N \n\t Set output size to N x N x P [pixels].\n"
        "\t This should be set big enough so that the PSF cone\n"
        "\t isn't cropped at the first and last slice of the PSF.\n");
    printf(" --nslice P\n\t Set output size to N x N x P [pixels].\n");
    printf("\t P has to be an odd number. Please not that P should be at\n"
           "\t least (2Z-1) where Z is the number of slices in the image\n"
           "\t to be deconvolved\n");
    printf(" --threads\n\t Set number of threads. Currently set to %d\n", s->nThreads);
    printf("Optional, related to quality\n");
    printf(" --samples Q\n\t Sets the samples to QxQxQ for pixel integration. "
           "Default: %d\n", s->samples_xy);
    printf("\t Q has to be an odd number.\n");
    printf("\t Unless --li or --pg, integration over Z is done to machine precision\n");
    printf(" --li\n\t Use Li's method for the BW integral.\n");
    printf("\t About 14x faster for the default PSF size.\n"
           "\t Not thoroughly tested.\n");
    printf("\t Possibly unstable for large large PSFs.");
    if(s->mode_int1 == MODE_INT1_LI)
    {
        printf(" Default.");
    }
    printf("\n");
    printf(" --gsl\n\t Use GSL for the BW integral as well as the integration over z.");
    if(s->mode_int1 == MODE_INT1_GSL)
    {
        printf(" Default.");
    }
    printf("\n");
    printf(" --pg\n"
           "\t Use something similar to the PSFGenerator for the BW integral. Will be removed in future versions");
    if(s->mode_int1 == MODE_INT1_PG)
    {
        printf(" Default.");
    }
    printf("\n");

    printf("Other options:\n");
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

void bw_argparsing(int argc, char ** argv, bw_conf * s)
{
    bw_conf * defaults = malloc(sizeof(bw_conf));
    memcpy(defaults, s, sizeof(bw_conf));

    if(argc < 2)
    {
        // Well it could run with the defaults but that does not make sense
        usage(argc, argv, defaults);
        exit(-1);
    }

    getCmdLine(argc, argv, s);

    struct option longopts[] =
        {
         { "version",      no_argument,       NULL, 'v'},
         { "help",         no_argument,       NULL, 'h'},
         // Settings
         { "lambda",       required_argument, NULL, 'l'},
         { "threads",      required_argument, NULL, 't'},
         { "verbose",      required_argument, NULL, 'p'},
         { "overwrite",    no_argument,       NULL, 'w'},
         { "NA",           required_argument, NULL, 'n'},
         { "ni",           required_argument, NULL, 'i'},
         { "resxy",        required_argument, NULL, 'x'},
         { "resz",         required_argument, NULL, 'z'},
         { "size",         required_argument, NULL, 'N'},
         { "nslice",       required_argument, NULL, 'P'},
         { "samples",      required_argument, NULL, 'q'},
         { "li",         no_argument,       NULL, 'f'},
         { "gsl",         no_argument,       NULL, 'g'},
         { "pg",           no_argument,       NULL, 'G'},
         { "testing",       no_argument, NULL, 'T'},
         { NULL,           0,                 NULL,  0 }
        };
    int Pset = 0;
    int ch;
    while((ch = getopt_long(argc, argv, "fqGvho:t:p:w:n:i:x:z:N:P:l:T", longopts, NULL)) != -1)
    {
        switch(ch) {
        case 'T':
            s->testing = 1;
            break;
        case 'f':
            s->mode_int1 = MODE_INT1_LI;
            break;
        case 'g':
            s->mode_int1 = MODE_INT1_GSL;
            break;
        case 'G':
            s->mode_int1 = MODE_INT1_PG;
            break;
        case 'q':
            s->samples_xy = atoi(optarg);
            s->samples_z = atoi(optarg);
            break;
        case 'v':
            bw_version(stdout);
            exit(0);
            break;
        case 'h':
            usage(argc, argv, defaults);
            exit(0);
        case 'o':
            s->outFile = malloc(strlen(optarg) + 1);
            sprintf(s->outFile, "%s", optarg);
            break;
        case 't':
            s->nThreads = atoi(optarg);
            break;
        case 'p':
            s->verbose = atoi(optarg);
            break;
        case 'w':
            s->overwrite = 1;
            break;
        case 'n':
            s->NA = atof(optarg);
            break;
        case 'i':
            s->ni = atof(optarg);
            break;
        case 'x':
            s->resLateral = atof(optarg);
            break;
        case 'z':
            s->resAxial = atof(optarg);
            break;
        case 'l':
            s->lambda = atof(optarg);
            break;
        case 'N':
            s->M = atoi(optarg);
            s->N = s->M;
            break;
        case 'P':
            s->P = atoi(optarg);
            Pset = 1;
            break;
        default:
            exit(EXIT_FAILURE);
            break;
        }
    }

    free(defaults);

    if(s->lambda < 50)
    {
        fprintf(stderr, "Error: lambda has be be at least 50 nm\n");
        exit(-1);
    }

    if(s->lambda > 10000)
    {
        fprintf(stderr, "Error lambda can be at most 1000 nm\n");
        exit(-1);
    }

    if(s->samples_xy % 2 == 0)
    {
        fprintf(stderr, "The number of integration points in x-y (per dimension) has to be odd\n");
        exit(-1);
    }

    if(s->M % 2 == 0)
    {
        printf("Warning: The size has to be odd, 1, 3, ...\n");
        printf("Changing from %d to %d.\n", s->M, s->M+1);
        s->M++;
        s->N++;
    }

    if(s->P % 2 == 0)
    {
        printf("Warning: The number of slices has to be odd, 1, 3, ...\n");
        printf("          Changing from %d to %d.\n", s->P, s->P+1);
        s->P++;
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
        sprintf(s->outFile, "%s", argv[argc-1]);
    }

    if(s->outFile == NULL)
    {
        s->outFile = malloc(100*sizeof(char));
        sprintf(s->outFile, "PSFBW_%.2f_%.2f_%.1f_%.1f_%.1f.tif",
                s->NA, s->ni, s->lambda, s->resLateral, s->resAxial);
    }


    s->logFile = malloc(strlen(s->outFile) + 10);
    sprintf(s->logFile, "%s.log.txt", s->outFile);
    return;
}

double cabs2(double complex v)
{
    return pow(creal(v), 2) + pow(cimag(v), 2);
}

void * BW_thread(void * data)
{
    // Entry point for pthread_create
    bw_conf * conf = (bw_conf *) data;

    if(conf->verbose > 3)
    {
        pthread_mutex_lock(&stdout_mutex);
        printf("-> From thread %d\n", conf->thread);
        bw_conf_printf(stdout, conf);
        pthread_mutex_unlock(&stdout_mutex);
    }

    float * V = conf->V;
    int M = conf->M;
    int N = conf->N;
    int P = conf->P;
    size_t MN = M*N;

    for (int z = conf->thread; z <= (P-1)/2; z+=conf->nThreads) {
        float defocus = conf->resAxial * (z - (P - 1.0) / 2.0);
        BW_slice(V + z*MN, defocus, conf);
    }
    return NULL;
}

void BW(bw_conf * conf)
{
    assert(conf->V == NULL);
    assert(conf->nThreads > 0);
    conf->V = malloc(conf->M*conf->N*conf->P*sizeof(float));

    float * V = conf->V;
    int M = conf->M;
    int N = conf->N;
    int P = conf->P;
    size_t MN = M*N;

    int nThreads = conf->nThreads;

    if(nThreads == 1)
    {
        for (int z = 0; z <= (P-1)/2; z++) {
            float defocus = conf->resAxial * (z - (P - 1.0) / 2.0);
            printf("P = %d, z = %d, defocus = %f\n", P, z, defocus);
            BW_slice(conf->V + z*M*N, defocus, conf);
        }
    } else {
        pthread_t * threads = malloc(nThreads*sizeof(pthread_t));
        bw_conf ** confs = malloc(nThreads*sizeof(bw_conf*));

        for(int kk = 0; kk<nThreads; kk++)
        {
            confs[kk] = (bw_conf*) malloc(sizeof(bw_conf));
            memcpy(confs[kk], conf, sizeof(bw_conf));
            confs[kk]->thread = kk;
            pthread_create(&threads[kk], NULL, BW_thread, (void *) confs[kk]);
        }

        for(int kk = 0; kk<nThreads; kk++)
        {
            pthread_join(threads[kk], NULL);
            free(confs[kk]);
        }
        free(confs);
        free(threads);
    }

    // symmetry in Z
    for(int z = 0; z<(P-1)/2; z++)
    {
        memcpy(V+(P-z-1)*MN, V+z*MN, MN*sizeof(float));
    }
}


float complex integrand(float rho, float r, float z, bw_conf * conf) {
     assert(rho<=1); assert(rho>=0);

    float k0 = 2.0 * M_PI / conf->lambda;
    // j0 has smaller errors but I don't think that we need that precision.
    // j0f is much faster

#if 0
    // This is how it looks in the Java code
    float BesselValue = j0f(k0 * conf->NA * r * rho);
    // Optical path difference
    float OPD = pow(conf->NA,2) * z * pow(rho,2) / (2.0 * conf->ni);
#endif

    // This is how it looks on their web page
    // as well as in Optics by Born and Wolf,
    // p. 437 (section 8.8), 6th edition.
    // "Three-dimensional light distribution near focus"
    float q = conf->NA / conf->ni;
    float BesselValue = j0f(k0 * q * r * rho);
    float OPD = pow(q,2) * z * pow(rho,2) / (2.0);

    // Phase aberrations.
    float W = k0 * OPD;
    //printf("k0=%e, OPD=%e, W=%e\n", k0, OPD, W);
    float complex y = BesselValue*cos(W)*rho - I*BesselValue*sin(W)*rho;
    return y;
}


double complex calculate(float r, float defocus, bw_conf * conf) {
    // Simpson approximation for the Kirchhoff diffraction integral

    float curDifference = conf->TOL; // Stopping criterion

    complex float value = 0;
    float curI = 0.0;
    float prevI = 0.0;

    // Initialization of the Simpson sum (first iteration)
    int64_t N = 2; // number of sub-intervals
    double del = 1.0/N; // sub-interval length

    // Initial 3-point estimation
    //printf("r = %f, defocus = %f\n", r, defocus);
    complex float sumOddIndex = integrand(0.5, r, defocus, conf);
    //printf("Initial: %f, %fi\n", creal(sumOddIndex), cimag(sumOddIndex) );
    complex float sumEvenIndex = 0;
    complex float valueX0 = integrand(0.0, r, defocus, conf);
    float complex  valueXn = integrand(1.0, r, defocus, conf);
    double complex sum = valueX0 + 2*sumEvenIndex + 4*sumOddIndex + valueXn;
    curI = cabs2(sum) * pow(del,2);
    prevI = curI;

    // (almost) double the number of points until desired precision
    size_t k = 0; // number of consecutive successful approximations
    size_t iteration = 1;
    size_t max_iterations = 15;

    while (k < conf->K && iteration <= max_iterations) {
        iteration++;
        N *= 2;
        del /= 2;

        sumEvenIndex = sumEvenIndex + sumOddIndex;
        sumOddIndex = 0 + I*0;

        for (int64_t n = 1; n < N; n = n + 2) {
            double rho = n * del;
            value = integrand(rho, r, defocus, conf);
            sumOddIndex += value;
        }

        sum = valueX0 + 2*sumEvenIndex + 4*sumOddIndex + valueXn;
        //printf("N=%ld, sum=%f %fi\n", N, creal(sum), cimag(sum));

        curI = cabs2(sum) * pow(del, 2);

        // Relative error between consecutive approximations
        if (prevI == 0.0)
        {
            curDifference = fabs((prevI - curI) / 1E-5);
        } else {
            curDifference = fabs((prevI - curI) / curI);
        }
        if (curDifference <= conf->TOL)
        {
            k++;
        } else {
            k = 0;
        }

        // This gives a quick finish for close-to-zero pixels
        // the relative error is not that interesting for those
        if(fabs(curI-prevI) < 1e-12)
        { break; }

        prevI = curI;
    }

    //printf("Final: N=%ld, sum=%f %fi, ZW=%ld\n", N, creal(sum), cimag(sum), (N-1)*2+2);
    return sum/((double) (N)*3.0);
}


double simp1_weight(size_t k, size_t N)
{
    /* Weights for 1D-Simpson
     * [1, 4, 1] N=3
     * [1, 4, 2, 4, 1] N=5
     * [1, 4, 2, 4, 2, 4, 1] N=7, */
    assert(N%2 == 1);
    if(k == 0 || k == N-1)
        return 1.0;
    if(k%2 == 0)
        return 2.0;
    return 4.0;
}

/* OVER_SAMPLING: Number of radial samples per dx
 *  xc, yc
 *
 * x0, x1, y0, y1 : integration region given in _pixels_
 * N: number of integration points per dimension
 */

double simp2(int OVER_SAMPLING,
             double xc, double yc,
             double * H, double * R,
             double x0, double x1,
             double y0, double y1, int N)
{

    double dx = (x1-x0)/(double) (N-1);
    double dy = (y1-y0)/(double) (N-1);
    double v = 0;
    size_t sW = 0; /* Sum of weights */
    for(int kk = 0; kk<N; kk++)
    {
        double wx = simp1_weight(kk, N);
        for(int ll = 0; ll<N; ll++)
        {
            double wy = simp1_weight(ll, N);
            double w = wx*wy;
            double x = x0+(double) kk*dx;
            double y = y0+(double) ll*dy;
            double r = sqrt(pow(x-xc, 2) + pow(y-yc, 2));
            int index = (int) floor(r*OVER_SAMPLING);
            double blend = (r-R[index])*OVER_SAMPLING;

            sW+=w;
            //double ival = (H[index] + (H[index+1] - H[index]) * blend);
            // printf("w=%f: x=%f, y=%f r=%f, ival=%f, index = %d, blend = %f\n", w, x, y, r, ival, index, blend);

            double ival_l3 = lanczos3(H, 10000, r*(double) OVER_SAMPLING);
            //printf("Linear/Lanczos3: %f %f %e\n", ival, ival_l3, fabs(ival-ival_l3));

            //v+=w*ival;
            v+=w*ival_l3;
        }
    }
    v = v/sW;

    return v;
}


void BW_slice(float * V, float z, bw_conf * conf)
{
    // Calculate the PSF for a single plane/slice
    //printf("z = %e\n", z);

    // V is a pointer to the _slice_ not the whole volume
    // The center of the image in units of [pixels]
    double x0 = (conf->M - 1) / 2.0;
    double y0 = (conf->N - 1) / 2.0;

    //int maxRadius = (int) round(sqrt(pow(conf->M - x0, 2) + pow(conf->N - y0, 2))) + 1;
    int maxRadius = (int) ceil(fabs(conf->M - x0))+1;
    int OVER_SAMPLING = conf->oversampling_R;

    size_t nr = maxRadius*conf->oversampling_R;
    double * r = malloc(nr * sizeof(double));

    // abs(x)^2
    double * hReal = malloc(nr * sizeof(double));

    /* Calculate the radial profile
       possibly by integrating over Z */

    for(size_t n = 0; n < nr; n++)
    {
        r[n] = ((double) n) / ((double) OVER_SAMPLING);
    }


    /* Calculate 1D integral using something similar to the PSFGenerator */
    if(conf->mode_int1 == MODE_INT1_PG)
    {
        int NS = conf->samples_z;
        for (size_t n = 0; n < nr; n++) {
            hReal[n] = 0;

            double complex now_complex = 0;
            double now_real = 0;

            double dz = (conf->resAxial)/(NS-1);
            double sW = 0;
            for(int kk = 0; kk < NS; kk++)
            {
                double zs = z + kk*dz -(NS-1)/2*dz;
                //printf("zs=%f ", zs);
                double w = simp1_weight(kk, NS);
                sW += w;
                //                printf("r[%zu] = %f, resLateral=%f, zs=%f", n, r[n], conf->resLateral, zs);
                double complex value = calculate(r[n] * conf->resLateral, zs, conf);
                now_complex += w*value;
                now_real += w*cabs2(value);

                //printf("(%f %fi)\n", creal(value), cimag(value));
            }
            //printf("sW = %f\n", sW);

            hReal[n] = now_real/sW;
        }
    }

    /* 1D integral using Bessel series, li */
    if(conf->mode_int1 == MODE_INT1_LI)
    {
        int NS = conf->samples_z;
        // NS = 1; Just for testing
        memset(hReal, 0, nr*sizeof(double));


        double dz = 0;
        if(NS > 1) {
            dz = (conf->resAxial)/(NS-1);
        }

        for(int kk = 0; kk < NS; kk++) // Over z
        {
            double zs = z + kk*dz - (NS-1)/2*dz;
            double w = simp1_weight(kk, NS);
            //printf("z = %e, dz = %e, zs = %e\n", z, dz, zs); fflush(stdout);
            li_conf * L = li_new(zs);
            L->lambda = conf->lambda;
            L->NA = conf->NA;
            L->ni = conf->ni;
            for(size_t n = 0; n < nr; n++) // over r
            {
                complex double v = li_calc(L, r[n]*conf->resLateral);
                hReal[n] += w*cabs2(v);
            }
            li_free(&L);
        }
        double wtot = 0;
        for(int kk = 0; kk < NS; kk++)
            wtot += simp1_weight(kk, NS);

        for(size_t n = 0; n < nr; n++)
        {
            hReal[n]/=wtot;
        }
    }


    /* BW integral and Z integral by GSL */
    if(conf->mode_int1 == MODE_INT1_GSL)
    {
        struct my_params gsl_params;
        gsl_params.ncalls = 0;
        gsl_params.NA = conf->NA;
        gsl_params.ni = conf->ni;

        gsl_params.lambda = conf->lambda;
        gsl_params.ninterval = 10000;
        gsl_integration_workspace * gsl_workspace1 =
            gsl_integration_workspace_alloc(gsl_params.ninterval);
        gsl_integration_workspace * gsl_workspaceZ =
            gsl_integration_workspace_alloc(gsl_params.ninterval);

        /* Integration over Z for a equidistant set of radii*/
        gsl_params.z0 = z - conf->resAxial / 2.0;
        gsl_params.z1 = z + conf->resAxial / 2.0;
        //printf("Integration over z: [%f, %f]\n", gsl_params.z0, gsl_params.z1);

        for(size_t n = 0; n < nr; n++) // over r
        {
            gsl_params.r = r[n]*conf->resLateral;
            double value = integrate_bw_gsl_z(&gsl_params,
                                              gsl_workspace1, gsl_workspaceZ);
            hReal[n] = value;
            //printf("r = %f -> %f\n", gsl_params.r, value); /* Looks ok */
        }

        gsl_integration_workspace_free(gsl_workspace1);
        gsl_integration_workspace_free(gsl_workspaceZ);
    }

    assert(r[0] == 0);

    /*
     * For each pixel in the quadrant integrate over x and y
     */

    for (int x = 0; 2*x <= conf->M; x++) {
        for (int y = 0; 2*y <= conf->N; y++) {

            // radius of the current pixel in units of [pixels]
            double rPixel = sqrt( pow((double) x - x0, 2)
                                  + pow((double) y - y0, 2));
            double pIntensity = 0;

            /* Pixel integration */
            assert(conf->samples_xy > 0);

            if(rPixel <= maxRadius)
            {
                pIntensity = simp2(OVER_SAMPLING,
                                   x0, y0, hReal, r,
                                   x - 0.5, x + 0.5,
                                   y - 0.5, y + 0.5, conf->samples_xy);
            } else {
                pIntensity = 0;
            }


            // Use symmetry to fill the output image plane
            int xf = conf->M-x-1;
            int yf = conf->N-y-1;
            V[x + conf->M * y] = pIntensity;
            V[y + conf->M * x] = pIntensity;
            V[xf + conf->M * y] = pIntensity;
            //      V[yf + conf->M * x] = v;
            V[x + conf->M * yf] = pIntensity;
            //      V[y + conf->M * xf] = v;
            V[xf + conf->M * yf] = pIntensity;
            V[yf + conf->M * xf] = pIntensity;
        }
    }

    free(r);
    free(hReal);
    return;
}

void unit_tests(bw_conf * conf)
{
    printf("Values of PG-integrand for different values of rho\n");
    double r = 0; double z = 0;
    printf("r = %f, z = %f, NA=%f, ni=%f, lambda=%f nm\n", r, z, conf->NA, conf->ni, conf->lambda);
    for(double rho = 0; rho<=1; rho+=0.1)
    {
        double complex value = integrand(rho, r, z, conf);
        printf("rho=%f, %f %fi\n", rho, creal(value), cimag(value));
    }

    r = 130; z = 130;
    printf("r = %f, z = %f, NA=%f, ni=%f, lambda=%f nm\n", r, z, conf->NA, conf->ni, conf->lambda);
    for(double rho = 0; rho<=1; rho+=0.1)
    {
        double complex value = integrand(rho, r, z, conf);
        printf("rho=%f, %f %fi\n", rho, creal(value), cimag(value));
    }

    z = 0;
    li_conf * L = li_new(z);
    L->lambda = conf->lambda;
    L->NA = conf->NA;
    L->ni = conf->ni;

    li_show(L);
    assert(L != NULL);

    for(double r = 0; r<500; r+=130)
    {
        double complex v = li_calc(L, r);
        double vreal = pow(creal(v),2) + pow(cimag(v),2);
        printf("PSF(r=%f, z=%f) = (%f %fi) %f dV\n", r, z, creal(v), cimag(v), vreal);
        double complex v2 = calculate(r, z, conf);
        double vreal2 = pow(creal(v2),2) + pow(cimag(v2),2);
        printf("PG(r=%f, z=%f) = (%f %fi) %f dV\n", r, z, creal(v2), cimag(v2), vreal2);
    }

    L = li_free(&L);
    assert(L == NULL);

}

int main(int argc, char ** argv)
{
    struct timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);

    // Use defaults
    bw_conf * conf = bw_conf_new();
    bw_argparsing(argc, argv, conf);

    if( conf->overwrite == 0 && file_exist(conf->outFile))
    {
        printf("%s already exist. Doing nothing\n", conf->outFile);
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
    fprint_time(conf->log);
    bw_conf_printf(conf->log, conf);

    // Run
    if(conf->testing)
    {
        unit_tests(conf);
        exit(0);
    }

    /* Turn off error handling in GSL. Otherwise GSL will close the program
     * if the tolerances can't be achieved. We set it back later on.
     * Ideally we should have a custom error handler that simply ignores
     * GSL_EROUND */
    gsl_error_handler_t * old_handler = gsl_set_error_handler(dw_bw_gsl_err_handler);

    /* Do the calculations */
    BW(conf);

    gsl_set_error_handler(old_handler);

    // Write to disk
    if(conf->verbose > 0) {
        printf("Writing as 32-bit floats to %s\n", conf->outFile);

    }

    ttags * T = ttags_new();
    char * swstring = malloc(1024);
    sprintf(swstring, "deconwolf %s", deconwolf_version);
    ttags_set_software(T, swstring);
    ttags_set_imagesize(T, conf->M, conf->N, conf->P);
    ttags_set_pixelsize(T, conf->resLateral, conf->resLateral, conf->resAxial);
    free(swstring);

    fim_tiff_write_float(conf->outFile, conf->V, T, conf->M, conf->N, conf->P, conf->log);
    ttags_free(&T);

    fprint_time(conf->log);
    clock_gettime(CLOCK_REALTIME, &tend);
    fprintf(conf->log, "Took: %f s\n", timespec_diff(&tend, &tstart));
    fprintf(conf->log, "done!\n");

    // Clean up
    free(conf->V);
    fclose(conf->log);
    free(conf->outFile);
    free(conf->logFile);
    free(conf->cmd);
    free(conf);

    return EXIT_SUCCESS;
}
