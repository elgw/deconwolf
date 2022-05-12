/* (C) 2022 Erik L. G. Wernersson
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

#include "dw_psf.h"

typedef struct{
    double NA;
    double ni;
    double dx;
    double dz;
    double lambda;
    double lambda2;
    double pinhole;
} optical_t;

typedef struct{
    int overwrite;
    int verbose;
    char * outfile; /* Where to write the tsv output */
    char * logfile;
    FILE * log;
    int nthreads;
    int M;
    int P;
    optical_t optical;
    /* How many times larger the calculation grid should be compared to images dimensions */
    int xgrid;
    /* Process at x times higher resolution and then average */
    int oversampling;
} opts;

static opts * opts_new();
static void opts_free(opts * s);
static void usage(__attribute__((unused)) int argc, char ** argv);
static void argparsing(int argc, char ** argv, opts * s);
static int file_exist(char * fname);

static opts * opts_new()
{
    opts * s = malloc(sizeof(opts));

    s->overwrite = 0;
    s->verbose = 1;
    s->outfile = NULL;
    s->logfile = NULL;
    s->log = NULL;
    s->nthreads = dw_get_threads();
    s->optical.NA = 1.45;
    s->optical.ni = 1.515;
    s->optical.dx = 65;
    s->optical.dz = 200;
    s->optical.lambda = 460;
    s->optical.lambda2 = 0;
    s->optical.pinhole = 1; /* Airy Units, 1.22*lambda/NA */
    s->M = 181;
    s->P = 181;
    s->xgrid = 7;
    s->oversampling = 1;
    return s;
}


static void opts_free(opts * s)
{
    dw_nullfree(s->outfile);
    dw_nullfree(s->logfile);
    if(s->log != NULL)
    {
        fclose(s->log);
    }
    free(s);
}

static void opts_print(FILE * f, opts * s)
{
    fprintf(f, "overwrite = %d\n", s->overwrite);
    fprintf(f, "verbose = %d\n", s->verbose);
    fprintf(f, "nthreads = %d\n", s->nthreads);
    fprintf(f, "outfile: %s\n", s->outfile);
    fprintf(f, "Out size: [%d x %d x %d] pixels\n", s->M, s->M, s->P);
    fprintf(f, "NA=%f\n", s->optical.NA);
    fprintf(f, "ni=%f\n", s->optical.ni);
    fprintf(f, "dx=%f nm\n", s->optical.dx);
    fprintf(f, "dz=%f nm\n", s->optical.dz);
    fprintf(f, "Emission lambda=%f nm\n", s->optical.lambda);
    if(s->optical.lambda2 > 0)
    {
        fprintf(f, "Excitation lambda: %f nm\n", s->optical.lambda2);
        fprintf(f, "Pinhole size %f AU\n", s->optical.pinhole);
        fprintf(f, "Pinhole shape: square\n");
    }
    return;
}

static void usage(__attribute__((unused)) int argc, char ** argv)
{
    opts * s = opts_new();
    printf("usage: %s [<options>] output.tif\n", argv[0]);
    printf("Options:\n");
    printf(" --NA NA\n\t Set numerical aperture\n");
    printf(" --ni ni\n\t Set refractive index\n");
    printf(" --dx dx\n\t Lateral pixel size\n");
    printf(" --dz dz\n\t Axial pixel size\n");
    printf(" --lambda l\n\t Emission wave length\n");
    printf(" --lambda2 l\n\t Excitation wave length\n");
    printf(" --pinhole p\n\t Pinhole size in AU\n");
    printf(" --size s\n\t Set the lateral size (pixels). Default: %d\n", s->M);
    printf(" --nslice p\n\t Set the number of planes (pixels). Default: %d\n", s->P);
    printf(" --xgrid x\n\t Set the factor for lateral image padding. Default: %d\n", s->xgrid);
    printf(" --oversample x\n\t Process the PSF at x times higher resolution. Default: %d\n", s->oversampling);
    printf("General:\n");
    printf(" --overwrite\n\t Overwrite existing files\n");
    printf(" --help\n\t Show this message\n");
    printf(" --verbose v\n\t Verbosity level\n");
    printf("\n");
    opts_free(s);
    return;
}

static void argparsing(int argc, char ** argv, opts * s)
{
    if(argc == 1)
    {
        printf("Deconwolf %s PSF generator.\n", deconwolf_version);
        printf("See `dw psf --help` or `man dw`.\n");
        exit(EXIT_SUCCESS);
    }
    int hasNA = 0;
    int hasni = 0;
    int hasdx = 0;
    int hasdz = 0;
    int haslambda = 0;

    struct option longopts[] = {
        {"oversample", required_argument, NULL, 'O'},
        {"xgrid", required_argument, NULL, 'g'},
        {"lambda", required_argument, NULL, 'L'},
        {"lambda2", required_argument, NULL, '2'},
        {"NA",     required_argument, NULL, 'n'},
        {"dx",     required_argument, NULL, 'x'},
        {"dz",     required_argument, NULL, 'z'},
        {"ni",     required_argument, NULL, 'i'},
        {"help", no_argument, NULL, 'h'},
        {"overwrite", no_argument, NULL, 'o'},
        {"nslice", required_argument, NULL, 'p'},
        {"pinhole", required_argument, NULL, 'P'},
        {"size", required_argument, NULL, 's'},
        {"threads", required_argument, NULL, 't'},
        {"verbose", required_argument, NULL, 'v'},
        {NULL, 0, NULL, 0}};
    int ch;
    while((ch = getopt_long(argc, argv,
                            "O:g:L:2:n:x:z:i:hop:s:v:P:", longopts, NULL)) != -1)
    {
        switch(ch){
        case 'O':
            s->oversampling = atoi(optarg);
            break;
        case 'g':
            s->xgrid = atoi(optarg);
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
            break;
        case 'i':
            s->optical.ni = atof(optarg);
            hasni = 1;
            break;
        case 'L':
            s->optical.lambda = atof(optarg);
            haslambda = 1;
            break;
        case '2':
            s->optical.lambda2 = atof(optarg);
            break;
        case 'o':
            s->overwrite = 1;
            break;
        case 't':
            s->nthreads = atoi(optarg);
            break;
        case 'v':
            s->verbose = atoi(optarg);
            break;
        case 'x':
            s->optical.dx = atof(optarg);
            hasdx = 1;
            break;
        case 'z':
            s->optical.dz = atof(optarg);
            hasdz = 1;
            break;
        case 'n':
            s->optical.NA = atof(optarg);
            hasNA = 1;
            break;
        case 'p':
            s->P = atoi(optarg);
            break;
        case 'P':
            s->optical.pinhole = atof(optarg);
            break;
        case 's':
            s->M = atoi(optarg);
            break;
        }
    }
    int ok = 1;
    if(haslambda == 0)
    {
        fprintf(stderr, "Emission wave length not specified (--lambda)\n");
        ok = 0;
    }
    if(hasni == 0)
    {
        fprintf(stderr, "Refractive index not specified (--ni)\n");
        ok = 0;
    }
    if(hasNA == 0)
    {
        fprintf(stderr, "Numerical Aperture not specified (--NA)\n");
        ok = 0;
    }

    if(hasdx == 0)
    {
        fprintf(stderr, "Lateral pixel size not specified (--dx)\n");
        ok = 0;
    }

    if(hasdz == 0)
    {
        fprintf(stderr, "Axial pixel size not specified (--dz)\n");
        ok = 0;
    }

    if(ok == 0)
    {
        exit(EXIT_FAILURE);
    }

    if(s->xgrid < 1)
    {
        printf("%d is an invalid value for --xgrid, please use a value >= 1\n", s->xgrid);
        exit(EXIT_FAILURE);
    }
    if(s->xgrid % 2 == 0)
    {
        s->xgrid++;
    }
    /* We work with odd number of pixels */
    if(s->M % 2 == 0)
    {
        s->M++;
    }

    if(s->P % 2 == 0)
    {
        s->P++;
    }

    if(optind == argc)
    {
        s->outfile = malloc(1024);
        if(s->optical.lambda2 == 0)
        {
        sprintf(s->outfile,
                "psf_NA%.2f_ni%.2f_lambda%.0f_dx%.0f_dz%.0f.tif",
                s->optical.NA, s->optical.ni, s->optical.lambda,
                s->optical.dx, s->optical.dz);
        } else {
            sprintf(s->outfile,
                    "cpsf_NA%.2f_ni%.2f_lambda%.0f_lambda2%.0f_dx%.0f_dz%.0f.tif",
                    s->optical.NA, s->optical.ni,
                    s->optical.lambda, s->optical.lambda2,
                    s->optical.dx, s->optical.dz);
        }
        if(s->verbose > 0)
        {
            printf("Will write to %s\n", s->outfile);
        }
    } else {
        s->outfile = strdup(argv[optind]);
    }

    if(s->overwrite == 0)
    {
        if(file_exist(s->outfile))
        {
            printf("%s exists, leaving\n", s->outfile);
            printf("if you don't care, use --overwrite\n");
            exit(EXIT_FAILURE); /* TODO */
        }
    }

    s->logfile = malloc(strlen(s->outfile) + 32);
    sprintf(s->logfile, "%s.log.txt", s->outfile);
    s->log = fopen(s->logfile, "w");
    if(s->log == NULL)
    {
        fprintf(stderr, "Failed to open %s for writing\n", s->logfile);
        exit(EXIT_FAILURE);
    }

    for(int kk = 0; kk< argc; kk++)
    {
        fprintf(s->log, "%s", argv[kk]);
        if(kk+1 != argc)
        {
            fprintf(s->log, " ");
        } else {
            fprintf(s->log, "\n");
        }
    }

    opts_print(s->log, s);
    if(s->verbose > 1)
    {
        opts_print(stdout, s);
    }
    return;
}


static int file_exist(char * fname)
{
    if( access( fname, F_OK ) != -1 ) {
        return 1; // File exist
    } else {
        return 0;
    }
}

void shift_double(double * V, int N, int stride, int shift)
{
    double * buff = malloc(N*sizeof(double));
    for(int kk = 0; kk<N; kk++)
    {
        buff[kk] = V[kk*stride];
    }
    shift = shift - N;
    for(int kk = 0; kk<N; kk++)
    {
        int to = kk;
        int from = ((kk - shift) % N);
        //printf("%d -> %d\n", from, to); fflush(stdout);
        assert(to >= 0);
        assert(to < N);
        assert(from >= 0);
        assert(from < N);
        V[to*stride] = buff[from];
    }

    free(buff);
}

void shift_float(float * V, int N, int stride, int shift)
{
    float * buff = malloc(N*sizeof(float));
    for(int kk = 0; kk<N; kk++)
    {
        buff[kk] = V[kk*stride];
    }
    shift = shift - N;
    for(int kk = 0; kk<N; kk++)
    {
        int to = kk;
        int from = ((kk - shift) % N);
        //printf("%d -> %d\n", from, to); fflush(stdout);
        assert(to >= 0);
        assert(to < N);
        assert(from >= 0);
        assert(from < N);
        V[to*stride] = buff[from];
    }

    free(buff);
}


void shift_complex(fftw_complex * V, int N, int stride, int shift)
{
    fftw_complex * buff = malloc(N*sizeof(fftw_complex));
    for(int kk = 0; kk<N; kk++)
    {
        buff[kk][0] = V[kk*stride][0];
        buff[kk][1] = V[kk*stride][1];
    }
    shift = shift - N;
    for(int kk = 0; kk<N; kk++)
    {
        int to = kk;
        int from = ((kk - shift) % N);
        //printf("%d -> %d\n", from, to); fflush(stdout);
        V[to*stride][0] = buff[from][0];
        V[to*stride][1] = buff[from][1];
    }

    free(buff);
}


void ifftshift2_double(double * ph, int M)
{
    /* Assuming ph is MxM */
    assert(M%2 == 1);
    /* Shift in dimension 1 */
    for(int kk = 0; kk<M; kk++)
    {
        shift_double(ph+kk*M, M, 1, (M-1)/2);
    }
    /* Shift in dimension 2 */
    for(int kk = 0; kk<M; kk++)
    {
        shift_double(ph+kk, M, M, (M-1)/2);
    }
}



void ifftshift2_float(float * ph, int M)
{
    /* Assuming ph is MxM */
    assert(M%2 == 1);
    /* Shift in dimension 1 */
    for(int kk = 0; kk<M; kk++)
    {
        shift_float(ph+kk*M, M, 1, (M-1)/2);
    }
    /* Shift in dimension 2 */
    for(int kk = 0; kk<M; kk++)
    {
        shift_float(ph+kk, M, M, (M-1)/2);
    }
}


void ifftshift2_complex(fftw_complex * ph, int M)
{
    assert(M%2 == 1);

    /* Shift in dimension 1 */
    for(int kk = 0; kk<M; kk++)
    {
        shift_complex(ph+kk*M, M, 1, (M-1)/2);
    }
    /* Shift in dimension 2 */
    for(int kk = 0; kk<M; kk++)
    {
        shift_complex(ph+kk, M, M, (M-1)/2);
    }

}

void fftshift2_complex(fftw_complex * ph, int M)
{
    /* Shift in dimension 1 */
    for(int kk = 0; kk<M; kk++)
    {
        shift_complex(ph+kk*M, M, 1, (M-1)/2);
    }
    /* Shift in dimension 2 */
    for(int kk = 0; kk<M; kk++)
    {
        shift_complex(ph+kk, M, M, (M-1)/2);
    }

}

void fftshift2_float(float * ph, int M)
{
    /* Assuming ph is MxM */
    assert(M%2 == 1);
    /* Shift in dimension 1 */
    for(int kk = 0; kk<M; kk++)
    {
        shift_float(ph+kk*M, M, 1, (M+1)/2);
    }
    /* Shift in dimension 2 */
    for(int kk = 0; kk<M; kk++)
    {
        shift_float(ph+kk, M, M, (M+1)/2);
    }
}



static void downsample_integrate(float * A, float * B, int nA, int nB)
{
    int factor = nB/nA;
    int mf = (factor-1)/2; /* Offset to hit the middle of the factor x factor regions */
    //printf("Integrating over %dx%d pixels\n", factor, factor);
    /* Reset A */
    memset(A, 0, nA*nA*sizeof(float));

    float * K = malloc(sizeof(float)*factor);
    for(int kk = 0; kk<factor; kk++)
    {
        K[kk] = 1;
    }

    /* Convolution */
    fim_convn1(B, nB, nB, 1,
              K, factor, 0, 0);
    fim_convn1(B, nB, nB, 1,
              K, factor, 1, 0);


    /* Subsample */
    for(int ll = 0; ll<nA; ll++)
    {
        for(int kk = 0; kk<nA; kk++)
        {
            A[kk + nA*ll] = B[kk*factor+mf + nB*ll*factor +nB+mf];
        }
    }
    free(K);
}

fim_t * gen_psf(opts * s, double lambda)
    {
    if(s->verbose > 2) { printf("Allocating for output\n"); fflush(stdout); }

    double optical_dx = s->optical.dx/s->oversampling;

    fim_tiff_init();
    fim_t * PSF = fimt_zeros(s->M*s->oversampling, s->M*s->oversampling, s->P);

    double Fs = 1/optical_dx;
    double Fn = Fs/2;

    /* Larger lateral grid than PSF */
    int grid_factor = s->xgrid;
    size_t XM = PSF->M*grid_factor;

    /* Phase of Phase Propagator */
    if(s->verbose > 2) { printf("Setting up phase propagator\n"); fflush(stdout); }
    double * ph = malloc(XM*XM*sizeof(double));
    memset(ph, 0, XM*XM*sizeof(double));
    for(int kk = 0; kk< (int) XM; kk++)
    {
        for(int ll = 0; ll< (int) XM; ll++)
        {
            /* Max radius in pixels */
            float rmax = ( (float)XM-1.0)/2.0;
            float dx = (float) kk - rmax;
            float dy = (float) ll - rmax;
            float r = Fn*sqrt(pow(dx/rmax, 2) + pow(dy/rmax, 2));
            float f = pow(s->optical.ni/lambda, 2) - pow(r, 2);
            if(f > 0)
            {
                ph[kk+XM*ll] = sqrt(f);
            }
        }
    }
    ifftshift2_double(ph, XM);

    if(s->verbose > 2) { printf("Setting up pupil function\n"); fflush(stdout); }

    /* Pupil function */
    fftw_complex * pu = fftw_malloc(XM*XM*sizeof(fftw_complex));

    double pur2 = pow(s->optical.NA/lambda, 2);
    for(int kk = 0; kk< (int) XM; kk++)
    {
        for(int ll = 0; ll< (int) XM; ll++)
        {
            float rmax = ( (float)XM-1.0)/2.0;
            float dx = (float) kk - rmax;
            float dy = (float) ll - rmax;
            float r2 = pow(Fn,2)*(pow(dx/rmax, 2) + pow(dy/rmax, 2));

            if(r2 <= pur2 )
            {
                pu[kk+XM*ll][0] = 1;
                pu[kk+XM*ll][1] = 0;
            } else {
                pu[kk+XM*ll][0] = 0;
                pu[kk+XM*ll][1] = 0;
            }
        }
    }


    ifftshift2_complex(pu, XM);

    /* TODO: apply moments */

    if(s->verbose > 2) { printf("Calculating\n"); fflush(stdout); }

    fftw_complex * h = fftw_malloc(XM*XM*sizeof(fftw_complex));
    fftw_complex * H = fftw_malloc(XM*XM*sizeof(fftw_complex));
    fftw_plan plan = fftw_plan_dft_2d(XM, XM,
                                      h, /* In */
                                      H, /* Out */
                                      FFTW_FORWARD, FFTW_ESTIMATE);

    printf("\n");
    for(int zz = 0; zz < (int) PSF->P; zz++)
    {
        double z = (zz - ((int) PSF->P-1)/2)*s->optical.dz;
        if(s->verbose > 0)
        {
            printf("\rplane %d/%zu (z = %.0f nm)        ", zz+1, PSF->P, z);
            fflush(stdout);
        }
        for(size_t kk = 0; kk < XM*XM; kk++)
        {
            double x = 2*M_PI*ph[kk]*z;
            h[kk][0] = pu[kk][0]*cos(x) - pu[kk][1]*sin(x);
            h[kk][1] = pu[kk][0]*sin(x) + pu[kk][1]*cos(x);
        }

        fftw_execute(plan);

        fftshift2_complex(H, XM);
        /* Copy the central MxM region to the PSF */
        int skip = (XM - PSF->M)/2;
        for(int kk = 0; kk < (int) PSF->M; kk++)
        {
            for(int ll = 0; ll < (int) PSF->N; ll++)
            {
                int kkx = kk+skip;
                int llx = ll+skip;
                PSF->V[zz*PSF->M*PSF->N + kk + ll*PSF->M] =
                    pow(H[kkx + XM*llx][0], 2) + pow(H[kkx + XM*llx][1], 2);
            }
        }
    }
    if(s->verbose > 0)
    {
        printf("\r done.                         ");
    }
    fftw_destroy_plan(plan);

    if(s->oversampling > 1)
    {
        float * Vout = malloc(s->M*s->M*s->P*sizeof(double));
        for(int pp = 0; pp< (int) PSF->P; pp++)
        {
            downsample_integrate(Vout + pp*s->M*s->M,
                       PSF->V + pp*PSF->M*PSF->M,
                       s->M, PSF->M);
        }
        free(PSF->V);
        PSF->V = Vout;
        PSF->M /= s->oversampling;
        PSF->N /= s->oversampling;
    }


    return PSF;
    }

static double max_double(double a, double b)
{
    if(a>b)
        return a;
    return b;
}

static double min_double(double a, double b)
{
    if(a<b)
        return a;
    return b;
}
double square_overlap(double x0, double x1, double y0, double y1,
                      double X0, double X1, double Y0, double Y1)
{
    /* Overlap between two squares, first covers [x0,x1]x[y0,y1] */


    double sx = max_double(x0, X0);
    double ex = min_double(x1, X1);
    double sy = max_double(y0, Y0);
    double ey = min_double(y1, Y1);

    if( (sx < ex) & (sy < ey) )
    {
        double AI = (ex-sx)*(ey-sy);

        return AI;
    } else {
        return 0;
    }
}

/* Convolve two images of the same size */
static float * conv2d_float(float * A, float * B,
                            size_t M, size_t N)
{
    fftwf_complex * FA = fftw_malloc(M*N*sizeof(fftwf_complex));
    fftwf_complex * FB = fftw_malloc(M*N*sizeof(fftwf_complex));

    fftwf_plan plan1 = fftwf_plan_dft_r2c_2d(M, N,
                                      A, /* In */
                                      FA, /* Out */
                                      FFTW_ESTIMATE);
    fftwf_execute(plan1);
    fftwf_destroy_plan(plan1);
    fftwf_plan plan2 = fftwf_plan_dft_r2c_2d(M, N,
                                             B, /* In */
                                             FB, /* Out */
                                             FFTW_ESTIMATE);
    fftwf_execute(plan2);
    fftwf_destroy_plan(plan2);

    for(size_t kk = 0; kk<M*N; kk++)
    {
        float ra = FA[kk][0];
        float ca = FA[kk][1];
        float rb = FB[kk][0];
        float cb = FB[kk][1];
        FA[kk][0] = ra*rb - ca*cb;
        FA[kk][1] = ra*cb + rb*ca;
    }

    float * out = fftwf_malloc(M*N*sizeof(float));
    fftwf_plan plan3 = fftwf_plan_dft_c2r_2d(M, N,
                                             FA, /* In */
                                             out, /* Out */
                                             FFTW_ESTIMATE);
    fftwf_execute(plan3);
    fftwf_destroy_plan(plan3);
    return out;
}


static void pinhole_convolution(opts * s, fim_t * PSF)
{
    /* Convolve PSF by a square filter */
    double pinhole_nm = s->optical.pinhole*1.22*s->optical.lambda/s->optical.NA;
    double pinhole_px = pinhole_nm / s->optical.dx;
    printf("Pinhole size: %.2f AU / %.0f nm / %.1f pixels\n",
           s->optical.pinhole, pinhole_nm, pinhole_px);
    float * P = malloc(PSF->M*PSF->N*sizeof(float));
    for(int aa = 0; aa < (int) PSF->M; aa++)
    {
        for(int bb = 0; bb < (int) PSF->N; bb++)
        {
            double apos = (double) aa - (PSF->M -1)/2.0;
            double bpos = (double) bb - (PSF->M -1)/2.0;
            double x0 = apos*s->optical.dx - s->optical.dx/2.0;
            double x1 = x0 + s->optical.dx;
            double y0 = bpos*s->optical.dx - s->optical.dx/2.0;
            double y1 = y0 + s->optical.dx;

            P[bb+aa*PSF->M] = square_overlap(-pinhole_nm/2.0, pinhole_nm/2.0,
                                             -pinhole_nm/2.0, pinhole_nm/2.0,
                                             x0, x1, y0, y1);
        }
    }

    //fim_tiff_write_float("pinhole.tif", P, NULL, PSF->M, PSF->N, 1);

    fftshift2_float(P, PSF->M);;
    //fim_tiff_write_float("pinhole_shifted.tif", P, NULL, PSF->M, PSF->N, 1);

    /* We assume that the image close to the edges will not be used */
    for(size_t pp = 0; pp<PSF->P; pp++)
    {
        float * plane = PSF->V + pp*PSF->M*PSF->N;
        float * C = conv2d_float(plane, P, PSF->M, PSF->N);

        for(size_t kk = 0; kk<PSF->M*PSF->N; kk++)
        {
            plane[kk] = C[kk];
        }
        free(C);
    }



    return;
}

static void dw_psf(opts * s)
{
    /* Internal calculations in double precision,
     * output in single precision */
    if(s->verbose > 0)
    {
        printf("Calculating emission PSF\n");
    }

    fim_t * PSF = gen_psf(s, s->optical.lambda);
    if(s->optical.lambda2 != 0)
    {
        /* We will generate a PSF for a confocal microscope with
        * square-shaped pin-hole */

        /* Convolve PSF with the pinhole */

        //printf("Convolving with pinhole\n");
        pinhole_convolution(s, PSF);

        if(s->verbose > 0)
        {
            printf("Calculating excitation PSF\n");
        }
        fim_t * PSF2 = gen_psf(s, s->optical.lambda2);
        //printf("Multiplying PSFs\n");
        for(size_t kk = 0; kk<fimt_nel(PSF); kk++)
        {
            PSF->V[kk] *= PSF2->V[kk];
        }
        fim_free(PSF2);
    }


    if(s->verbose > 0)
    {
        printf("Writing to %s\n", s->outfile);
    }
    float sum = fimt_sum(PSF);
    fim_mult_scalar(PSF->V, fimt_nel(PSF), 1.0/sum);

    ttags * T = ttags_new();
    char * swstring = malloc(1024);
    sprintf(swstring, "deconwolf %s", deconwolf_version);
    ttags_set_software(T, swstring);
    ttags_set_imagesize(T, PSF->M, PSF->N, PSF->P);
    ttags_set_pixelsize(T, s->optical.dx, s->optical.dx, s->optical.dz);
    free(swstring);


    fim_tiff_write_float(s->outfile, PSF->V,
                         T,
                         PSF->M, PSF->N, PSF->P);
    ttags_free(&T);
    fim_free(PSF);
}

int dw_psf_cli(int argc, char ** argv)
{
    /* Get default settings */
    opts * s = opts_new();
    /* Parse command line and open log file etc */
    argparsing(argc, argv, s);
    if(s->verbose > 2)
    {
        printf("Command lined parsed, continuing\n"); fflush(stdout);
    }
    /* Ready to go */
    dw_psf(s);
    /* Clean up */
    opts_free(s);
    return 0;
}

#ifdef STANDALONE
int main(int argc, char ** argv)
{
    return dw_psf_cli(argc, argv);
}
#endif
