#include "dw_psf.h"

typedef struct{
    double NA;
    double ni;
    double dx;
    double dz;
    double lambda;
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
    s->nthreads = dw_get_threads();
    s->optical.NA = 1.45;
    s->optical.ni = 1.515;
    s->optical.dx = 65;
    s->optical.dz = 200;
    s->optical.lambda = 460;
    s->M = 7;
    s->P = 7;
    return s;
}


static void opts_free(opts * s)
{
    dw_nullfree(s->outfile);
    dw_nullfree(s->logfile);
    fclose(s->log);
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
    fprintf(f, "lambda=%f nm\n", s->optical.lambda);
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
    printf("General:\n");
    printf(" --overwrite\n\t Overwrite existing files\n");
    printf(" --help\n\t Show this message\n");
    printf(" --notext\n\t Don't write any text on the image\n");
    printf(" --verbose v\n\t Verbosity level\n");
    printf("\n");
    opts_free(s);
    return;
}

static void argparsing(int argc, char ** argv, opts * s)
{
    struct option longopts[] = {
        {"lambda", required_argument, NULL, 'l'},
        {"NA",     required_argument, NULL, 'n'},
        {"dx",     required_argument, NULL, 'x'},
        {"dz",     required_argument, NULL, 'z'},
        {"ni",     required_argument, NULL, 'i'},
        {"help", no_argument, NULL, 'h'},
        {"overwrite", no_argument, NULL, 'o'},
        {"threads", required_argument, NULL, 't'},
        {"verbose", required_argument, NULL, 'v'},
        {NULL, 0, NULL, 0}};
    int ch;
    while((ch = getopt_long(argc, argv,
                            "l:n:x:z:i:hov:v:", longopts, NULL)) != -1)
    {
        switch(ch){
        case 'h':
            usage(argc, argv);
            exit(0);
            break;
        case 'i':
            s->optical.ni = atof(optarg);
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
            break;
        case 'z':
            s->optical.dz = atof(optarg);
            break;
        case 'n':
            s->optical.NA = atof(optarg);
            break;
        }
    }

    if(optind == argc)
    {
        s->outfile = malloc(1024);
        sprintf(s->outfile,
                "psf_NA%.2f_ni%.2f_lambda%.0f_dx%.0f_dz%.0f.tif",
                s->optical.NA, s->optical.ni, s->optical.lambda,
                s->optical.dx, s->optical.dz);
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

    // Shift in dimension 1
    for(int kk = 0; kk<M; kk++)
    {
        shift_double(ph+kk*M, M, 1, (M-1)/2);
    }
    // Shift in dimension 2
    for(int kk = 0; kk<M; kk++)
    {
        shift_double(ph+kk, M, M, (M-1)/2);
    }
}

void ifftshift2_complex(fftw_complex * ph, int M)
{
    // Shift in dimension 1
    for(int kk = 0; kk<M; kk++)
    {
        shift_complex(ph+kk*M, M, 1, (M-1)/2);
    }
    // Shift in dimension 2
    for(int kk = 0; kk<M; kk++)
    {
        shift_complex(ph+kk, M, M, (M-1)/2);
    }

}

void fftshift2_complex(fftw_complex * ph, int M)
{
    // Shift in dimension 1
    for(int kk = 0; kk<M; kk++)
    {
        shift_complex(ph+kk*M, M, 1, (M-1)/2);
    }
    // Shift in dimension 2
    for(int kk = 0; kk<M; kk++)
    {
        shift_complex(ph+kk, M, M, (M-1)/2);
    }

}


static void dw_psf(opts * s)
{
    /* Internal calculations in double precision,
     * output in single precision */
    if(s->verbose > 2) { printf("Allocating for output\n"); fflush(stdout); }

    fim_tiff_init();
    fim_t * PSF = fimt_zeros(s->M, s->M, s->P);


    double Fs = 1/s->optical.dx;
    double Fn = Fs/2;

    /* Larger lateral grid than PSF */
    int grid_factor = 1;
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
            float f = pow(s->optical.ni/s->optical.lambda, 2) - pow(r, 2);
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

    double pur2 = pow(s->optical.NA/s->optical.lambda, 2);
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
    for(int zz = 0; zz < (int) PSF->P; zz++)
    {
        double z = (zz - ((int) PSF->P-1)/2)*s->optical.dz;
        printf("z = %f\n", z);
        for(size_t kk = 0; kk < XM*XM; kk++)
        {
            double x = 2*M_PI*ph[kk]*z;
            h[kk][0] = pu[kk][0]*cos(x) - pu[kk][1]*sin(x);
            h[kk][1] = pu[kk][0]*sin(x) + pu[kk][1]*cos(x);
        }


        fftw_plan plan = fftw_plan_dft_2d(XM, XM,
                                          h, /* In */
                                          H, /* Out */
                                          FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_execute(plan);
        fftw_destroy_plan(plan);

        fftshift2_complex(H, XM);
        for(int kk = 0; kk < (int) PSF->M; kk++)
        {
            for(int ll = 0; ll < (int) PSF->N; ll++)
            {
                PSF->V[zz*PSF->M*PSF->N + kk + ll*PSF->M] =
                    pow(H[kk + XM*ll][0], 2) + pow(H[kk + XM*ll][1], 2);
            }
        }
    }

    ttags * TTAGS = NULL; /* TODO */

    if(s->verbose > 0)
    {
        printf("Writing to %s\n", s->outfile);
    }
    fim_tiff_write_float(s->outfile, PSF->V,
                         TTAGS,
                         PSF->M, PSF->N, PSF->P);

}

int dw_psf_cli(int argc, char ** argv)
{
    /* Get default settings */
    opts * s = opts_new();
    /* Parse command line and open log file etc */
    argparsing(argc, argv, s);
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
