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
    s->M = 181;
    s->P = 183;
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

static void dw_psf(opts * s)
{
    // TODO from here :)
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
