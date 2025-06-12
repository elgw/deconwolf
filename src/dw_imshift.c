#include <inttypes.h>
#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fim.h"
#include "fim_tiff.h"
#include "dw_version.h"

#include "dw_imshift.h"



typedef struct{
    int overwrite;
    int verbose;
    int optpos;
    char * refImage;
    char * inFile;
    char * outFile;
    double dx;
    double dy;
    double dz;
    int nthreads;
} opts;

static opts * opts_new();
static void opts_free(opts * s);
static void usage(__attribute__((unused)) int argc, char ** argv);
static void argparsing(int argc, char ** argv, opts * s);

int dw_imshift(int argc, char ** argv);

static opts * opts_new()
{
    opts * s = calloc(1, sizeof(opts));
    assert(s != NULL);
    s->overwrite = 0;
    s->verbose = 1;
    s->optpos = -1;
    s->refImage = NULL;
    s->inFile = NULL;
    s->dx = 0;
    s->dy = 0;
    s->dz = 0;
    s->nthreads = dw_get_threads();
    return s;
}

static void opts_free(opts * s)
{
    free(s->refImage);
    free(s->inFile);
    free(s->outFile);
    free(s);
}

static void usage(__attribute__((unused)) int argc, char ** argv)
{
    printf("Usage: %s [<options>] --image input1.tif\n", argv[0]);
    printf("Shifts a tif stack using linear interpolation according to dx dy \n"
           "and dz. What is outside of the image is interpreted as 0.\n");
    printf("\n");
    printf("Options:\n");
    printf("  --dx dx\n"
           "\tthe shift in the first dimension (non-strided)\n");
    printf("  --dy dy\n"
           "\tshift along the 2nd dimension\n");
    printf("  --dz dz\n"
           "\tshift along the 3rd dimension\n");
    printf("  --overwrite\n\t"
           "Overwrite existing files\n");
    printf("  --help\n\t"
           "Show this message\n");
    printf("  --verbose v\n"
           "\tset verbosity level\n");
    printf("  --threads t\n"
           "\tnumber of threads to use\n");
    printf("  --ref ref.tif\n"
           "\tspecify a reference image to align with using normalized\n"
           "\tcross correlation\n");
    printf("  --out file.tif\n\t"
           "Specify output file name\n");
    return;
}


static void argparsing(int argc, char ** argv, opts * s)
{
    struct option longopts[] = {
        {"help", no_argument, NULL, 'h'},
        {"image", required_argument, NULL, 'i'},
        {"overwrite", no_argument, NULL, 'o'},
        {"out", required_argument, NULL, 'O'},
        {"ref",     required_argument, NULL, 'r'},
        {"threads", required_argument, NULL, 't'},
        {"verbose", required_argument, NULL, 'v'},
        {"dx", required_argument, NULL, 'x'},
        {"dy", required_argument, NULL, 'y'},
        {"dz", required_argument, NULL, 'z'},
        {NULL, 0, NULL, 0}};

    int ch;
    while((ch = getopt_long(argc, argv, "hi:oO:r:v:x:y:z:", longopts, NULL)) != -1)
    {
        switch(ch){
        case 'o':
            s->overwrite = 1;
            break;
        case 'O':
            free(s->outFile);
            s->outFile = strdup(optarg);
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
            break;
        case 'i':
            s->inFile = strdup(optarg);
            break;
        case 'r':
            s->refImage = strdup(optarg);
            break;
        case 't':
            s->nthreads = atoi(optarg);
            break;
        case 'v':
            s->verbose = atoi(optarg);
            break;
        case 'x':
            s->dx = atof(optarg);
            break;
        case 'y':
            s->dy = atof(optarg);
            break;
        case 'z':
            s->dz = atof(optarg);
            break;
        }
    }
    s->optpos = optind;

#ifdef _OPENMP
    if(s->nthreads > 0)
    {
        omp_set_num_threads(s->nthreads);
    }
#endif

    return;
}


#ifdef STANDALONE
int main(int argc, char ** argv)
{
    return dw_imshift(argc, argv);
}
#endif

static void perform_imshift(const opts * s){
    fim_tiff_info info = {0};
    if(fim_tiff_get_info(s->inFile, &info))
    {
        fprintf(stderr, "Failed to read %s\n", s->inFile);
        return;
    }

    int64_t M = 0, N = 0, P = 0;
    float * A = fim_tiff_read(s->inFile, NULL, &M, &N, &P, s->verbose);
    if(A == NULL)
    {
        fprintf(stderr, "Failed to read %s\n", s->inFile);
        return;
    }
    fim_shift(A, M, N, P, s->dx, s->dy, s->dz);

    if(info.BPS <= 16)
    {
        fim_tiff_write_opt(s->outFile, A,
                           NULL,
                           M, N, P,
                           1.0);
    } else {
        fim_tiff_write_float(s->outFile, A,
                             NULL,
                             M, N, P);
    }
    free(A);
}

static void perform_xcorr(opts * s){
    printf("Using reference image: %s\n", s->refImage);
    int64_t M = 0, N = 0, P = 0;
    printf("reading images\n");
    float * A = fim_tiff_read(s->inFile,
                              NULL, &M, &N, &P, s->verbose);
    int64_t rM = 0, rN = 0, rP = 0;
    float * R = fim_tiff_read(s->refImage,
                              NULL, &rM, &rN, &rP,
                              s->verbose);

    if(M!=rM || N != rN || P != rP)
    {
        printf("Image dimensions does not match\n");
        exit(EXIT_FAILURE);
    }

    if(s->verbose > 1)
    {
        printf("Image size: [%" PRId64 ", %" PRId64 ", %" PRId64 "]\n", M, N, P);
        printf("Creating max projections\n");
    }
    //float * mA = fim_maxproj(A, M, N, P);
    //float * mR = fim_maxproj(R, M, N, P);
    float * mA = fim_sumproj(A, M, N, P);
    float * mR = fim_sumproj(R, M, N, P);
    fim_tiff_write_float("mA.tif", mA, NULL, M, N, 1);
    free(R);

    if(s->verbose > 1)
    {
        printf("Calculating normalized cross correlation\n");
    }
    float * XC = fim_xcorr2(mA, mR, M, N);
    free(mA);
    free(mR);

    //fim_tiff_write_float("XC.tif", XC, NULL, 2*M-1, 2*N-1, 1);

    int64_t aM=0;
    int64_t aN=0;
    int64_t aP = 0;

    fim_argmax(XC, 2*M-1, 2*N-1, 1, &aM, &aN, &aP);

    if(s->verbose > 1)
    {
        printf("Found max at (%" PRId64 ", %" PRId64 ", %" PRId64 ")\n",
               aM, aN, aP);
    }
    s->dx = (float) M-aM-1.0;
    s->dy = (float) N-aN-1.0;
    s->dz = 0;
    if(s->verbose > 0)
    {
        float corr = fim_max(XC, (2*M-1)*(2*N-1));
        printf("Normalized cross correlation: %f\n", corr);
        printf("Shifting image by %f, %f, %f\n", s->dx, s->dy, s->dz);
    }
    free(XC);
    fim_shift(A, M, N, P, s->dx, s->dy, s->dz);
    if(s->verbose > 0)
    {
        printf("Writing to %s\n", s->outFile);
    }
    fim_tiff_write(s->outFile, A, NULL, M, N, P);
    free(A);
    return;
}

int dw_imshift(int argc, char ** argv)
{

    fim_tiff_init();
    opts * s = opts_new();

    argparsing(argc, argv, s);

    if(s->inFile == NULL)
    {
        printf("No image specified (use the --image argument)\n");
        exit(EXIT_FAILURE);
    }

    if(!dw_isfile(s->inFile))
    {
        printf("Can't open %s!\n", s->inFile);
        exit(1);
    }

    if(s->outFile == NULL)
    {
        s->outFile = dw_prefix_file(s->inFile, "sh");
    }

    if(s->verbose > 1)
    {
        fprintf(stdout, "Output file: %s\n", s->outFile);
    }

    if(s->overwrite == 0 && dw_isfile(s->outFile))
    {
        printf("%s exists, skipping.\n", s->outFile);
        exit(EXIT_SUCCESS);
    }

    if(s->verbose > 0)
    {
        printf("%s -> %s\n", s->inFile, s->outFile);
    }

    if(s->refImage == NULL)
    {
        perform_imshift(s);
    } else {
        myfftw_start(s->nthreads, s->verbose, stdout);
        perform_xcorr(s);
        myfftw_stop();
    }

    opts_free(s);

    return 0;
}
