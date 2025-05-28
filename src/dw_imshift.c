#include "dw_imshift.h"

#include <omp.h>

typedef struct{
    int overwrite;
    int verbose;
    int optpos;
    char * refImage;
    char * image;
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
    s->image = NULL;
    s->dx = 0;
    s->dy = 0;
    s->dz = 0;
    s->nthreads = dw_get_threads();
    return s;
}

static void opts_free(opts * s)
{
    free(s->refImage);
    free(s->image);

    free(s);
}

static void usage(__attribute__((unused)) int argc, char ** argv)
{
    printf("Usage: %s [<options>] --image input1.tif\n", argv[0]);
    printf("Shifts a tif stack using linear interpolation according to dx dy \n"
           "and dz. What is outside of the image is interpreted as 0.\n");
    printf("The method is only efficient for small shifts.\n");
    printf("\n");
    printf("Options:\n");
    printf(" --dx dx\n"
           "\tthe shift in the first dimension (non-strided)\n");
    printf("--dy dy\n"
           "\tshift along the 2nd dimension\n");
    printf("--dz dz\n"
           "\tshift along the 3rd dimension\n");
    printf(" --overwrite\n\t Overwrite existing files\n");
    printf(" --help\n\t Show this message\n");
    printf(" --verbose v\n"
           "\tset verbosity level\n");
    printf(" --threads t\n"
           "\tnumber of threads to use\n");
    printf("--ref ref.tif\n"
           "\tspecify a reference image to align with using normalized\n"
           "\tcross correlation\n");
    // TODO: does the normalized cross correlation feature work? Then
    // document that as well.
}


static void argparsing(int argc, char ** argv, opts * s)
{
    struct option longopts[] = {
        {"help", no_argument, NULL, 'h'},
        {"image", required_argument, NULL, 'i'},
        {"overwrite", no_argument, NULL, 'o'},
        {"ref",     required_argument, NULL, 'r'},
        {"threads", required_argument, NULL, 't'},
        {"verbose", required_argument, NULL, 'v'},
        {"dx", required_argument, NULL, 'x'},
        {"dy", required_argument, NULL, 'y'},
        {"dz", required_argument, NULL, 'z'},
        {NULL, 0, NULL, 0}};

    int ch;
    while((ch = getopt_long(argc, argv, "hi:or:v:x:y:z:", longopts, NULL)) != -1)
    {
        switch(ch){
        case 'o':
            s->overwrite = 1;
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
            break;
        case 'i':
            s->image = strdup(optarg);
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
    return;
}


#ifdef STANDALONE
int main(int argc, char ** argv)
{
    return dw_imshift(argc, argv);
}
#endif


int dw_imshift(int argc, char ** argv)
{

    fim_tiff_init();
    opts * s = opts_new();

    argparsing(argc, argv, s);

    if(s->image == NULL)
    {
        printf("No image specified (use the --image argument)\n");
        exit(EXIT_FAILURE);
    }

    myfftw_start(s->nthreads, s->verbose, stdout);
#ifdef _OPENMP
    omp_set_num_threads(s->nthreads);
#endif

    char * inFile = s->image;
    char * outFile = NULL;

    double dx = s->dx;
    double dy = s->dy;
    double dz = s->dz;

    if(!dw_isfile(inFile))
    {
        printf("Can't open %s!\n", inFile);
        exit(1);
    }

    outFile = malloc(strlen(inFile) + 20);

    char * dname = dw_dirname(inFile);
    char * fname = dw_basename(inFile);

    if(s->verbose > 1)
    {
        fprintf(stdout, "Input file: %s\n", inFile);
        fprintf(stdout, "'%s' /  '%s'\n", dname, fname);
    }

    sprintf(outFile, "%s/sh_%s", dname, fname);
    free(dname);
    free(fname);
    if(s->verbose > 1)
    {
        fprintf(stdout, "Ouput file: %s\n", outFile);
    }

    if(s->overwrite == 0 && dw_isfile(outFile))
    {
        printf("%s exists, skipping.\n", outFile);
        exit(EXIT_SUCCESS);
    }
    if(s->refImage == NULL)
    {
        printf("%s -> %s\n", inFile, outFile);
        int64_t M = 0, N = 0, P = 0;
        float * A = fim_tiff_read(inFile, NULL, &M, &N, &P, s->verbose);
        fim_shift(A, M, N, P, dx, dy, dz);
        fim_tiff_write(outFile, A, NULL, M, N, P);
        free(A);
    }
    if(s->refImage != NULL)
    {
        printf("%s -> %s\n", inFile, outFile);
        printf("Using reference image: %s\n", s->refImage);
        int64_t M = 0, N = 0, P = 0;
        printf("reading images\n");
        float * A = fim_tiff_read(inFile,
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
            printf("Writing to %s\n", outFile);
        }
        fim_tiff_write(outFile, A, NULL, M, N, P);
        free(A);

    }

    free(outFile);
    opts_free(s);
    myfftw_stop();
    return 0;
}
