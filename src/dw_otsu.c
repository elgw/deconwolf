#include "dw_otsu.h"

typedef struct{
    int overwrite;
    int verbose;
    int optpos;
    char * image;
    char * out;
    int nthreads;
} opts;

static opts * opts_new();
static void opts_free(opts * s);
static void usage(__attribute__((unused)) int argc, char ** argv);
static void argparsing(int argc, char ** argv, opts * s);
static int file_exist(char * fname);

static int dw_get_threads(void)
{
    int nThreads = 4;
#ifndef WINDOWS
/* Reports #threads, typically 2x#cores */
    nThreads = sysconf(_SC_NPROCESSORS_ONLN)/2;
#endif
#ifdef OMP
/* Reports number of cores */
    nThreads = omp_get_num_procs();
#endif
    return nThreads;
}


static opts * opts_new()
{
    opts * s = malloc(sizeof(opts));

    s->overwrite = 0;
    s->verbose = 1;
    s->optpos = -1;
    s->image = NULL;
    s->out = NULL;
    s->nthreads = dw_get_threads();
    return s;
}

static void opts_free(opts * s)
{
    if(s->out != NULL)
    {
        free(s->out);
    }
    if(s->image != NULL)
    {
        free(s->image);
    }
    free(s);
}

static void usage(__attribute__((unused)) int argc, char ** argv)
{
    printf("usage: %s [<options>] --file input.tif --out output.tif \n", argv[0]);
    printf("Options:\n");
    printf(" --overwrite\n\t Overwrite existing files\n");
    printf(" --help\n\t Show this message\n");
}


static void argparsing(int argc, char ** argv, opts * s)
{
    struct option longopts[] = {
        {"file", required_argument, NULL, 'f'},
        {"help", no_argument, NULL, 'h'},
        {"image", required_argument, NULL, 'i'},
        {"overwrite", no_argument, NULL, 'o'},
        {"out",     required_argument, NULL, 'p'},
        {"threads", required_argument, NULL, 't'},
        {"verbose", required_argument, NULL, 'v'},
        {NULL, 0, NULL, 0}};
    int ch;
    while((ch = getopt_long(argc, argv, "hi:op:r:v:", longopts, NULL)) != -1)
    {
        switch(ch){
        case 'f':
            s->image = strdup(optarg);
            break;
        case 'o':
            s->overwrite = 1;
            break;
        case 'p':
            s->out = strdup(optarg);
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
            break;
        case 'i':
            s->image = strdup(optarg);
            break;
        case 't':
            s->nthreads = atoi(optarg);
            break;
        case 'v':
            s->verbose = atoi(optarg);
            break;
        }
    }
    s->optpos = optind;
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


#ifdef STANDALONE
int main(int argc, char ** argv)
{
    return dw_otsu(argc, argv);
}
#endif


int dw_otsu(int argc, char ** argv)
{

    fim_tiff_init();
    opts * s = opts_new();

    argparsing(argc, argv, s);

    if(s->image == NULL)
    {
        printf("No --image specified\n");
        exit(EXIT_FAILURE);
    }

    if(s->out == NULL)
    {
        printf("No --out specified\n");
        exit(EXIT_FAILURE);
    }

    omp_set_num_threads(s->nthreads);

    char * inFile = s->image;
    char * outFile = s->out;

    if(!file_exist(inFile))
    {
        printf("Can't open %s!\n", inFile);
        exit(1);
    }

    if(s->verbose > 1)
    {
        fprintf(stdout, "Input file: %s\n", inFile);
    }

    if(s->verbose > 1)
    {
        fprintf(stdout, "Ouput file: %s\n", outFile);
    }

    if(s->overwrite == 0 && file_exist(outFile))
    {
        printf("%s exists, skipping.\n", outFile);
        exit(EXIT_SUCCESS);
    }

    printf("%s -> %s\n", inFile, outFile);
    int64_t M = 0, N = 0, P = 0;
    float * A = fim_tiff_read(inFile, NULL, &M, &N, &P, s->verbose);
    if(s->verbose > 0)
    {
        printf("Max projection\n");
    }
    float * Mproj = fim_maxproj(A, M, N, P);
    free(A);
    if(s->verbose > 0)
    {
        printf("Thresholding\n");
    }
    float * B = fim_otsu(Mproj, M, N);
    free(Mproj);
    if(s->verbose > 0)
    {
        printf("Labelling\n");
    }

    int * L = fim_conncomp6(B, M, N);
    float * fL = malloc(M*N*sizeof(float));
    for(size_t kk = 0; kk< (size_t) M*N; kk++)
    {
        fL[kk] = L[kk];
    }
    free(L);

    fim_tiff_write_noscale(outFile, fL, NULL, M, N, 1);

    free(fL);

    opts_free(s);
    return 0;
}
