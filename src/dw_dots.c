#include "dw_dots.h"

typedef struct{
    int overwrite;
    int verbose;
    int optpos;
    char * image;
    char * out;
    int nthreads;
    float asigma; /* Axial sigma */
    float lsigma; /* Lateral sigma */
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
    s->lsigma = 1;
    s->asigma = 1;
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
        {"asigma", required_argument, NULL, 'a'},
        {"file", required_argument, NULL, 'f'},
        {"help", no_argument, NULL, 'h'},
        {"image", required_argument, NULL, 'i'},
        {"lsigma", required_argument, NULL, 'l'},
        {"overwrite", no_argument, NULL, 'o'},
        {"out",     required_argument, NULL, 'p'},
        {"threads", required_argument, NULL, 't'},
        {"verbose", required_argument, NULL, 'v'},
        {NULL, 0, NULL, 0}};
    int ch;
    while((ch = getopt_long(argc, argv, "a:hi:l:op:r:v:", longopts, NULL)) != -1)
    {
        switch(ch){
        case 'a':
            s->asigma = atof(optarg);
            break;
        case 'f':
            s->image = strdup(optarg);
            break;
        case 'l':
            s->lsigma = atof(optarg);
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
    return dw_dots(argc, argv);
}
#endif


int dw_dots(int argc, char ** argv)
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
        printf("LoG filter, lsigma=%.2f asigma=%.2f\n", s->lsigma, s->asigma);
    }
    float * LoG = fim_LoG(A, M, N, P, s->lsigma, s->asigma);
    free(A);

    fim_tiff_write_float(outFile, LoG, NULL, M, N, P);

    /* Detect local maxima */
    fim_table_t * T = fim_lmax(LoG, M, N, P);

    printf("Writing dots to dots.tsv\n");
    FILE * fid = fopen("dots.tsv", "w");
    fprintf(fid, "x\ty\tz\tvalue\n");
    float th = 200;
    fim_table_sort(T, 3);
    size_t nwritten = 0;
    for(size_t kk = 0; kk<T->nrow; kk++)
    {
        float * row = T->T + kk*T->ncol;
        if(T->T[kk*T->ncol + 3] > th)
        {
            fprintf(fid, "%f\t%f\t%f\t%f\n", row[0], row[1], row[2], row[3]);
            nwritten++;
        }
    }
    printf("Wrote %zu dots\n", nwritten);
    fclose(fid);
    fim_table_free(T);

    free(LoG);

    opts_free(s);
    return 0;
}
