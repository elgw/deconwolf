#include "dw_imshift.h"

typedef struct{
    int overwrite;
    int verbose;
    int optpos;
    char * refImage;
} opts;

static opts * opts_new();
static void opts_free(opts * s);
static void usage(__attribute__((unused)) int argc, char ** argv);
static void argparsing(int argc, char ** argv, opts * s);
static int file_exist(char * fname);
int dw_imshift(int argc, char ** argv);


static opts * opts_new()
{
    opts * s = malloc(sizeof(opts));

    s->overwrite = 0;
    s->verbose = 1;
    s->optpos = -1;
    s->refImage = NULL;
    return s;
}

static void opts_free(opts * s)
{
    if(s->refImage != NULL)
    {
        free(s->refImage);
    }
    free(s);
}

static void usage(__attribute__((unused)) int argc, char ** argv)
{
    printf("usage: %s [<options>] -- input1.tif dx dy dz \n", argv[0]);
    printf("Options:\n");
    printf(" --overwrite\n\t Overwrite existing files\n");
    printf(" --help\n\t Show this message\n");
}


static void argparsing(int argc, char ** argv, opts * s)
{
    struct option longopts[] = {
        {"help", no_argument, NULL, 'h'},
        {"overwrite", no_argument, NULL, 'o'},
        {"ref",     required_argument, NULL, 'r'},
        {"verbose", required_argument, NULL, 'v'},
        {NULL, 0, NULL, 0}};
    int ch;
    while((ch = getopt_long(argc, argv, "hor:v:", longopts, NULL)) != -1)
    {
        switch(ch){
        case 'o':
            s->overwrite = 1;
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
            break;
        case 'r':
            s->refImage = strdup(optarg);
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
    return dw_imshift(argc, argv);
}
#endif


int dw_imshift(int argc, char ** argv)
{

    fim_tiff_init();
    opts * s = opts_new();

    argparsing(argc, argv, s);

    if(s->refImage == NULL)
    {
        if(argc < 4)
        {
            usage(argc, argv);
            exit(1);
        }
    } else {
        if(argc < 4)
        {
            usage(argc, argv);
            exit(1);
        }
    }


    char * outFile;
    char * inFile;


    inFile = argv[s->optpos];
    double dx = -atof(argv[s->optpos+1]);
    double dy = -atof(argv[s->optpos+2]);
    double dz = -atof(argv[s->optpos+3]);

    if(!file_exist(inFile))
    {
        printf("Can't open %s!\n", inFile);
        exit(1);
    }

    outFile = malloc(strlen(inFile) + 20);

    char * _dname = strdup(inFile);
    char * dname = dirname(_dname);
    char * _fname = strdup(inFile);
    char * fname = basename(_fname);

    if(s->verbose > 1)
    {
        fprintf(stdout, "Input file: %s\n", inFile);
        fprintf(stdout, "'%s' /  '%s'\n", dname, fname);
    }

    sprintf(outFile, "%s/sh_%s", dname, fname);
    free(_dname);
    free(_fname);
    if(s->verbose > 1)
    {
        fprintf(stdout, "Ouput file: %s\n", outFile);
    }

    if(s->overwrite == 0 && file_exist(outFile))
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

        // TODO: Max projections

        // TODO: xcorr2

        // TODO: argmax

        fim_shift(A, M, N, P, dx, dy, dz);
        fim_tiff_write(outFile, A, NULL, M, N, P);
        free(A);
        free(R);
    }

    free(outFile);
    opts_free(s);
    return 0;
}
