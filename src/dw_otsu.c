#include "dw_otsu.h"

typedef struct{
    int overwrite;
    int verbose;
    int optpos;
    char * image;
    char * out;
    /* If annotated data exist, this should be used for training */
    char * anno_label;
    char * anno_image;
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
    s->anno_image = NULL;
    s->anno_label = NULL;
    return s;
}

static void nullfree(void * p)
{
    if(p != NULL)
    {
        free(p);
    }
}

static void opts_free(opts * s)
{
    if(s == NULL)
    {
        return;
    }
    nullfree(s->out);
    nullfree(s->image);
    nullfree(s->anno_image);
    nullfree(s->anno_label);
    free(s);
}


/* Read png image, when more green that red, set to 1, when more red
 * than green set to 2. Everything else set to 0.
 * http://www.libpng.org/pub/png/libpng-manual.txt
 */
fim_t * fim_png_read_green_red(char * fname)
{
    fim_t * F = NULL;

    png_image image;
    memset(&image, 0, sizeof(image));
    image.version = PNG_IMAGE_VERSION;
    if(png_image_begin_read_from_file(&image, fname))
    {
        png_bytep buffer;
        image.format = PNG_FORMAT_RGB;
        buffer = malloc(PNG_IMAGE_SIZE(image));
        if(buffer == NULL)
        {
            fprintf(stderr, "Unable to allocate memory for the image buffer\n");
            exit(EXIT_FAILURE);
        }
        if(png_image_finish_read(&image, NULL, buffer, 0, NULL))
        {
            size_t M = image.width;
            size_t N = image.height;
            printf("M=%zu, N=%zu image buffer %u b\n", M, N, PNG_IMAGE_SIZE(image));

            F = malloc(sizeof(fim_t));
            F->M = M;
            F->N = N;
            F->P = 1;
            F->V = malloc(M*N*sizeof(float));
            printf("Reading image...\n"); fflush(stdout);
            size_t nfg = 0;
            size_t nbg = 0;
            for(size_t kk = 0; kk<M*N; kk++)
            {
                F->V[kk] = 0;
                if(buffer[3*kk+1] > buffer[3*kk+0])
                {
                    F->V[kk] = 1;
                    nfg++;
                }
                if(buffer[3*kk+1] < buffer[3*kk+0])
                {
                    F->V[kk] = 2;
                    nbg++;
                }
            }
            printf("%zu foreground (green->2) and %zu background (red->1) pixels\n",
                   nfg, nbg);
            printf("done\n"); fflush(stdout);
        }

        free(buffer);
    }

    png_image_free(&image);
    return F;
}

static void usage(__attribute__((unused)) int argc, char ** argv)
{
    printf("usage: %s [<options>] --file input.tif --out output.tif \n", argv[0]);
    printf("Options:\n");
    printf(" --alabel file.png\n\t Annotated image\n");
    printf(" --aimage file.tif\n\t Annotated raw image\n");
    printf(" --overwrite\n\t Overwrite existing files\n");
    printf(" --help\n\t Show this message\n");
}


static void argparsing(int argc, char ** argv, opts * s)
{
    struct option longopts[] = {
        {"alabel", required_argument, NULL, 'a'},
        {"aimage", required_argument, NULL, 'b'},
        {"file", required_argument, NULL, 'f'},
        {"help", no_argument, NULL, 'h'},
        {"image", required_argument, NULL, 'i'},
        {"overwrite", no_argument, NULL, 'o'},
        {"out",     required_argument, NULL, 'p'},
        {"threads", required_argument, NULL, 't'},
        {"verbose", required_argument, NULL, 'v'},
        {NULL, 0, NULL, 0}};
    int ch;
    while((ch = getopt_long(argc, argv, "a:b:hi:op:r:v:", longopts, NULL)) != -1)
    {
        switch(ch){
        case 'a':
            s->anno_label = strdup(optarg);
            break;
        case 'b':
            s->anno_image = strdup(optarg);
            break;
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
    if(s->anno_label != NULL && s->anno_image != NULL)
    {
        /* Read annotated image */
        fim_t * anno = fim_png_read_green_red(s->anno_label);

        /* Read corresponding raw image */
        int64_t M = 0; int64_t N = 0; int64_t P = 0;
        float * _anno_raw = fim_tiff_read(s->anno_image, NULL, &M, &N, &P, 0);
        fim_t * anno_raw = malloc(sizeof(fim_t));
        anno_raw->V = _anno_raw;
        anno_raw->M = M;
        anno_raw->N = N;
        anno_raw->P = P;
        if(P != 1)
        {
            fprintf(stderr, "The raw annotated image should have only one slice\n");
            exit(EXIT_FAILURE);
        }
        /* Extract features */
        ftab_t * features = fim_features_2d(anno_raw);

        // Insert annotations into table.
        // ftab_insert(features, ftab_get_col("class"), anno->V);

        /* Train classifier */
        // Use prf_forest from pixel_random_forest
        fprintf(stderr, "Pipeline not finished!\n");
        exit(EXIT_FAILURE);
    }
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
