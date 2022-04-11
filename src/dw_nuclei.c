#include "dw_nuclei.h"

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
    int ntree; /* Number of trees in random forest classifier */
} opts;

static opts * opts_new();
static void opts_free(opts * s);
static void usage(__attribute__((unused)) int argc, char ** argv);
static void argparsing(int argc, char ** argv, opts * s);

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
    s->ntree = 50;
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


float * transpose(const float * X, size_t M, size_t N)
{
    float * Y = malloc(M*N*sizeof(float));
    size_t pos = 0;
    for(size_t nn = 0; nn<N; nn++)
    {
        for(size_t mm = 0; mm<M; mm++)
        {
            Y[pos++] = X[nn+mm*N];
        }
    }
    return Y;
}

/* Append a column with nrow elements to a table with nrow*ncol values
 * A should not be freed after this call as it is replaced by the
 * returned pointer */
float * cm_append_column(float * A, size_t nrow, size_t ncol, const float * B)
{
    A = realloc(A, nrow*(ncol+1)*sizeof(float));
    float * C = A+nrow*ncol; /* last, still unset column */
    for(size_t kk = 0; kk<nrow; kk++)
    {
        C[kk] = B[kk];
    }
    return A;
}

float * subset_cm(float * A,
                  float * V,
                  size_t nrow, size_t ncol,
                  size_t * nrow_out, size_t * ncol_out)
{
    /* Count how many of the rows to use */
    size_t rows = 0;
    for(size_t kk = 0; kk<nrow; kk++)
    {
        V[kk] > 0 ? rows++ : 0;
    }

    if(rows == 0)
    {
        *nrow_out = 0;
        *ncol_out = 0;
        return NULL;
    }

    float * S = malloc(rows*ncol*sizeof(float));
    size_t writepos = 0;
    for(size_t cc = 0; cc<ncol; cc++)
    {
        for(size_t rr = 0; rr<nrow; rr++)
        {
            V[rr] > 0 ? S[writepos++] = A[cc*nrow + rr] : 0;
        }
    }
    *nrow_out = rows;
    *ncol_out = ncol;
    return S;
}


void random_forest_pipeline(opts * s)
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

    /* Transpose to column-major */
    float * features_cm = transpose(features->T,
                                    features->nrow, features->ncol);

    /* Append one columns for the annotations */
    features_cm = cm_append_column(features_cm, features->nrow, features->ncol, anno->V);
    size_t cols = features->ncol+1;

    /* Extract the data there anno->V > 0 to a separate array */
    size_t ncol_train = 0;
    size_t nrow_train = 0;
    float * features_cm_train = subset_cm(features_cm,
                                          anno->V,
                                          features->nrow, cols,
                                          &nrow_train, &ncol_train);


    /* Train classifier */
    // Use prf_forest from pixel_random_forest
    PrfForest * F = prf_forest_new(s->ntree);
    F->nthreads = s->nthreads;

    if(prf_forest_train(F, features_cm_train, nrow_train, ncol_train))
    {
        printf("dw_nuclei: Failed to train the random forest\n");
        exit(EXIT_FAILURE);
    }

    /* Now apply it to all pixels */
    int * class = prf_forest_classify_table(F, features_cm, anno->M*anno->N,
                                            features->ncol);
    prf_forest_free(F);

    float * result = malloc(anno->M*anno->N*sizeof(float));
    for(size_t kk = 0; kk < anno->M*anno->N; kk++)
    {
        if(class[kk] == 2)
        {
            result[kk] = 1;
        } else {
            result[kk] = 0;
        }
    }
    fim_tiff_write_noscale("training_classification.tif", result, NULL,
                           anno->M, anno->N, 1);
    return;
}


#ifdef STANDALONE
int main(int argc, char ** argv)
{
    return dw_nuclei(argc, argv);
}
#endif


int dw_nuclei(int argc, char ** argv)
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

    if(!dw_file_exist(inFile))
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

    if(s->overwrite == 0 && dw_file_exist(outFile))
    {
        printf("%s exists, skipping.\n", outFile);
        exit(EXIT_SUCCESS);
    }

    printf("%s -> %s\n", inFile, outFile);
    int64_t M = 0, N = 0, P = 0;
    float * A = fim_tiff_read(inFile, NULL, &M, &N, &P, s->verbose);
    if(s->anno_label != NULL && s->anno_image != NULL)
    {
        random_forest_pipeline(s);
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
