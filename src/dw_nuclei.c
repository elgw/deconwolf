#include "dw_nuclei.h"

typedef enum {
    REDU_MAX,
    REDU_MEAN,
    REDU_FOCUS
} reduction_type;

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
    int train_loop;
    reduction_type redu;
} opts;

static opts * opts_new();
static void opts_free(opts * s);
static void opts_print(FILE * f, opts * s);
static void usage(__attribute__((unused)) int argc, char ** argv);
static void argparsing(int argc, char ** argv, opts * s);

static opts * opts_new()
{
    opts * s = malloc(sizeof(opts));
    assert(s != NULL);

    s->overwrite = 0;
    s->verbose = 1;
    s->optpos = -1;
    s->image = NULL;
    s->out = NULL;
    s->nthreads = dw_get_threads();
    s->anno_image = NULL;
    s->anno_label = NULL;
    s->ntree = 50;
    s->train_loop = 0;
    s->redu = REDU_FOCUS;
    return s;
}

static void opts_print(FILE * f, opts * s)
{
    fprintf(f, "Overwrite: %d\n", s->overwrite);
    switch(s->redu)
    {
    case REDU_MAX:
        fprintf(f, "Reduction: Max projection\n");
        break;
    case REDU_FOCUS:
        fprintf(f, "Reduction: Most in focus slice\n");
        break;
    case REDU_MEAN:
        fprintf(f, "Reduction: Mean projection\n");
        break;
    }
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
fimo * fim_png_read_green_red(char * fname)
{
    fimo * F = NULL;

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
            //printf("M=%zu, N=%zu image buffer %u b\n", M, N, PNG_IMAGE_SIZE(image));

            F = malloc(sizeof(fimo));
            assert(F != NULL);
            F->M = M;
            F->N = N;
            F->P = 1;
            F->V = malloc(M*N*sizeof(float));
            assert(F->V != NULL);
            printf("Reading %s\n", fname); fflush(stdout);
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
            printf("Found %zu foreground (green->2) and %zu background (red->1) pixels\n",
                   nfg, nbg);
        }

        free(buffer);
    }

    png_image_free(&image);
    return F;
}

static void usage(__attribute__((unused)) int argc, char ** argv)
{
    printf("usage: %s [<options>] --aimage input.tif --alabel labels.png file1.tif ... \n", argv[0]);
    printf("Options:\n");
    printf(" --aimage file.tif\n\t (raw) annotated image\n");
    printf(" --alabel file.png\n\t "
           "Annotated image where nuclei is marked green, background is marked red\n");
    printf(" --overwrite\n\t Overwrite existing files\n");
    printf(" --loop\n\t Enter training loop\n");
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
        {"loop", no_argument, NULL, 'l'},
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
        case 'l':
            s->train_loop = 1;
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
    assert(Y != NULL);
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

/* Return a newly allocated array AA = [A ; B] */
float * cm_append_column(float * A, size_t nrow, size_t ncol, const float * B)
{
    float * AA = malloc(nrow*(ncol+1)*sizeof(float));
    assert(AA != NULL);
    memcpy(AA, A,
           nrow*ncol*sizeof(float));
    memcpy(AA+nrow*ncol, B, nrow*sizeof(float));
    return AA;
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
    assert(S != NULL);
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

static int float_arg_max(const float * v, size_t N)
{
    float max = v[0];
    int argmax = 0;
    for(size_t kk = 0; kk<N; kk++)
    {
        if(v[kk] > max)
        {
            max = v[kk];
            argmax = kk;
        }
    }
    return argmax;
}

fimo * get_reduction(opts * s, char * file)
{


    int64_t M = 0; int64_t N = 0; int64_t P = 0;
    if(s->verbose > 0)
    {
        printf("Reading tif: %s\n", file);
    }
    float * I = fim_tiff_read(file, NULL, &M, &N, &P, 0);
    if(I == NULL)
    {
        printf("%s could not be read, unrecognizable format.\n", file);
        return NULL;
    }

    /* Determine scaling */
    if(s->verbose > 1)
    {
        printf("Looking for a scaling factor in the dw log file\n");
        fflush(stdout);
    }
    float scaling = dw_read_scaling(file);
    if(s->verbose > 0)
    {
        printf("Scaling by %f\n", 1.0/scaling);
        fflush(stdout);
    }
    fim_mult_scalar(I, M*N*P, 1.0/scaling);

    if(s->redu == REDU_MAX)
    {
        float * maxI = fim_maxproj(I, M, N, P);
        free(I);
        fimo * result = fim_image_from_array(maxI, M, N, 1);
        free(maxI);
        return result;
    }

    if(s->redu == REDU_FOCUS)
    {
        if(s->verbose > 0)
        {
            printf("Finding focus\n"); fflush(stdout);
        }
        fimo * II = fim_image_from_array(I, M, N, P);
        free(I);
        float sigma = 1;
        float * gm = fim_focus_gm(II, sigma);
        int slice = float_arg_max(gm, II->P);
        fimo * result = malloc(sizeof(fimo));
        assert(result != NULL);
        result->M = M;
        result->N = N;
        result->P = 1;
        result->V = malloc(M*N*sizeof(float));
        assert(result->V != NULL);
        memcpy(result->V, II->V+slice*M*N, M*N*sizeof(float));
        fimo_free(II);
        if(s->verbose > 0)
        {
            printf("Returning slice %d\n", slice); fflush(stdout);
        }
        return result;
    }

    if(s->redu == REDU_MEAN)
    {
        fprintf(stderr, "Reduction: Mean projection is not implemented\n");
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "Unknown reduction type\n");
    exit(EXIT_FAILURE);
    return NULL;
}

void segment_image_rf(opts * s, PrfForest * F, char * file)
{
    if(!dw_isfile(file))
    {
        printf("%s does not exist\n", file);
        return;
    }
    char * outfile = malloc(strlen(file) + 32);
    assert(outfile != NULL);
    sprintf(outfile, "%s.mask.tif", file);
    printf("%s -> %s", file, outfile);
    if(dw_isfile(outfile) && s->overwrite == 0)
    {
        printf(" File exist. Skipping\n");
        return;
    }
    printf("\n");

    /* Read the image */
    fimo * redu = get_reduction(s, file);
    assert(redu != NULL);
    printf("1\n"); fflush(stdout);
    /* Extract features */
    ftab_t * features = fim_features_2d(redu);
    assert(features != NULL);
    printf("2\n"); fflush(stdout);
    /* Transpose to column-major */
    float * features_cm = transpose(features->T,
                                    features->nrow, features->ncol);
    printf("3\n"); fflush(stdout);
    /* Classify */
    int * class = prf_forest_classify_table(F, features_cm,
                                            redu->M*redu->N, features->ncol);
    ftab_free(features);
    free(features_cm);
    printf("4\n"); fflush(stdout);
    float * result = malloc(redu->M*redu->N*sizeof(float));
    assert(result != NULL);
    for(size_t kk = 0; kk < (size_t) (redu->M*redu->N); kk++)
    {
        if(class[kk] == 2)
        {
            result[kk] = 0;
        } else {
            result[kk] = 1;
        }
    }
    printf("5\n"); fflush(stdout);
    int fill_holes = 1;
    float max_hole_size = pow(25, 2);
    if(fill_holes)
    {
        printf("Filling holes up to %.1f pixels large\n", max_hole_size);
        float * result2 = fim_fill_holes(result, redu->M, redu->N, max_hole_size);
        free(result);
        result = result2;
    }

    /* remove debris */
    int remove_small = 1;
    /* something like everything < 1/4 of the expected
     *  might make sense */
    float min_size = pow(40, 2);
    printf("Removing objects with < %.1f pixels\n", min_size-1);
    if(remove_small)
    {
        float * result2 = fim_remove_small(result, redu->M, redu->N, min_size);
        free(result);
        result = result2;
    }

    /* TODO: split large? */

    /* TODO: export properties to tsv */

    /* label */
    int * L = fim_conncomp6(result, redu->M, redu->N);
    for(size_t kk = 0; kk<(size_t) (redu->M*redu->N); kk++)
    {
        result[kk] = L[kk];
    }
    free(L);

    fim_tiff_write_noscale(outfile, result, NULL,
                           redu->M, redu->N, 1);
    fimo_free(redu);
    return;
}

PrfForest * loop_training_data(opts * s, float * features_cm,
                               size_t nsamples, size_t nfeatures)
{

    PrfForest * F = NULL;
    int done = 0;
    while(done == 0)
    {
        if(F != NULL)
        {
            free(F);
        }

        /* Read annotated image */
        fimo * anno = fim_png_read_green_red(s->anno_label);
        assert(anno != NULL);

        /* Append one columns for the annotations */

        float * features_cma = cm_append_column(features_cm,
                                        nsamples, nfeatures,
                                        anno->V);
        /* Extract the data there anno->V > 0 to a separate array */
        size_t ncol_train = 0;
        size_t nrow_train = 0;
        float * features_cma_train = subset_cm(features_cma,
                                               anno->V,
                                               nsamples, nfeatures+1,
                                               &nrow_train, &ncol_train);


        /* Train classifier */
        // Use prf_forest from pixel_random_forest
        F = prf_forest_new(200); //s->ntree);
        F->nthreads = s->nthreads;

        if(prf_forest_train(F, features_cma_train, nrow_train, ncol_train))
        {
            printf("dw_nuclei: Failed to train the random forest\n");
            exit(EXIT_FAILURE);
        }

        printf("Validating the training data\n");
        int * class = prf_forest_classify_table(F, features_cma_train, nrow_train,
                                                ncol_train);

        size_t ncorrect = 0;
        for(size_t kk = 0; kk<nrow_train; kk++)
        {
            if(class[kk] == features_cma_train[nrow_train*(ncol_train-1)+kk])
            {
                ncorrect++;
            }
        }
        printf("%zu / %zu training pixels correctly classified\n", ncorrect, nrow_train);

        /* Now apply it to all pixels */
        class = prf_forest_classify_table(F, features_cm, anno->M*anno->N,
                                          nfeatures);


        float * result = malloc(anno->M*anno->N*sizeof(float));
        assert(result != NULL);
        for(size_t kk = 0; kk < anno->M*anno->N; kk++)
        {
            if(class[kk] == 2)
            {
                result[kk] = 0;
            } else {
                result[kk] = 1;
            }
        }

        if(s->train_loop)
        {
            printf("Writing to training_classification.tif\n");
            fim_tiff_write_noscale("training_classification.tif", result, NULL,
                                   anno->M, anno->N, 1);
            printf("Please inspect the image.\n");
            printf("Satisfied? type y<enter>\n");
            printf("Else, modify %s and press <enter> to try again\n",
                   s->anno_label);
            int ret = fgetc(stdin);
            if(ret == 'y')
            {
                done = 1;
            } else {
                printf("Trying again\n");
            }
        } else {
            done = 1;
        }
        free(features_cma);
        free(features_cma_train);
        fimo_free(anno);
        free(result);
    }

    return F;
}

void random_forest_pipeline(opts * s, int argc, char ** argv)
{

    /* Read raw image */
    fimo * anno_raw = get_reduction(s, s->anno_image);

    /* Extract features */
    ftab_t * features = fim_features_2d(anno_raw);

    /* Transpose to column-major */
    float * features_cm = transpose(features->T,
                                    features->nrow, features->ncol);

    PrfForest * F = loop_training_data(s,
                                       features_cm,
                                       features->nrow, features->ncol);

    for(size_t kk = s->optpos; kk < (size_t) argc; kk++)
    {
        segment_image_rf(s, F, argv[kk]);
    }

    prf_forest_free(F);

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

    #ifdef _OPENMP
    omp_set_num_threads(s->nthreads);
#endif

    if(s->anno_image == NULL)
    {
        printf("No --aimage specified, can't continue\n");
        goto done;
    }

    if(s->anno_label == NULL)
    {
        printf("No --alabel image specified.\n");
        printf("Please create one and run again.\n");
        goto done;
    }

    if(s->verbose > 1)
    {
        opts_print(stdout, s);
    }

    random_forest_pipeline(s, argc, argv);

 done:
    opts_free(s);
    return 0;
}
