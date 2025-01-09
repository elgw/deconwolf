#include "dw_nuclei.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <stdint.h>


#include <omp.h>

#include "fim.h"
#include "fim_tiff.h"
#include "dw_util.h"
#include "dw_version.h"
#include "trafo/src/trafo.h"
#include "dw_png.h"

/* Suggested command line interface:
 * dw_nuclei --fit model.name [more options] file1.tif annotations1. tif file2.tif annotations2.tif ...
 * dw_nuclei --classify model.name [more options] file1.tif file2.tif ...
 *
 */

typedef uint32_t u32;
typedef uint64_t u64;
typedef int64_t i64;
typedef uint8_t u8;

typedef enum {
    REDU_MAX,
    REDU_MEAN,
    REDU_FOCUS
} reduction_type;

typedef enum {
    /* Create 2D PNG images for a set of input images */
    NUC_INIT,
    NUC_UNSET,
    NUC_FIT,
    NUC_CLASSIFY,
    /* Test a classifier on annotated data (like fit but load a
       model) */
    NUC_EVALUATE,
} purpose_type;

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

    /* File to read/write model to  */
    purpose_type purpose;
    char * modelfile;
    reduction_type redu;
} opts;

static opts * opts_new();
static void opts_free(opts * s);
static void opts_print(FILE * f, opts * s);
static void usage(__attribute__((unused)) int argc, char ** argv);
static void argparsing(int argc, char ** argv, opts * s);

static opts * opts_new()
{
    opts * s = calloc(1, sizeof(opts));
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
    s->purpose = NUC_UNSET;
    return s;
}

static void opts_print(FILE * f, opts * s)
{
    fprintf(f, "Overwrite: %d\n", s->overwrite);
    fprintf(f, "3D->2D reduction: ");
    switch(s->redu)
    {
    case REDU_MAX:
        fprintf(f, "Max projection\n");
        break;
    case REDU_FOCUS:
        fprintf(f, "Most in focus slice\n");
        break;
    case REDU_MEAN:
        fprintf(f, "Mean projection\n");
        break;
    }
    fprintf(f, "Number of trees: %d\n", s->ntree);
}


static void opts_free(opts * s)
{
    if(s == NULL)
    {
        return;
    }
    free(s->out);
    free(s->image);
    free(s->anno_image);
    free(s->anno_label);
    free(s->modelfile);
    free(s);
}

int fimo_to_png(fimo * I, const char * outname)
{
    if(I->P != 1)
    {
        fprintf(stderr, "fimo_to_png: ERROR: can only write 2D images\n");
        return EXIT_FAILURE;
    }

    uint8_t * img_data = calloc(3*I->M*I->N*3, sizeof(uint8_t));
    assert(img_data != NULL);

    float imax = fim_max(I->V, I->M*I->N);
    for(u64 kk = 0; kk < I->M*I->N; kk++)
    {
        float v = I->V[kk]/imax*255.0;
        u8 v8 = (u8) round(v);
        if(v > 255)
        {
            v8 = 255;
        }
        img_data[3*kk] = v8; // Red
        img_data[3*kk+1] = v8; // Green
        img_data[3*kk+2] = v8; // Blue
    }


    int status = rgb_to_png(img_data, I->M, I->N, outname);
    free(img_data);

    return status;
}

double * double_from_float(const float * source, size_t n)
{
    assert(source != NULL);
    double * target = calloc(n, sizeof(double));
    assert(target != NULL);
    for(size_t kk = 0; kk < n; kk++)
    {
        target[kk] = source[kk];
    }
    return target;
}

/* Read png image, when more green that red, set to 1, when more red
 * than green set to 2. Everything else set to 0.
 * http://www.libpng.org/pub/png/libpng-manual.txt
 */
fimo * fim_png_read_green_red(const char * fname)
{
    printf("Reading %s\n", fname); fflush(stdout);

    u32 height, width;
    u8 * png_data = rgb_from_png(fname, &width, &height);
    if(png_data == NULL)
    {
        return NULL;
    }

    fimo * F = NULL;
    F = calloc(1, sizeof(fimo));
    assert(F != NULL);
    F->M = height;
    F->N = width;
    F->P = 1;
    F->V = calloc(F->M*F->N, sizeof(float));
    assert(F->V != NULL);

    size_t nfg = 0;
    size_t nbg = 0;
    for(size_t kk = 0; kk<F->M*F->N; kk++)
    {
        F->V[kk] = 0;
        if(png_data[3*kk+1] > png_data[3*kk+0])
        {
            F->V[kk] = 1;
            nfg++;
        }
        if(png_data[3*kk+0] > png_data[3*kk+1])
        {
            F->V[kk] = 2;
            nbg++;
        }
    }
    printf("Found %zu foreground (green->2) and %zu background (red->1) pixels\n",
           nfg, nbg);

    return F;
}

static void usage(__attribute__((unused)) int argc, char ** argv)
{
    printf("This module can produce a random forest classifier based on\n"
           "annotated images and classify non-annotated images\n"
           "\n");
    printf("usage: %s [<options>] --aimage input.tif --alabel labels.png file1.tif ... \n", argv[0]);
    printf("\n");
    printf("Options:\n");
    printf(" --aimage file.tif\n\t input image\n");
    printf(" --alabel file.png\n\t "
           "Corresponding annotation where nuclei is marked green, background is marked red\n");

    printf("--init\n\t"
           "png images for drawing annotations\n");
    printf(" --fit model.trf\n\t"
           "Fit a model to the supplied training images\n");
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
        {"fit", required_argument, NULL, 'F'},
        {"help", no_argument, NULL, 'h'},
        {"image", required_argument, NULL, 'i'},
        {"init", no_argument, NULL, 'I'},
        {"loop", no_argument, NULL, 'l'},
        {"overwrite", no_argument, NULL, 'o'},
        {"out",     required_argument, NULL, 'p'},
        {"predict", required_argument, NULL, 'P'},
        {"threads", required_argument, NULL, 't'},
        {"verbose", required_argument, NULL, 'v'},
        {NULL, 0, NULL, 0}};
    int ch;
    while((ch = getopt_long(argc, argv, "a:b:F:hi:Iop:P:r:v:", longopts, NULL)) != -1)
    {
        switch(ch){
        case 'a':
            free(s->anno_label);
            s->anno_label = strdup(optarg);
            break;
        case 'b':
            free(s->anno_image);
            s->anno_image = strdup(optarg);
            break;
        case 'f':
            free(s->image);
            s->image = strdup(optarg);
            break;
        case 'F':
            free(s->modelfile);
            s->modelfile = strdup(optarg);
            s->purpose = NUC_FIT;
            break;
        case 'I':
            s->purpose = NUC_INIT;
            break;
        case 'l':
            s->train_loop = 1;
            break;
        case 'o':
            s->overwrite = 1;
            break;
        case 'p':
            free(s->out);
            s->out = strdup(optarg);
            break;
        case 'P':
            s->purpose = NUC_CLASSIFY;
            free(s->modelfile);
            s->modelfile = strdup(optarg);
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
            break;
        case 'i':
            free(s->image);
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


/* Transpose a dense table/image */
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
static float *
cm_append_column(float * A, size_t nrow, size_t ncol, const float * B)
{
    float * AA = malloc(nrow*(ncol+1)*sizeof(float));
    assert(AA != NULL);
    memcpy(AA, A,
           nrow*ncol*sizeof(float));
    memcpy(AA+nrow*ncol, B, nrow*sizeof(float));
    return AA;
}

/* A: Input table data
 * V:  Binary indicator if elements should be included
 */

static float *
subset_cm(const float * A,
          const float * V,
          const size_t nrow, const size_t ncol,
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

static int
float_arg_max(const float * v, size_t N)
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

static fimo *
get_reduction(opts * s, const char * file)
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
        fimo * result = calloc(1, sizeof(fimo));
        assert(result != NULL);
        result->M = M;
        result->N = N;
        result->P = 1;
        result->V = calloc(M*N, sizeof(float));
        assert(result->V != NULL);
        memcpy(result->V, II->V+slice*M*N, M*N*sizeof(float));
        fimo_free(II);
        if(s->verbose > 0)
        {
            printf("Returning slice %d\n", slice); fflush(stdout);
        }
        free(gm);
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

void
segment_image_rf(opts * s, trf * F, char * file)
{
    if(!dw_isfile(file))
    {
        printf("%s does not exist\n", file);
        return;
    }
    char * outfile = calloc(strlen(file) + 32, 1);
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

    double * features_cm_double = double_from_float(features_cm, features->nrow*features->ncol);
    free(features_cm);

    printf("3\n"); fflush(stdout);
    /* Classify */
    u32 * class = trafo_predict(F, features_cm_double, NULL, redu->M);

    ftab_free(features);
    free(features_cm_double);
    printf("4\n"); fflush(stdout);
    float * result = calloc(redu->M*redu->N, sizeof(float));
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

/* The loop is a feedback loop with the user  */
trf *
loop_training_data(opts * s, const float * features_cm,
                   const size_t n_sample, const size_t n_feature)
{

    trf * F = NULL;
    int done = 0;
    while(done == 0)
    {
        if(F != NULL)
        {
            trafo_free(F);
        }

        /* Read annotated image -> labels*/
        fimo * anno = fim_png_read_green_red(s->anno_label);
        assert(anno != NULL);
        size_t n_pixel = anno->M*anno->N;
        /* Select only annotated positions and convert to u32 */
        u32 * selected_labels = calloc(n_pixel, sizeof(u32));
        assert(selected_labels != NULL);

        size_t n_label = 0;
        for(size_t kk = 0; kk < anno->M*anno->N; kk++)
        {
            if(anno->V[kk] > 0)
            {
                selected_labels[n_label] = anno->V[kk];
            }
        }

        fimo_free(anno);


        /* Extract the data where anno->V > 0 to a separate array */
        size_t ncol_train = 0;
        size_t nrow_train = 0;
        float * features_cm_train = subset_cm(features_cm,
                                              anno->V,
                                              n_sample, n_feature+1,
                                              &nrow_train, &ncol_train);

        assert(n_label == nrow_train);

        double * features_cm_train_double = double_from_float(features_cm_train,
                                                              nrow_train*ncol_train);
        free(features_cm_train);

        /* Train classifier */
        trafo_settings tconf = {0};
        printf("TODO: THIS IS WRONG\n"); fflush(stdout);
        tconf.label = selected_labels;
        tconf.F_col_major = features_cm_train_double;
        tconf.n_feature = ncol_train;
        tconf.n_sample = nrow_train;
        tconf.n_tree = 200;

        F = trafo_fit(&tconf);

        if(F == NULL)
        {
            printf("dw_nuclei: Failed to train the random forest\n");
            exit(EXIT_FAILURE);
        }

        printf("Validating the training data\n");
        u32 * class = trafo_predict(F, features_cm_train_double, NULL, nrow_train);

        size_t ncorrect = 0;
        for(size_t kk = 0; kk<nrow_train; kk++)
        {
            if(class[kk] == selected_labels[kk])
            {
                ncorrect++;
            }
        }
        printf("%zu / %zu training pixels correctly classified\n",
               ncorrect, nrow_train);
        free(class);

        /* Now apply it to all pixels */

        double * features_cm_double =
            double_from_float(features_cm, n_feature*n_sample);

        class = trafo_predict(F, features_cm_double, NULL, anno->M);

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

        fimo_free(anno);
        free(result);
    }

    return F;
}

void
random_forest_pipeline(opts * s, int argc, char ** argv)
{

    /* Read raw image */
    fimo * anno_raw = get_reduction(s, s->anno_image);

    /* Extract features */
    ftab_t * features = fim_features_2d(anno_raw);

    /* Transpose to column-major */
    float * features_cm = transpose(features->T,
                                    features->nrow, features->ncol);

    trf * F = loop_training_data(s,
                                 features_cm,
                                 features->nrow, features->ncol);

    for(size_t kk = s->optpos; kk < (size_t) argc; kk++)
    {
        segment_image_rf(s, F, argv[kk]);
    }

    trafo_free(F);

    return;
}

static ftab_t *
read_png_labels(const char * name_image, u32 * height, u32 *width)
{
    fimo * anno = fim_png_read_green_red(name_image);
    *height = anno->M;
    *width = anno->N;
    assert(anno->P == 1);
    if(anno == NULL)
    {
        printf("Unable to read %s\n", name_image);
        return NULL;
    }

    ftab_t * tab_anno = ftab_new_from_data(anno->M*anno->N, 1, anno->V);
    fimo_free(anno);
    ftab_set_colname(tab_anno, 0, "class");
    return tab_anno;
}

/* Read annotations from one or several images and fit a model */
static int
fit(opts * s, int argc, char ** argv)
{
    if(s->overwrite == 0)
    {
        if(dw_isfile(s->modelfile))
        {
            printf("%s does already exist (--overwrite not specified)\n",
                   s->modelfile);
            return EXIT_FAILURE;
        }
    }
    int nimages = (argc - s->optpos);
    if( nimages == 0)
    {
        fprintf(stderr, "Error: No images given\n");
        return EXIT_FAILURE;
    }

    ftab_t * fit_features = NULL;
    ftab_t * fit_labels = NULL;

    printf("Loading data from %d image pairs\n", nimages);
    for(int kk = s->optpos; kk < argc; kk++)
    {
        const char * name_image = argv[kk];
        char * name_annotation = calloc(strlen(name_image)+16, sizeof(char));
        assert(name_annotation != NULL);
        sprintf(name_annotation, "%s.a.png", name_image);
        if(!dw_isfile(name_annotation))
        {
            printf("%s does not exist\n", name_annotation);
            free(name_annotation);
            continue;
        }
        printf("Processing %s (%s)\n", name_image, name_annotation);

        u32 height, width;
        ftab_t * image_labels = read_png_labels(name_annotation, &height, &width);
        free(name_annotation);
        if(image_labels == NULL)
        {
            continue;
        }

        /* Read raw image and reduce it to 2D*/
        fimo * image = get_reduction(s, name_image);
        if(image->M != height)
        {
            printf("Error: PNG and TIF image dimensions mismatch\n");
            exit(EXIT_FAILURE);


        }
        /* Extract features */
        ftab_t * image_features = fim_features_2d(image);
        fimo_free(image);

        assert(image_features->nrow == image_labels->nrow);

        /* Subselect rows where the class != 0 from both tables */
        u8 * labeled_pixels = calloc(image_features->nrow, sizeof(u8));
        for(size_t kk = 0; kk < image_features->nrow; kk++)
        {
            image_labels->T[kk] > 0 ? labeled_pixels[kk] = 1 : 0;
        }
        ftab_subselect_rows(image_labels, labeled_pixels);
        ftab_subselect_rows(image_features, labeled_pixels);
        free(labeled_pixels);

        /* Concatenate tables */
        ftab_t * fit_features2 = ftab_concatenate_rows(fit_features,
                                                       image_features);
        ftab_free(fit_features);
        ftab_free(image_features);
        fit_features = fit_features2;

        ftab_t * fit_labels2 = ftab_concatenate_rows(fit_labels,
                                                     image_labels);
        ftab_free(fit_labels);
        ftab_free(image_labels);
        fit_labels = fit_labels2;
    }

    /* Convert features to double and class labels to XXX */
    u32 * fit_labels4 = calloc(fit_labels->nrow, sizeof(u32));
    assert(fit_labels4 != NULL);
    double * fit_features8 = calloc(fit_features->nrow*fit_features->ncol,
                                    sizeof(double));
    assert(fit_features8 != NULL);
    for(u64 kk = 0; kk < fit_labels->nrow; kk++)
    {
        fit_labels4[kk] = (u32) fit_labels->T[kk];
    }
    ftab_free(fit_labels);
    for(u64 kk = 0; kk < fit_features->nrow*fit_features->ncol; kk++)
    {
        fit_features8[kk] = fit_features->T[kk];
    }
    u64 nfeature = fit_features->ncol;
    u64 nsample = fit_features->nrow;
    ftab_free(fit_features);


    /* Train classifier */
    trafo_settings tconf = {0};
    tconf.label = fit_labels4; // selected_labels;
    tconf.F_col_major = NULL; //features_cm_train_double;
    tconf.F_row_major = fit_features8;
    tconf.n_feature = nfeature; // ncol_train;
    tconf.n_sample = nsample; // nrow_train;
    tconf.n_tree = 200;

    trf * F = trafo_fit(&tconf);



    if(F == NULL)
    {
        printf("dw_nuclei: Failed to train the random forest\n");
        exit(EXIT_FAILURE);
    }

    printf("Validating the training data\n");
    u32 * class = trafo_predict(F, NULL, fit_features8, nsample);
    free(fit_features8);
    size_t ncorrect = 0;
    for(size_t kk = 0; kk<nsample; kk++)
    {
        if(class[kk] == fit_labels4[kk])
        {
            ncorrect++;
        }
    }
    printf("%zu / %zu training pixels correctly classified\n",
           ncorrect, nsample);
    free(fit_labels4);
    free(class);

    /* Save classifier */
    trafo_save(F, s->modelfile);
    trafo_free(F);

    return EXIT_SUCCESS;
}

/* 1-> Green
*  2->Red
*/
static void
labels_to_png(u32 * label, u32 M, u32 N, const char * name_out)
{
    u8 * pixels = calloc(M*N*3, sizeof(u8));
    for(u64 kk = 0; kk < M*N; kk++)
    {
        switch(label[kk])
        {
        case 1:
            pixels[3*kk+1] = 255;
            break;
        case 2:
            pixels[3*kk+0] = 255;
            break;
        default:
            break;
            ;
        }
    }
    if(rgb_to_png(pixels, N, M, name_out))
    {
        fprintf(stderr, "Failed to write to %s\n", name_out);
    }
    free(pixels);
    return;
}

static int
classify(opts * s, int argc, char ** argv)
{
    // Load classifier/model
    trf * RF = trafo_load(s->modelfile);
    if(RF == NULL)
    {
        fprintf(stderr, "Unable to load a random forest model from %s\n",
                s->modelfile);
        return EXIT_FAILURE;
    }

    for(int kk = s->optpos; kk < argc; kk++)
    {
        const char * name_image = argv[kk];
        char * name_out = calloc(strlen(name_image)+16, 1);
        sprintf(name_out, "%s.predict.png", name_image);
        if(s->overwrite == 0 && dw_isfile(name_out))
        {
            printf("%s does already exist (--overwrite not used)\n", name_out);
            free(name_out);
            continue;
        }
        /* Read raw image and reduce it to 2D*/
        fimo * image = get_reduction(s, name_image);
        if(image == NULL)
        {
            continue;
        }
        /* Extract features */
        ftab_t * image_features = fim_features_2d(image);
        u32 M = image->M;
        u32 N = image->N;
        fimo_free(image);

        double * image_features8 = calloc(image_features->nrow*image_features->ncol,
                                          sizeof(double));
        assert(image_features8 != NULL);
        for(u64 kk = 0; kk < image_features->nrow*image_features->ncol; kk++)
        {
            image_features8[kk] = image_features->T[kk];
        }
        u64 nsample = image_features->nrow;
        assert(nsample == M*N);
        ftab_free(image_features);
        u32 * class = trafo_predict(RF, NULL, image_features8, nsample);
        free(image_features8);

        // Write to disk ...
        labels_to_png(class, M, N, name_out);

        free(name_out);
        free(class);
    }

    trafo_free(RF);
    return EXIT_SUCCESS;
}

#ifdef STANDALONE
int main(int argc, char ** argv)
{
    return dw_nuclei(argc, argv);
}
#endif

static int
init(opts * s, int argc, char ** argv)
{
    // For each image: Load, create max projection, save as PNG
    for(int kk = s->optpos; kk < argc; kk++)
    {
        const char * image = argv[kk];
        char * out = calloc(strlen(image)+16, 1);

        assert(out != NULL);
        sprintf(out, "%s.a.png", image);
        printf("%s -> %s\n", image, out);
        if(s->overwrite == 0 && dw_isfile(out) == 1)
        {
            printf("%s does already exist (--overwrite not specified)\n", out);
            free(out);
            continue;
        }
        i64 M, N, P;
        fim_tiff_get_size(image, &M, &N, &P);
        printf("%lu x %lu\n", M, N);
        fimo * I = fimo_tiff_read(image);
        fimo * Iz = fimo_maxproj(I);
        fimo_free(I);
        fimo_to_png(Iz, out);
        fimo_free(Iz);
        free(out);
    }
    return EXIT_SUCCESS;
}


int
dw_nuclei(int argc, char ** argv)
{

    fim_tiff_init();
    opts * s = opts_new();
    int status = EXIT_SUCCESS;

    argparsing(argc, argv, s);

#ifdef _OPENMP
    omp_set_num_threads(s->nthreads);
#endif

    switch(s->purpose)
    {
    case NUC_INIT:
        status = init(s, argc, argv);
        break;
    case NUC_FIT:
        status = fit(s, argc, argv);
        break;
    case NUC_CLASSIFY:
        status = classify(s, argc, argv);
        break;
    default:
        ;
    }

    /* Old stuff to be removed */
    if(s->purpose == NUC_UNSET)
    {
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
    }
 done:
    opts_free(s);
    return status;
}
