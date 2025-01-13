#include "dw_nuclei.h"

#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dw_png.h"
#include "dw_util.h"
#include "dw_version.h"
#include "fim.h"
#include "fim_tiff.h"
#include "trafo/src/trafo.h"

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

    char * bg; // file name of background model image
    fimo * bg_model;
} opts;


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
    fimo_free(s->bg_model);
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
    if(s->purpose == NUC_FIT)
    {
        fprintf(f, "Number of trees: %d\n", s->ntree);
    }
    return;
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
    free(s->bg);
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

/* Read png image and decode colors to class IDs
 * unset -> 0
 * red -> 1
 * green -> 2
 * blue -> 3
 *
 * For reference:
 * http://www.libpng.org/pub/png/libpng-manual.txt
 */
fimo * fim_png_decode(const char * fname)
{
    u32 height, width;
    u8 * png_data = rgb_from_png(fname, &width, &height);
    if(png_data == NULL)
    {
        fprintf(stderr, "Error reading %s\n", fname); fflush(stdout);
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
        F->V[kk] = -1;
        u8 r = png_data[3*kk+0];
        u8 g = png_data[3*kk+1];
        u8 b = png_data[3*kk+2];
        if(r > g && r > b)
        {
            F->V[kk] = 0;
            nfg++;
        }
        if(g > r && g > b)
        {
            F->V[kk] = 1;
            nbg++;
        }
        if(b > r && b > g)
        {
            F->V[kk] = 2;
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
    printf("To create png images (.a.png) for drawing annotations:\n");
    printf("usage: %s [<options>] --init file1.tif ... \n", argv[0]);
    printf("Create a classifier based on annotated images:\n");
    printf("usage: %s [<options>] --fit my_model.trf file1.tif ... \n", argv[0]);
    printf("Predict/classify images:\n");
    printf("usage: %s [<options>] --predict my_model.trf file1.tif ... \n", argv[0]);
    printf("\n");
    printf("Options:\n");
    printf(" --init\n\t"
           "create png images for drawing annotations\n");
    printf(" --fit model.trf\n\t"
           "Fit a model to the supplied training images\n");
    printf(" --predict model.trf\n\t"
           "Predict/classify images");
    printf(" --rmax\n\t"
           "use max projection over z as 3D->2D reduction\n");
    printf(" --rfocus\n\t"
           "use the plane most in focus as 3D->2D reduction\n");
    printf(" --rmean\n\t"
           "use the mean value over z as 3D->2D reduction\n");
    printf(" --bg\n\t"
           "specify a background model used for vignetting correction\n");
    printf(" --overwrite\n\t Overwrite existing files\n");
    printf(" --help\n\t Show this message\n");
}

static void argparsing(int argc, char ** argv, opts * s)
{
    struct option longopts[] = {
        {"init", no_argument, NULL, 'I'},
        {"fit", required_argument, NULL, 'F'},
        {"predict", required_argument, NULL, 'P'},

        {"help", no_argument, NULL, 'h'},
        {"overwrite", no_argument, NULL, 'o'},

        {"threads", required_argument, NULL, 't'},
        {"verbose", required_argument, NULL, 'v'},
        {"rmax", no_argument, NULL, '1'},
        {"rfocus", no_argument, NULL, '2'},
        {"rmean", no_argument, NULL, '3'},
        {"bg", required_argument, NULL, 'b'},
        {NULL, 0, NULL, 0}};
    int ch;
    while((ch = getopt_long(argc, argv, "IF:P:hot:v:123b", longopts, NULL)) != -1)
    {
        switch(ch){
        case 'F':
            free(s->modelfile);
            s->modelfile = strdup(optarg);
            s->purpose = NUC_FIT;
            break;
        case 'I':
            s->purpose = NUC_INIT;
            break;
        case 'P':
            s->purpose = NUC_CLASSIFY;
            free(s->modelfile);
            s->modelfile = strdup(optarg);
            break;
        case 'o':
            s->overwrite = 1;
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
            break;
        case 't':
            s->nthreads = atoi(optarg);
            break;
        case 'v':
            s->verbose = atoi(optarg);
            break;
        case '1':
            s->redu = REDU_MAX;
            break;
        case '2':
            s->redu = REDU_FOCUS;
            break;
        case '3':
            s->redu = REDU_MEAN;
            break;
        case 'b':
            free(s->bg);
            s->bg = strdup(optarg);
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

/* A reduction is a mapping from 3D to 2D */
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
        printf("Scaling image by %f\n", 1.0/scaling);
    }
    fim_mult_scalar(I, M*N*P, 1.0/scaling);


    if(s->redu == REDU_MAX)
    {
        float * maxI = fim_maxproj(I, M, N, P);
        free(I);
        fimo * result = fim_image_from_array(maxI, M, N, 1);
        free(maxI);
        if(s->bg_model != NULL)
        {
            printf("Applying background model\n");
            if(fimo_div_image(result, s->bg_model))
            {
                printf("Background correction failed\n");
            }
        }

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
        if(s->bg_model != NULL)
        {
            fimo_div_image(result, s->bg_model);
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
        case 0:
            pixels[3*kk+0] = 255;
            break;
        case 1:
            pixels[3*kk+1] = 255;
            break;
        case 2:
            pixels[3*kk+2] = 255;
            break;
        default:
            pixels[3*kk+0] = 255;
            pixels[3*kk+1] = 255;
            pixels[3*kk+2] = 255;
            break;
            ;
        }
    }
    if(rgb_to_png(pixels, M, N, name_out))
    {
        fprintf(stderr, "Failed to write to %s\n", name_out);
    }
    free(pixels);
    return;
}


static ftab_t *
read_png_labels(const char * name_image, u32 * height, u32 *width)
{
    fimo * anno = fim_png_decode(name_image);
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

    char * logfile = calloc(strlen(s->modelfile) + 16, 1);
    assert(logfile != NULL);
    sprintf(logfile, "%s.log.txt", s->modelfile);
    FILE * log = fopen(logfile, "w");
    for(int kk = 0; kk < argc; kk++)
    {
        fprintf(log, "%s ", argv[kk]);
    }
    fprintf(log, "\n");
    fim_tiff_set_log(log);

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
        if(image->M != width && image->N != height)
        {
            printf("Error: PNG and TIF image dimensions mismatch\n");
            exit(EXIT_FAILURE);
        }
        /* Extract features */
        ftab_t * image_features = fim_features_2d(image);
        // TODO: Collect features names (or get the in some other way)
        fimo_free(image);

        assert(image_features->nrow == image_labels->nrow);

        /* Subselect rows where the class != 0 from both tables */
        u8 * labeled_pixels = calloc(image_features->nrow, sizeof(u8));
        for(size_t kk = 0; kk < image_features->nrow; kk++)
        {
            image_labels->T[kk] > -1 ? labeled_pixels[kk] = 1 : 0;
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
        // TODO: Print feature names
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

    if(s->verbose > 0)
    {
        double * fimp = trafo_importance(F);
        printf("Feature importance\n");
        for(u32 kk = 0; kk < nfeature; kk++)
        {
            printf("%d : %f\n", kk+1, fimp[kk]);
        }
        free(fimp);
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

    fim_tiff_set_log(stdout);
    fclose(log);
    return EXIT_SUCCESS;
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


static int
init(opts * s, int argc, char ** argv)
{
    if(s->optpos >= argc)
    {
        printf("No files given\n");
        return EXIT_FAILURE;
    }
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
        // printf("%lu x %lu\n", M, N);

        fimo * Ir = get_reduction(s, image);
        //fimo * I = fimo_tiff_read(image);
        //fimo * Iz = fimo_maxproj(I);

        //fimo_free(I);
        fimo_to_png(Ir, out);
        fimo_free(Ir);
        free(out);
    }
    return EXIT_SUCCESS;
}

static int
load_background(opts * s)
{

    s->bg_model = fimo_tiff_read(s->bg);
    if(s->bg_model == NULL)
    {
        printf("ERROR: Unable to load a background model from %s\n", s->bg);
        return 1;
    }
    if(s->bg_model->P != 1)
    {
        printf("ERROR: The background image is not 2D\n");
        return 1;
    }
    float max = fimo_max(s->bg_model);
    if(max != 1.0)
    {
        fimo_mult_scalar(s->bg_model, 1.0/max);
    }
    return 0;
}

int
dw_nuclei(int argc, char ** argv)
{

    fim_tiff_init();
    FILE * out = fopen("/dev/null", "w");
    fim_tiff_set_log(out);
    opts * s = opts_new();
    int status = EXIT_SUCCESS;

    argparsing(argc, argv, s);

#ifdef _OPENMP
    omp_set_num_threads(s->nthreads);
#endif

    if(s->verbose > 0)
    {
        opts_print(stdout, s);
    }

    if(s->bg != NULL)
    {
        if(load_background(s))
        {
            printf("Unable to use file %s for background correction\n", s->bg);
            opts_free(s);
            return EXIT_FAILURE;
        }
    }

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

    opts_free(s);
    return status;
}

#ifdef STANDALONE
int main(int argc, char ** argv)
{
    return dw_nuclei(argc, argv);
}
#endif
