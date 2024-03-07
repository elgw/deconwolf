#include "dw_dots.h"

typedef struct{
    int overwrite;
    int verbose;
    int optpos;
    int LoG;
    char * image; /* Image to analyze */
    char * outfile; /* Where to write the tsv output */
    char * fout; /* Where to (optionally) write the filtered image */
    int nthreads;
    float asigma; /* Axial sigma */
    float lsigma; /* Lateral sigma */
    int ndots; /* Number of dots to export */
    char * logfile;
    FILE * log;
    int fwhm;
    float th;
    int optind;
} opts;

static opts * opts_new();
static void opts_free(opts * s);
static void usage(__attribute__((unused)) int argc, char ** argv);
static void argparsing(int argc, char ** argv, opts * s);

static opts * opts_new()
{
    opts * s = calloc(1, sizeof(opts));
    assert(s != NULL);
    s->verbose = 1;
    s->optpos = -1;
    s->nthreads = dw_get_threads();
    s->ndots = -1; /* Auto */
    s->log = NULL;
    s->logfile = NULL;
    s->fout = NULL;
    s->fwhm = 1;
    return s;
}


static void opts_free(opts * s)
{
    free(s->outfile);
    free(s->image);
    free(s->logfile);
    if(s->log != NULL)
    {
        fclose(s->log);
    }
    free(s);
}

static void opts_print(FILE * f, opts * s)
{
    assert(s != NULL);
    assert(f != NULL);
    fprintf(f, "overwrite = %d\n", s->overwrite);
    fprintf(f, "verbose = %d\n", s->verbose);
    if(s->image != NULL)
    {
        fprintf(f, "image: %s\n", s->image);
    }
    if(s->outfile != NULL)
    {
        fprintf(f, "outfile: %s\n", s->outfile);
    }
    if(s->logfile != NULL)
    {
        fprintf(f, "logfile: %s\n", s->logfile);
    }
    fprintf(f, "nthreads = %d\n", s->nthreads);
    fprintf(f, "Lateral sigma: %.2f\n", s->lsigma);
    fprintf(f, "Axial sigma: %.2f\n", s->asigma);
    fprintf(f, "Laplacian of Gaussian: %d\n", s->LoG);
    if(s->ndots > 0)
    {
        fprintf(f, "Will write at most %d dots\n", s->ndots);
    } else {
        fprintf(f, "Automatic number of dots in the output\n");
    }
    if(s->fwhm)
    {
        fprintf(f, "Will calculate FWHM on the filtered image\n");
    } else {
        fprintf(f, "No FWHM calculations\n");
    }
    return;
}

static void usage(__attribute__((unused)) int argc, char ** argv)
{
    opts * s = opts_new();
    printf("usage: %s [<options>] input.tif input2.tif ...\n", argv[0]);
    printf("Options:\n");
    printf(" --LoG\n\t Use Laplacian of Gaussian (recommended for "
           "non-deconvolved images, default=%d).\n", s->LoG);
    printf(" --lsigma s\n\t Lateral sigma (default %.1f)\n", s->lsigma);
    printf(" --asigma s\n\t Axial sigma (default %.1f)\n", s->asigma);
    printf(" --overwrite\n\t Overwrite existing files (default %d)\n", s->overwrite);
    printf(" --help\n\t Show this message\n");
    printf(" --ndots n\n\t Number of dots to export (default %d)\n", s->ndots);
    printf(" --logfile file.txt\n\t Specify where the log file should be written\n");
    printf(" --fwhm\n\t Include FWHM in the output (default %d). Will be "
           "based on the filtered image.\n", s->fwhm);
    printf(" --verbose v\n\t Verbosity level (default %d)\n", s->verbose);
    printf(" --nthreads n\n\t Set the number of computational threads\n");
    printf(" --fout file.tif\n\t Write filtered image -- for debugging\n");
    printf("\n");
    printf("The log messages will be written to input.tif.log.txt\n");
    printf("Dots will be exported to input.tif.dots.tsv\n");
    free(s);
}

ftab_t * ftab_insert_col(ftab_t * T, float * C, const char * cname)
{
    /* Create a new table with one extra column */
    ftab_t * T2 = ftab_new(T->ncol + 1);
    for(size_t cc = 0; cc<T->ncol; cc++)
    {
        ftab_set_colname(T2, cc, T->colnames[cc]);
    }
    float * row = malloc((T->nrow+1)*sizeof(float));
    assert(row != NULL);
    for(size_t rr = 0; rr<T->nrow; rr++)
    {
        memcpy(row, T->T+rr*T->ncol, T->ncol*sizeof(float));
        row[T->ncol] = C[rr];
        ftab_insert(T2, row);
    }
    free(row);
    ftab_set_colname(T2, T->ncol, cname);
    ftab_free(T);
    return T2;
}

static void argparsing(int argc, char ** argv, opts * s)
{
    struct option longopts[] = {
        {"LoG", no_argument, NULL, '2'},
        {"logfile", required_argument, NULL, 'L'},
        {"out", required_argument, NULL, 'O'},
        {"asigma", required_argument, NULL, 'a'},
        {"fwhm", no_argument, NULL, 'f'},
        {"help", no_argument, NULL, 'h'},
        {"lsigma", required_argument, NULL, 'l'},
        {"ndots",   required_argument, NULL, 'n'},
        {"overwrite", no_argument, NULL, 'o'},
        {"fout",     required_argument, NULL, 'p'},
        {"threads", required_argument, NULL, 't'},
        {"verbose", required_argument, NULL, 'v'},
        {NULL, 0, NULL, 0}};
    int ch;
    while((ch = getopt_long(argc, argv, "2L:a:hi:l:n:op:r:v:", longopts, NULL)) != -1)
    {
        switch(ch){
        case '2':
            s->LoG = 1;
            break;
        case 'L':
            free(s->logfile);
            s->logfile = strdup(optarg);
            assert(s->logfile != NULL);
            break;
        case 'O':
            free(s->outfile);
            s->outfile = strdup(optarg);
            assert(s->outfile != NULL);
            break;
        case 'a':
            s->asigma = atof(optarg);
            break;
        case 'f':
            s->fwhm = 1;
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
            break;
        case 'l':
            s->lsigma = atof(optarg);
            break;
        case 'n':
            s->ndots = atoi(optarg);
            break;
        case 'o':
            s->overwrite = 1;
            break;
        case 'p':
            s->fout = strdup(optarg);
            assert(s->fout != NULL);
            break;
        case 't':
            s->nthreads = atoi(optarg);
            break;
        case 'v':
            s->verbose = atoi(optarg);
            break;
        default:
            exit(EXIT_FAILURE);
        }
    }


    s->optpos = optind;
    s->optpos = optind;
    return;
}

static ftab_t * append_fwhm(opts * s, ftab_t * T, fim_t * V,
                            const char * cname_lateral,
                            const char * cname_axial)
{
    int xcol = ftab_get_col(T, "x");
    int ycol = ftab_get_col(T, "y");
    int zcol = ftab_get_col(T, "z");

    struct timespec tstart, tend;


    if(s->verbose > 1)
    {
        printf("Calculating fwhm..."); fflush(stdout);
    }
    dw_gettime(&tstart);
    float * fwhm_vals_lateral = malloc(T->nrow*sizeof(float));
    assert(fwhm_vals_lateral != NULL);
    float * fwhm_vals_axial = malloc(T->nrow*sizeof(float));
    assert(fwhm_vals_axial != NULL);

#pragma omp parallel for
    for(size_t kk = 0; kk<T->nrow; kk++)
    {
        float * row = T->T + kk*T->ncol;
        fwhm_vals_lateral[kk] = fwhm_lateral(V, row[xcol], row[ycol], row[zcol], s->verbose);
        fwhm_vals_axial[kk] = fwhm_axial(V, row[xcol], row[ycol], row[zcol], s->verbose);
    }
    dw_gettime(&tend);
    double t = timespec_diff(&tend, &tstart);
    if(s->verbose > 1)
    {
        printf(" took %.4f s\n", t);
    }

    T = ftab_insert_col(T, fwhm_vals_lateral, cname_lateral);
    free(fwhm_vals_lateral);
    T = ftab_insert_col(T, fwhm_vals_axial, cname_axial);
    free(fwhm_vals_axial);

    return T;
}


#ifdef STANDALONE
int main(int argc, char ** argv)
{
    return dw_dots(argc, argv);
}
#endif

void detect_dots(opts * s, char * inFile)
{
    if(!dw_file_exist(inFile))
    {
        fprintf(stderr, "Can't open %s!\n", inFile);
        exit(1);
    }

    assert(inFile != NULL);
    free(s->image);
    s->image = strdup(inFile);
    assert(s->image != NULL);

    /* Set up the log file for this image */
    free(s->logfile);
    assert(s != NULL);
    assert(s->image != NULL);
    s->logfile = malloc(strlen(s->image) + 64);
    assert(s->logfile != NULL);
    sprintf(s->logfile, "%s.dots.log.txt", s->image);
    s->log = fopen(s->logfile, "w");
    opts_print(s->log, s);


    /* Set up the name for the output file */
    free(s->outfile);
    s->outfile = malloc(strlen(s->image) + 64);
    assert(s->outfile != NULL);
    sprintf(s->outfile, "%s.dots.tsv", s->image);

    char * outFile = s->outfile;

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
        return;
    }

    fprintf(s->log, "Reading %s\n", inFile);
    int64_t M = 0, N = 0, P = 0;
    float * A = fim_tiff_read(inFile, NULL, &M, &N, &P, s->verbose);
    float scaling = dw_read_scaling(inFile);
    if(scaling != 1)
    {
        printf("Scaling by %f\n", 1.0/scaling);
        fim_mult_scalar(A, M*N*P, 1.0/scaling);
    }

    if(s->ndots < 0)
    {
        s->ndots = 0.025 * (float) M * (float) N;
        fprintf(s->log, "Setting the max number of output dots to %d\n",
                s->ndots);
    }

    float * feature = NULL;
    if(s->LoG == 1)
    {
        if(s->verbose > 1)
        {
            printf("LoG filter, lsigma=%.2f asigma=%.2f\n", s->lsigma, s->asigma);
        }
        fim_set_verbose(2);
        feature = fim_LoG_S(A, M, N, P, s->lsigma, s->asigma);
        fim_set_verbose(0);

    } else {
        if(s->verbose > 1)
        {
            printf("Low pass filtering with lsigma=%f, asigma=%f\n",
                   s->lsigma, s->asigma);
        }
        feature = malloc(M*N*P*sizeof(float));
        assert(feature != NULL);
        memcpy(feature, A, M*N*P*sizeof(float));
        fim_gsmooth_aniso(feature, M, N, P, s->lsigma, s->asigma);
    }


    if(s->fout != NULL)
    {
        if(s->verbose > 1)
        {
            printf("Writing filtered image to %s\n", s->fout);
        }
        fim_tiff_write_float(s->fout, feature, NULL, M, N, P);
    }

    /* Detect local maxima */
    if(s->verbose > 1)
    {
        printf("Detecting local maxima\n");
    }

    ftab_t * T = fim_lmax(feature, M, N, P);

    if(s->verbose > 1)
    {
        printf("Found %zu points\n", T->ncol);
        printf("Sorting the local maxima by value\n");
    }

    /* Sort with highest value first */
    ftab_sort(T, ftab_get_col(T, "value"));

    if(s->verbose > 1)
    {
        printf("Looking for a threshold\n");
    }
    /* Get a threshold suggestion */
    float * values = malloc(sizeof(float)*T->nrow);
    assert(values != NULL);
    for(size_t kk = 0; kk<T->nrow; kk++)
    {
        values[kk] = T->T[kk*T->ncol + 3];
    }
    float mean = fim_mean(values, T->nrow);
    float std = fim_std(values, T->nrow);
    if(s->verbose > 1)
    {
        printf("Extrated %zu dots: Mean=%f, std=%f\n", T->nrow, mean, std);
    }
    fprintf(s->log, "Extracted %zu dots: Mean=%f, std=%f\n", T->nrow, mean, std);

    fim_histogram_t * H = fim_histogram(values, T->nrow);
    free(values);
    //s->th = fim_histogram_otsu(H);
    fim_histogram_log(H);
    s->th = fim_histogram_otsu(H);
    if(s->verbose > 1)
    {
        printf("Suggested threshold (from %zu dots): %f\n", T->nrow, s->th);
    }
    fprintf(s->log, "Suggested threshold (from %zu dots): %f\n", T->nrow, s->th);
    fim_histogram_free(H);

    float * use = malloc(T->nrow*sizeof(float));
    assert(use != NULL);
    int vcol = ftab_get_col(T, "value");
    for(size_t kk = 0; kk<T->nrow; kk++)
    {
        use[kk] = T->T[T->ncol*kk + vcol] > s->th;
    }
    T = ftab_insert_col(T, use, "use");

    fim_t * fI = fim_image_from_array(feature, M, N, P);

    if(s->fwhm)
    {
        T = append_fwhm(s, T, fI, "fwhm_lateral_filtered", "fwhm_axial_filtered");
    }
    free(A);
    free(fI);
    free(feature);

    /* Discard unwanted dots */
    ftab_head(T, s->ndots);

    /* Write to file */
    if(s->verbose > 0)
    {
        printf("Writing dots to %s\n", s->outfile);
    }

    if(ftab_write_tsv(T, s->outfile))
    {
        fprintf(stderr, "Failed to write to %s\n", s->outfile);
        exit(EXIT_FAILURE);
    }

    fprintf(s->log, "Wrote %zu dots. Done!\n", T->nrow);
    //fprintf(s->log, "Optimal detection when \sigma = 0.425*FWHM\n");

    ftab_free(T);
    fclose(s->log);
    s->log = NULL;
    return;
}

int dw_dots(int argc, char ** argv)
{

    fim_tiff_init();
    opts * s = opts_new();


    argparsing(argc, argv, s);
    #ifdef _OPENMP
    omp_set_num_threads(s->nthreads);
    #endif
    if(s->verbose > 0)
    {
        opts_print(stdout, s);
    }


    for(size_t kk = s->optpos; kk < (size_t) argc; kk++)
    {
        if(s->verbose > 0)
        {
            printf("Processing file %d/%d : %s\n",
                   (int) (kk - s->optpos + 1),
                   argc - (int) s->optpos,
                   argv[kk]);
        }
        detect_dots(s, argv[kk]);
    }

    opts_free(s);
    return EXIT_SUCCESS;
}
