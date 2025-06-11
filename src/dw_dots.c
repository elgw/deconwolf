#include "dw_dots.h"

// TODO:
// - '--multiscale' as the only option for multiscaling ...
// --dots option for nd2tool
/* Possible improvements:
 * - Code path for 2D images (max projections etc)
 * - Separate Dot vs edge/ridge using the Hessian or the structure tensor.
 * - Use GPU
 */

typedef struct{
    int overwrite;
    int verbose;
    int optpos;
    char * image; /* Image to analyze */
    char * outfile; /* Where to write the tsv output */
    char * fout; /* Where to (optionally) write the filtered image */
    int nthreads;

    /* Diffraction limited size of dots*/
    float fit_asigma;
    float fit_lsigma;

    /* Size of the LoG filter */
    float log_asigma; /* Axial sigma for LoG filter */
    float log_lsigma; /* Lateral sigma */

    /* Optical configuration */
    float NA;
    float ni;
    float lambda;
    float dx;
    float dz;

    int ndots; /* Number of dots to export */
    char * logfile;
    FILE * log;
    int fitting; /* Set to 1 to enable fitting */
    float th;
    int optind;
    char * cmdline;
    int write_csv;

    int circularity; /* Set to 1 to enable circularity estimation */

    /* For multi scale dot detection.
     *For single scale paths,
     * scales[0] is also used to scale up filter sizes if set. */
    int nscale;
    float max_rel_scale;
    float * scales;
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
    s->fitting = 0;
    s->nscale = 1;
    s->max_rel_scale = 2;
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
    free(s->cmdline);
    free(s->scales);
    free(s);
}

static void opts_print(FILE * f, opts * s)
{
    assert(s != NULL);
    assert(f != NULL);
    fprintf(f, "CMD: %s\n", s->cmdline);
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

    fprintf(f, "Assuming that the dots look like Gaussians with\n"
            "a lateral sigma of %f and an axial sigma of %f\n",
            s->fit_lsigma, s->fit_asigma);
    fprintf(f, "Dot detection method: Laplacian of Gaussian (LoG)\n");
    fprintf(f, "   Lateral sigma: %.2f pixels\n", s->log_lsigma);
    fprintf(f, "   Axial sigma: %.2f pixels\n", s->log_asigma);

    if(s->nscale < 2)
    {
        fprintf(f, "Multiscale: No\n");
    } else {
        fprintf(f, "Multiscale: Yes, %d scales\n", s->nscale);
    }

    if(s->ndots > 0)
    {
        fprintf(f, "Will write at most %d dots\n", s->ndots);
    } else {
        fprintf(f, "Automatic number of dots in the output\n");
    }

    if(s->fitting)
    {
        fprintf(f, "Will fit the dots:\n");
        fprintf(f, "   Lateral sigma: %.2f\n", s->fit_lsigma);
        fprintf(f, "   Axial sigma: %.2f\n", s->fit_asigma);
    } else {
        fprintf(f, "Fitting not enabled\n");
    }

    fprintf(f, "Output format: ");
    if(s->write_csv)
    {
        fprintf(f, "CSV\n");
    } else {
        fprintf(f, "TSV\n");
    }
    fprintf(f, "NA: %f, ni: %f, lambda: %f, dx: %f, dz: %f\n", s->NA, s->ni, s->lambda, s->dx, s->dz);
    return;
}

static void usage(__attribute__((unused)) int argc, char ** argv)
{
    opts * s = opts_new();
    printf("Detection of diffraction limited dots in 3D images\n"
           "using a Laplacian of Gaussian filter. Optional fitting\n"
           "using a Gaussian model\n"
           "Limited functionality for 2D images\n");
    printf("\n");
    printf("usage: %s [<options>] input.tif input2.tif ...\n", argv[0]);
    printf("\n");
    printf("Recommended/required arguments:\n");
    printf("  --NA NA\n\t"
           "Set numerical aperture\n");
    printf("  --ni ni\n\t"
           "Set refractive index\n");
    printf("  --dx dx\n\t"
           "Lateral pixel size [nm]\n");
    printf("  --dz dz\n\t"
           "Axial pixel size [nm]\n");
    printf("  --lambda l\n\t"
           "Emission wave length [nm]\n");
    printf("  --ndots n\n\t"
           "Number of dots to export (default M x N x 0.005)\n");
    printf("\n");
    printf("Additional options\n");
    printf("  --nscale n\n\t"
           "set the number of scales to use\n");
    printf("  --swell f\n\t"
           "Tell the program how much larger the dots are compared to\n"
           "the diffraction limit. Default = 1, i.e. diffraction limited dots\n"
           "For some experiments values up to 2 makes sense\n");
    printf("  --overwrite\n\t"
           "Overwrite existing files (default %d)\n",
           s->overwrite);
    printf("  --help\n\t"
           "Show this message\n");
    printf("  --logfile file.txt\n\t"
           "Specify where the log file should be written\n");

    printf("  --verbose v\n\t"
           "Verbosity level (default %d)\n", s->verbose);
    printf("  --nthreads n\n\t"
           "Set the number of computational threads\n");
    printf("  --fout file.tif\n\t"
           "Write filtered image -- for debugging\n");
    printf("\n");
    printf("If you want to control the filter sizes, skip the optical parameters\n"
           "above and set the filter sizes manually by:\n");
    printf("  --log_ls s\n\t"
           "Lateral sigma (location of zero-crossing)\n");
    printf("  --log_as s\n\t"
           "Axial sigma (location of zero-crossing)\n");
    printf("  --fit_ls\n\t"
           "Lateral sigma, initial guess for the dot fitting\n");
    printf("  --fit_as\n\t"
           "Axial sigma, initial guess for the dot fitting");
    printf("\n");
    printf("Notes:\n");
    printf("  - Log messages will be written to [input file].log.txt\n");
    printf("  - Dots will be exported to [input file].dots.tsv\n");
    printf("\n");
    printf(" See the man page for more information.\n");
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
    float * row = malloc((T->ncol+1)*sizeof(float));
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
    size_t cmdline_size = 2 + 3*argc;
    for(int kk = 0; kk < argc; kk++)
    { cmdline_size += strlen(argv[kk]); }
    s->cmdline = calloc(cmdline_size, 1);
    assert(s->cmdline != NULL);
    for(int kk = 0; kk < argc; kk++)
    {
        strncat(s->cmdline, "'", cmdline_size);
        strncat(s->cmdline, argv[kk], cmdline_size);
        strncat(s->cmdline, "' ", cmdline_size);
    }


    struct option longopts[] = {
        {"log_as", required_argument, NULL, 'a'},
        {"fit_as", required_argument, NULL, 'A'},
        {"csv",    no_argument, NULL, 'c'},
        {"circularity", no_argument, NULL, 'C'},
        {"logfile", required_argument, NULL, 'w'},
        {"out", required_argument, NULL, 'O'},
        {"log_ls", required_argument, NULL, 'L'},
        {"fit_ls", required_argument, NULL, 'l'},
        {"fitting", no_argument, NULL, 'F'},
        {"help", no_argument, NULL, 'h'},
        {"max_scale", required_argument, NULL, 'm'},
        {"ndots",   required_argument, NULL, 'n'},
        {"nscale",  required_argument, NULL, 'N'},
        {"overwrite", no_argument, NULL, 'o'},
        {"fout",     required_argument, NULL, 'p'},
        {"swell",   required_argument, NULL, 's'},
        {"threads", required_argument, NULL, 't'},
        {"verbose", required_argument, NULL, 'v'},
        {"lambda", required_argument, NULL, '1'},
        {"NA",     required_argument, NULL, '3'},
        {"dx",     required_argument, NULL, '4'},
        {"dz",     required_argument, NULL, '5'},
        {"ni",     required_argument, NULL, '6'},
        {NULL, 0, NULL, 0}};
    int ch;
    while((ch = getopt_long(argc, argv, "1:3:4:5:6:L:a:A:cCF:hi:l:L:m:n:N:op:r:s:v:w:", longopts, NULL)) != -1)
    {
        switch(ch){
        case '1':
            s->lambda = atof(optarg);
            break;
        case '3':
            s->NA = atof(optarg);
            break;
        case '4':
            s->dx = atof(optarg);
            break;
        case '5':
            s->dz = atof(optarg);
            break;
        case '6':
            s->ni = atof(optarg);
            break;
        case 'a':
            s->log_asigma = atof(optarg);
            break;
        case 'A':
            s->fit_asigma = atof(optarg);
            break;
        case 'c':
            s->write_csv = 1;
            break;
        case 'C':
            s->circularity = 1;
            break;
        case 'F':
            s->fitting = 1;
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
            break;
        case 'L':
            s->log_lsigma = atof(optarg);
            break;
        case 'l':
            s->fit_lsigma = atof(optarg);
            break;
        case 'm':
            s->max_rel_scale = atof(optarg);
            break;
        case 'n':
            s->ndots = atoi(optarg);
            break;
        case 'N':
            s->nscale = atoi(optarg);
            break;
        case 'o':
            s->overwrite = 1;
            break;
        case 'O':
            free(s->outfile);
            s->outfile = strdup(optarg);
            assert(s->outfile != NULL);
            break;
        case 'p':
            s->fout = strdup(optarg);
            assert(s->fout != NULL);
            break;
        case 's':
            free(s->scales);
            s->scales = calloc(1, sizeof(float));
            assert(s->scales != 0);
            s->scales[0] = atof(optarg);
            s->nscale = 1;
            break;
        case 't':
            s->nthreads = atoi(optarg);
            break;
        case 'v':
            s->verbose = atoi(optarg);
            break;
        case 'w':
            free(s->logfile);
            s->logfile = strdup(optarg);
            assert(s->logfile != NULL);
            break;
        default:
            exit(EXIT_FAILURE);
        }
    }

    if(s->verbose > 1)
    {
        printf("VERBOSE>1: Printing out the settings just before validation:\n");
        opts_print(stdout, s);
    }


    if(s->NA*s->ni*s->dx*s->dz*s->lambda > 0)
    {
        float fwhm_pixels = abbe_res_xy(s->lambda, s->NA)/s->dx;
        float fwhm_pixels_z = abbe_res_z(s->lambda, s->NA) / s->dz;
        /* To determine the initial guess for the fitting,
         * we assume a Gaussian signal and convert the
         * fwhm to sigma
         */
        s->fit_lsigma = fwhm_pixels / (2.0*sqrt(2.0*log(2.0)));
        s->fit_asigma = fwhm_pixels_z / (2.0*sqrt(2.0*log(2.0)));
        /* The LoG achieve the maximum response when sigma = r/sqrt(2)
         * for disks. For Gaussians it seems that the max response is given
         * when sigma_LoG = sqrt(2)*sigma_sigmal
         */
        s->log_lsigma = s->fit_lsigma*sqrt(2.0);
        s->log_asigma = s->fit_asigma*sqrt(2.0);
    }


    if( (s->log_asigma)*(s->log_lsigma) <= 0 )
    {
        fprintf(stderr, "ERROR: "
                "Not enough parameters specified to determine the LoG filter size\n");
        fprintf(stderr, "Please specify\n");
        fprintf(stderr, "          --log_ls and --log_as\n");
        fprintf(stderr, "      OR\n");
        fprintf(stderr, "          --NA, --ni, --lambda, --dx and --dz\n");
        exit(EXIT_FAILURE);
    }

    if(s->fitting)
    {
        if( (s->fit_asigma)*(s->fit_lsigma) <= 0)
        {
            fprintf(stderr, "ERROR: "
                    "Not enough parameters specified to determine the initial spot size for fitting\n");
            fprintf(stderr, "Please specify\n");
            fprintf(stderr, "Please specify\n");
            fprintf(stderr, "          --fit_ls and --fit_as\n");
            fprintf(stderr, "      OR\n");
            fprintf(stderr, "          --NA, --ni, --lambda, --dx and --dz\n");
            exit(EXIT_FAILURE);
        }
    }

    // TODO: if s->multiscale. Also scale factor as a parameter.
    if(s->nscale > 1)
    {
        /* scale == 1 means at the scale of diffraction limited dots
         * deduced from the optical parameters.
         *
         * The multi scale detection is performed with sizes with consecutive
         * ratio of sqrt(2) (same number again :)). This is because
         * for some reason we like to have LoG filter at four scales:
         * a) 1/sqrt(2), b) 1, c) sqrt(2), d) 2
         *
         * The LoG filter has a sigma sqrt(2) times larger than the dots.
         * it performs quite poor when it is below 1 i.e. when the dots
         * has a sigma < 1/sqrt(2) in the axial or lateral direction.
         * If that occurs, we start at at the first possible scale, and
         * continue with as many as possible until we reach double the
         * theoretical dot size.
         */

        /* What we would like to use */
        const float min_scale_ideal = 1.0/sqrt(2.0);
        float min_scale = min_scale_ideal;

        const float max_scale = s->max_rel_scale;
        const float f = sqrt(2.0); /* ratio between scales */

        /* See what lowest scale we can use ...  */
        float min_fit_sigma = 1.0/sqrt(2.0);
        if(s->fit_lsigma*min_scale < min_fit_sigma)
        {
            min_scale = min_fit_sigma / s->fit_lsigma;
        }
        if(s->fit_asigma*min_scale < min_fit_sigma)
        {
            min_scale = min_fit_sigma / s->fit_asigma;
        }
        if(min_scale > 1.0/sqrt(2.0))
        {
            if(s->verbose > 0)
            {
                printf("Warning: Couldn't start at scale %f due to the sampling of\n"
                       "the image. Starting at scale %f\n",
                       1.0/sqrt(2), min_scale);
            }
        }

        /* It is possible that the max scale is problematic, typically when
         * there are too few z planes. In that case it might be better to do 2D
         * dot detection ... */

        s->scales = calloc(s->nscale, sizeof(float));
        s->nscale = ceil( (log(max_scale) - log(min_scale)) / log(f));
        if(s->verbose > 0)
        {
            printf("Will use %d scales\n", s->nscale);
        }
        for(int kk = 0; kk<s->nscale; kk++)
        {
            s->scales[kk] = min_scale*pow(f, kk);
        }
    }

    s->optpos = optind;
    return;
}

static ftab_t *
append_circularity(opts * s, ftab_t * T, const float * restrict I,
                   size_t M, size_t N, size_t P)
{
    if(s->verbose > 2)
    {
        printf("append_circularity()\n");
    }
    int xcol = ftab_get_col(T, "f_x");
    int ycol = ftab_get_col(T, "f_y");
    int zcol = ftab_get_col(T, "f_z");
    int sigma_col = ftab_get_col(T, "f_sigma_lateral");
    if(xcol < 0)
    {
        printf("Error: no sub pixel locations available for circularity estimates\n");
        return T;
    }

    assert(xcol >= 0);
    assert(ycol >= 0);
    assert(zcol >= 0);

    ftab_t * TC = ftab_new(1);
    ftab_set_colname(TC, 0, "lat_circularity");
    free(TC->T);
    TC->nrow = T->nrow;
    TC->T = calloc(T->nrow, sizeof(float));
    assert(TC->T != NULL);

    #pragma omp parallel for
    for(size_t kk = 0; kk < T->nrow; kk++)
    {
        float * row = T->T + kk*T->ncol;
        double x = row[xcol];
        double y = row[ycol];
        double z = row[zcol];
        double sigma = row[sigma_col];
        TC->T[kk] = fim_dot_lateral_circularity(I,
                                                M, N, P,
                                                x,y,z,
                                                sigma);
    }


    ftab_t * TT = ftab_concatenate_columns(T, TC);

    ftab_free(T);
    ftab_free(TC);


    return TT;
}

static ftab_t * append_fitting(opts * s, ftab_t * T, float * I,
                               size_t M, size_t N, size_t P)
{

    /* 1. Extract the start coordinates */
    int xcol = ftab_get_col(T, "x");
    int ycol = ftab_get_col(T, "y");
    int zcol = ftab_get_col(T, "z");

    assert(xcol >= 0);
    assert(ycol >= 0);
    assert(zcol >= 0);

    double * X = calloc(3*T->nrow, sizeof(double));
    assert(X != NULL);

    for(size_t kk = 0; kk < T->nrow; kk++)
    {
        float * row = T->T + kk*T->ncol;
        X[3*kk+0] = row[xcol];
        X[3*kk+1] = row[ycol];
        X[3*kk+2] = row[zcol];
    }

    /* Per dot scaling */
    double * DS = NULL;
    int row_id_scale = ftab_get_col(T, "LoG_scale");
    if(s->verbose > 1)
    {
        printf("row_id_scale = %d\n", row_id_scale);
    }

    if(row_id_scale > -1)
    {
        DS = calloc(T->nrow, sizeof(double));
        assert(DS != NULL);

        for(size_t kk = 0; kk < T->nrow; kk++)
        {
            float * row = T->T + kk*T->ncol;
            DS[kk] = row[row_id_scale];
        }
    }

    float scaling = 1;
    if(s->scales != NULL)
    {
        scaling = s->scales[0];
    }

    /* 2. Run the fitting */
    gmlfit * config = gmlfit_new();
    assert(config != NULL);
    config->image = I;
    config->M = M;
    config->N = N;
    config->P = P;
    config->sigma_xy = scaling*s->fit_lsigma;
    config->sigma_z = scaling*s->fit_asigma;
    config->log = s->log;
    config->verbose = s->verbose;
    config->X = X;
    config->nX = T->nrow;
    config->DS = DS;
    if(s->verbose > 0)
    {
        printf("Fitting %zu dots\n", config->nX);
    }
    if(s->log != NULL)
    {
        fprintf(s->log, "Fitting %zu dots\n", config->nX);
    }

    double * F = gmlfit_run(config);
    free(config);
    free(X);
    free(DS);

    /* 3. Insert into the table */
    if(F == NULL)
    {
        fprintf(stderr, "Unable to run dot fitting. This is a bug. Please file a report\n");
        goto fail1;
    }

    printf("Concatenating tables\n");
    ftab_t * TF = ftab_new(11);
    free(TF->T);
    ftab_set_colname(TF, 0, "f_bg");
    ftab_set_colname(TF, 1, "f_signal_count");
    ftab_set_colname(TF, 2, "f_peak_signal");
    ftab_set_colname(TF, 3, "f_x");
    ftab_set_colname(TF, 4, "f_y");
    ftab_set_colname(TF, 5, "f_z");
    ftab_set_colname(TF, 6, "f_sigma_lateral");
    ftab_set_colname(TF, 7, "f_sigma_axial");
    ftab_set_colname(TF, 8, "f_status");
    ftab_set_colname(TF, 9, "f_error");
    ftab_set_colname(TF, 10, "f_corr");
    TF->nrow = T->nrow;
    TF->nrow_alloc=T->nrow;
    assert(TF->ncol == 11);
    TF->T = calloc(TF->nrow*TF->ncol, sizeof(float));
    assert(TF->T != NULL);
    for(size_t kk = 0; kk< TF->nrow*TF->ncol; kk++)
    {
        TF->T[kk] = F[kk];
    }
    free(F);

    ftab_t * TT = ftab_concatenate_columns(T, TF);
    ftab_free(TF);
    ftab_free(T);

    xcol = ftab_get_col(TT, "x");
    ycol = ftab_get_col(TT, "y");
    zcol = ftab_get_col(TT, "z");
    assert(xcol != ycol);
    assert(xcol != zcol);

    int fxcol = ftab_get_col(TT, "f_x");
    int fycol = ftab_get_col(TT, "f_y");
    int fzcol = ftab_get_col(TT, "f_z");
    assert(fxcol != fycol);
    assert(fxcol != fzcol);

    /* Make the fitted positions absolute */
    for(size_t kk = 0; kk < TT->nrow; kk++)
    {
        float * row = TT->T + kk*TT->ncol;
        row[fxcol]+=row[xcol];
        row[fycol]+=row[ycol];
        row[fzcol]+=row[zcol];
    }

    return TT;


fail1: ;
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
    if(!dw_isfile(inFile))
    {
        fprintf(stderr, "Can't open %s!\n", inFile);
        return;
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
    if(s->write_csv)
    {
        sprintf(s->outfile, "%s.dots.csv", s->image);
    } else {
        sprintf(s->outfile, "%s.dots.tsv", s->image);
    }

    char * outFile = s->outfile;

    if(s->verbose > 1)
    {
        fprintf(stdout, "Input file: %s\n", inFile);
    }

    if(s->verbose > 1)
    {
        fprintf(stdout, "Ouput file: %s\n", outFile);
    }

    if(s->overwrite == 0 && dw_isfile(outFile))
    {
        printf("%s exists, skipping.\n", outFile);
        return;
    }

    fprintf(s->log, "Reading %s\n", inFile);
    int64_t M = 0, N = 0, P = 0;
    float * A = fim_tiff_read(inFile, NULL, &M, &N, &P, s->verbose);
    {
        float scaling = dw_read_scaling(inFile);
        if(scaling != 1)
        {
            if(s->verbose > 0)
            {
                printf("Scaling by %f\n", 1.0/scaling);
                fprintf(s->log, "Scaling by %f\n", 1.0/scaling);
                fim_mult_scalar(A, M*N*P, 1.0/scaling);
            }
        } else {
            fprintf(s->log, "Could not read a scaling value from %s.log.txt\n", inFile);
        }
    }

    if(s->verbose > 1)
    {
        printf("Image size: %ld x %ld x %ld\n", M, N, P);
    }

    if(s->ndots < 0)
    {
        s->ndots = 0.005 * (float) M * (float) N;
        fprintf(s->log, "Setting the max number of output dots to %d\n",
                s->ndots);
    }

    float * feature = NULL;


    //fim_set_verbose(2);

    ftab_t * T = NULL;
    if(s->nscale < 2)
    {
        if(s->verbose > 0)
        {
            printf("Diffraction limited dots: sigma = %.2f, %.2f\n",
                   s->fit_lsigma, s->fit_asigma);
        }
        float scaling = 1;
        if(s->scales != NULL)
        {
            scaling = s->scales[0];
        }
        if(scaling != 1)
        {
            if(s->verbose > 0)
            {
                printf("Size multiplier: %.2f -> sigma = %.2f, %.2f\n",
                       scaling, scaling*s->fit_lsigma, scaling*s->fit_asigma);
            }
        }

        if(s->verbose > 1)
        {
            printf("LoG filter, LoG_lsigma=%.2f LoG_asigma=%.2f\n",
                   scaling*s->log_lsigma, scaling*s->log_asigma);

        }

        feature = fim_LoG_S2(A, M, N, P,
                             scaling*s->log_lsigma,
                             scaling*s->log_asigma);

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

        T = fim_lmax(feature, M, N, P);
        free(feature);
        feature = NULL;
    }


    if(s->nscale > 1)
    {
        /* use s->nscale, starting at factor 1 (for diffraction limited dots)
         * and continue up till max_rel_scale
         * f^(nscale-1) = max_rel_scale;
         *
         */
        float ** LoG = calloc(s->nscale, sizeof(float*));
        assert(LoG != NULL);
        for(int ss = 0; ss < s->nscale; ss++)
        {
            float scaling = s->scales[ss];
            if(s->verbose > 0)
            {
                printf("LoG filter %d/%d, sigma = %f, %f (%f x)\n",
                       ss+1, s->nscale,
                       scaling*s->log_lsigma,
                       scaling*s->log_asigma,
                       scaling);
            }
            assert(P > 0);
            LoG[ss] = fim_LoG_S2(A, M, N, P,
                                 scaling*s->log_lsigma,
                                 scaling*s->log_asigma);

            float s2 = scaling*scaling;
#pragma omp parallel for
            for(size_t kk = 0; kk < (int64_t) M*N*P; kk++)
            {
                LoG[ss][kk] *= s2;
            }

        }
        printf("Multiscale maxima detection\n");

        T = fim_lmax_multiscale(LoG, s->scales, s->nscale, M, N, P);
        for(int kk = 0 ; kk < s->nscale; kk++)
        {
            free(LoG[kk]);
        }
        free(LoG);
    }

    if(T == NULL)
    {
        printf("Unable to continue, table could not be read\n");
        exit(EXIT_FAILURE);
    }

    if(s->verbose > 1)
    {
        printf("Found %zu points\n", T->nrow);
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
        printf("Extracted %zu dots: Mean=%f, std=%f\n", T->nrow, mean, std);
    }
    fprintf(s->log, "Extracted %zu dots: Mean=%f, std=%f\n", T->nrow, mean, std);

    fim_histogram_t * H = fim_histogram(values, T->nrow);
    free(values);
    //s->th = fim_histogram_otsu(H);
    if(H != NULL)
    {
        fim_histogram_log(H);
        s->th = fim_histogram_otsu(H);
    } else {
        fprintf(stderr, "Warning: could not create a histogram for the dots\n");
        s->th = -1;
    }
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
    free(use);

    /* Discard unwanted dots before the computationally demanding fitting */
    ftab_head(T, s->ndots);

    if(s->fitting)
    {
        T = append_fitting(s, T,
                           A, M, N, P);
    }

    /* Will use sub pixel locations if fitting was performed */
    if(s->circularity)
    {
        T = append_circularity(s, T, A, M, N, P);
    }

    free(A);



    /* Write to file */
    if(s->verbose > 0)
    {
        printf("Writing dots to %s\n", s->outfile);
    }

    if(s->write_csv == 1)
    {
        if(ftab_write_csv(T, s->outfile))
        {
            fprintf(stderr, "Failed to write to %s\n", s->outfile);
            exit(EXIT_FAILURE);
        }
    } else {
        if(ftab_write_tsv(T, s->outfile))
        {
            fprintf(stderr, "Failed to write to %s\n", s->outfile);
            exit(EXIT_FAILURE);
        }
    }

    fprintf(s->log, "Wrote %zu dots. Done!\n", T->nrow);

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
    if(s->verbose > 1)
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

    if(s->verbose > 1)
    {
        fprint_peak_memory(stdout);
    }

    opts_free(s);

    return EXIT_SUCCESS;
}
