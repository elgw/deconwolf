#include "dw_dots.h"

typedef struct{
    int overwrite;
    int verbose;
    int optpos;
    char * image; /* Image to analyze */
    char * outfile; /* Where to write the tsv output */
    char * fout; /* Where to (optionally) write the filtered image */
    int nthreads;
    /* This determines the filter size */
    float log_asigma; /* Axial sigma for LoG filter */
    float log_lsigma; /* Lateral sigma */
    float fit_asigma; /* For dot fitting */
    float fit_lsigma;
    float NA;
    float ni;
    float lambda;
    float dx;
    float dz;
    float swell;
    int ndots; /* Number of dots to export */
    char * logfile;
    FILE * log;
    int fwhm;
    int fitting; /* Set to 1 to enable fitting */
    float th;
    int optind;
    char * cmdline;
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
    s->fwhm = 0;
    s->swell = 1;
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
    if(s->fwhm)
    {
        fprintf(f, "Will calculate FWHM\n");
    } else {
        fprintf(f, "No FWHM calculations\n");
    }
    if(s->fitting)
    {
        fprintf(f, "Will fit the dots:\n");
        fprintf(f, "   Lateral sigma: %.2f\n", s->fit_lsigma);
        fprintf(f, "   Axial sigma: %.2f\n", s->fit_asigma);
    } else {
        fprintf(f, "Fitting not enabled\n");
    }
    fprintf(f, "NA: %f, ni: %f, lambda: %f, dx: %f, dz: %f\n", s->NA, s->ni, s->lambda, s->dx, s->dz);
    return;
}

static void usage(__attribute__((unused)) int argc, char ** argv)
{
    opts * s = opts_new();
    printf("Detection of diffraction limited dots in 3D images\n"
           "using a Laplacian of Gaussian filter. Optional fitting\n"
           "using a Gaussian model\n");
    printf("\n");
    printf("usage: %s [<options>] input.tif input2.tif ...\n", argv[0]);
    printf("\n");
    printf("Recommended/required arguments:\n");
    printf(" --NA NA\n\t Set numerical aperture\n");
    printf(" --ni ni\n\t Set refractive index\n");
    printf(" --dx dx\n\t Lateral pixel size\n");
    printf(" --dz dz\n\t Axial pixel size\n");
    printf(" --lambda l\n\t Emission wave length\n");
    printf(" --ndots n\n\t Number of dots to export (default %d)\n", s->ndots);
    printf("\n");
    printf("Additional options\n");
    printf(" --nscale n\n"
           "\t set the number of scales to use\n");
    printf(" --swell f\n\t"
           "Tell the program how much larger the dots are compared to\n"
           "the diffraction limit. Default = 1, i.e. diffraction limited dots\n"
           "For some experiments values up to 2 makes sense\n");
    printf(" --overwrite\n\t Overwrite existing files (default %d)\n", s->overwrite);
    printf(" --help\n\t Show this message\n");
    printf(" --logfile file.txt\n\t Specify where the log file should be written\n");
    printf(" --fwhm\n"
           "\tInclude FWHM based on interpolation in the output (default %d)\n",
           s->fwhm);
    printf(" --verbose v\n\t Verbosity level (default %d)\n", s->verbose);
    printf(" --nthreads n\n\t Set the number of computational threads\n");
    printf(" --fout file.tif\n\t Write filtered image -- for debugging\n");
    printf("\n");
    printf("If you want to control the filter sizes, skip the optical parameters "
           "above and set the filter sizes manually by:\n");
    printf(" --dog_ls s\n\t Lateral sigma (location of zero-crossing)\n");
    printf(" --dog_as s\n\t Axial sigma (location of zero-crossing)\n");
    printf(" --fit_ls\n\t"
           "Lateral sigma, initial guess for the dot fitting");
    printf(" --fit_as\n\t"
           "Axial sigma, initial guess for the dot fitting");
    printf("\n");
    printf("Notes:\n");
    printf(" - Log messages will be written to [input file].log.txt\n");
    printf(" - Dots will be exported to [input file].dots.tsv\n");
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
        {"dog_as", required_argument, NULL, 'a'},
        {"fit_as", required_argument, NULL, 'A'},
        {"logfile", required_argument, NULL, 'w'},
        {"out", required_argument, NULL, 'O'},

        {"dog_ls", required_argument, NULL, 'l'},
        {"fit_ls", required_argument, NULL, 'l'},

        {"fwhm", no_argument, NULL, 'f'},
        {"fitting", no_argument, NULL, 'F'},
        {"help", no_argument, NULL, 'h'},

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
    while((ch = getopt_long(argc, argv, "1:3:4:5:6:L:a:A:fF:hi:l:n:N:op:r:s:v:w:", longopts, NULL)) != -1)
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
        case 'f':
            s->fwhm = 1;
            break;
        case 'F':
            s->fitting = 1;
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
            break;
        case 'L':
            s->fit_lsigma = atof(optarg);
            break;
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
            s->swell = atof(optarg);
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

        /* Apply swelling correction */
        s->fit_lsigma *= s->swell;
        s->fit_asigma *= s->swell;
        s->log_lsigma *= s->swell;
        s->log_asigma *= s->swell;
    }

    if(s->verbose > 1)
    {
        printf("VERBOSE>1: Printing out the settings just before validation:\n");
        opts_print(stdout, s);
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

    if(s->nscale > 1)
    {
        s->scales = calloc(s->nscale, sizeof(float));
        double f = exp( log(s->max_rel_scale)/ (double) (s->nscale -1.0));
        for(int kk = 0; kk<s->nscale; kk++)
        {
            s->scales[kk] = pow(f, kk);
        }
    }

    s->optpos = optind;
    s->optpos = optind;
    return;
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

    /* 2. Run the fitting */
    gmlfit * config = gmlfit_new();
    assert(config != NULL);
    config->image = I;
    config->M = M;
    config->N = N;
    config->P = P;
    config->sigma_xy = s->fit_lsigma;
    config->sigma_z = s->fit_asigma;
    config->log = s->log;
    config->verbose = s->verbose;
    config->X = X;
    config->nX = T->nrow;
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

    /* 3. Insert into the table */
    if(F == NULL)
    {
        fprintf(stderr, "Unable to run dot fitting. This is a bug. Please file a report\n");
        goto fail1;
    }

    printf("Concatenating tables\n");
    ftab_t * TF = ftab_new(10);
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
    TF->nrow = T->nrow;
    TF->nrow_alloc=T->nrow;
    assert(TF->ncol == 10);
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
    int fxcol = ftab_get_col(TT, "f_x");
    int fycol = ftab_get_col(TT, "f_y");
    int fzcol = ftab_get_col(TT, "f_z");

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
    {
        float scaling = dw_read_scaling(inFile);
        if(scaling != 1)
        {
            if(s->verbose > 0)
            {
                printf("Scaling by %f\n", 1.0/scaling);
                fim_mult_scalar(A, M*N*P, 1.0/scaling);
            }
        }
    }

    if(s->ndots < 0)
    {
        s->ndots = 0.025 * (float) M * (float) N;
        fprintf(s->log, "Setting the max number of output dots to %d\n",
                s->ndots);
    }

    float * feature = NULL;

    if(s->verbose > 1)
    {
        printf("LoG filter, lsigma=%.2f asigma=%.2f\n",
               s->log_lsigma, s->log_asigma);
    }
    //fim_set_verbose(2);

    ftab_t * T = NULL;
    if(s->nscale < 2)
    {
        feature = fim_LoG_S2(A, M, N, P, s->log_lsigma, s->log_asigma);


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
        float ** LoG = calloc(s->nscale, sizeof(float**));
        for(int ss = 0; ss < s->nscale; ss++)
        {
            float scaling = s->scales[ss];
            if(s->verbose > 0)
            {
                printf("LoG filter %d/%d, sigma = %f, %f\n",
                       ss+1, s->nscale,
                       scaling*s->log_lsigma,
                       scaling*s->log_asigma);
            }

            LoG[ss] = fim_LoG_S2(A, M, N, P,
                                 scaling*s->log_lsigma,
                                 scaling*s->log_asigma);

            float s2 = scaling*scaling;
            #pragma omp parallel for
            for(size_t kk = 0; kk < M*N*P; kk++)
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

    if(s->fwhm)
    {
        fim_t * fI = fim_image_from_array(A, M, N, P);
        T = append_fwhm(s, T, fI, "fwhm_lateral_filtered", "fwhm_axial_filtered");
        fimt_free(fI);
    }

    if(s->fitting)
    {
        T = append_fitting(s, T,
                           A, M, N, P);
    }

    free(A);



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
        printf("Peak memory usage: %zu kB\n",
               get_peakMemoryKB());
    }

    opts_free(s);

    return EXIT_SUCCESS;
}
