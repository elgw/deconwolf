#include "dw_render.h"

enum DOT_MODE { DOT_MODE_AUTO, DOT_MODE_TH, DOT_MODE_NDOTS};

typedef struct{
    int overwrite;
    int verbose;
    int optpos;
    char * image; /* Image to analyze */
    char * outfile; /* Where to write the tsv output */
    char * dotfile; /* TSV file */
    int nthreads;
    float lsigma; /* Lateral sigma */
    int ndots; /* Number of dots to export */
    /* Mapping of image intensities to black, white */
    enum DOT_MODE dot_mode;
    float plow; /* Lower percentile for image mapping */
    float phigh; /* Higher percentile for image mapping */
    float th; /* Dot threshold, only used when --th is specified */
    int drawtext;
} opts;

static opts * opts_new();
static void opts_free(opts * s);
static void usage(__attribute__((unused)) int argc, char ** argv);
static void argparsing(int argc, char ** argv, opts * s);
static int file_exist(char * fname);

static opts * opts_new()
{
    opts * s = malloc(sizeof(opts));

    s->overwrite = 0;
    s->verbose = 1;
    s->optpos = -1;
    s->image = NULL;
    s->outfile = NULL;
    s->dotfile = NULL;
    s->nthreads = dw_get_threads();

    s->lsigma = 0;

    s->ndots = -1; /* Auto */
    s->th = -1;
    s->dot_mode = DOT_MODE_AUTO;
    /* For stretching the image data */
    s->plow = 0.01;
    s->phigh = 1.0-0.005;
    s->drawtext = 1;
    return s;
}


static void opts_free(opts * s)
{
    dw_nullfree(s->outfile);
    dw_nullfree(s->image);
    dw_nullfree(s->dotfile);
    free(s);
}

static void opts_print(FILE * f, opts * s)
{
    fprintf(f, "overwrite = %d\n", s->overwrite);
    fprintf(f, "verbose = %d\n", s->verbose);
    fprintf(f, "image: %s\n", s->image);
    fprintf(f, "outfile: %s\n", s->outfile);

    fprintf(f, "nthreads = %d\n", s->nthreads);
    fprintf(f, "Lateral sigma: %.2f\n", s->lsigma);
    switch(s->dot_mode)
    {
    case DOT_MODE_AUTO:
        fprintf(f, "DOT_MODE: Automatic\n");
        break;
    case DOT_MODE_TH:
        fprintf(f, "DOT_MODE: Threshold (>= %f)\n", s->th);;
        break;
    case DOT_MODE_NDOTS:
        fprintf(f, "DOT_MODE: %d dots\n", s->ndots);
        break;
    }

    if(s->ndots > 0)
    {
        fprintf(f, "Will write at most %d dots\n", s->ndots);
    } else {
        fprintf(f, "Automatic number of dots in the output\n");
    }

    return;
}

static void usage(__attribute__((unused)) int argc, char ** argv)
{
    opts * s = opts_new();
    printf("usage: %s [<options>] --image input.tif\n", argv[0]);
    printf("Options:\n");
    printf(" --image im.tif\n\t Image to load\n");
    printf(" --dots file.tsv\n\t Specify TSV file to read, "
           "tries <image>.dots.tsv by default\n");
    printf("Image processing:\n");
    printf(" --lsigma s\n\t Lateral sigma (default %f)\n", s->lsigma);
    printf(" --plow p\n\t Lower percentile for image scaling (default %.4f)\n", s->plow);
    printf(" --phigh q\n\t Higher percentile for image scaling (default %.4f)\n", s->phigh);
    printf("Dot plotting\n");
    printf(" --ndots n\n\t Number of dots to export\n");
    printf(" --th th\n\t Set dot threshold\n");
    printf("General:\n");
    printf(" --overwrite\n\t Overwrite existing files\n");
    printf(" --help\n\t Show this message\n");
    printf(" --notext\n\t Don't write any text on the image\n");
    printf(" --verbose v\n\t Verbosity level\n");
    printf("\n");
    opts_free(s);
    return;
}

static void argparsing(int argc, char ** argv, opts * s)
{
    struct option longopts[] = {
        {"notext", no_argument, NULL, 'N'},
        {"out", required_argument, NULL, 'O'},
        {"th", required_argument, NULL, 'T'},
        {"dots",   required_argument, NULL, 'd'},
        {"help", no_argument, NULL, 'h'},
        {"image", required_argument, NULL, 'i'},
        {"lsigma", required_argument, NULL, 'l'},
        {"ndots",   required_argument, NULL, 'n'},
        {"overwrite", no_argument, NULL, 'o'},
        {"plow", required_argument, NULL, 'p'},
        {"phigh", required_argument, NULL, 'q'},
        {"threads", required_argument, NULL, 't'},
        {"verbose", required_argument, NULL, 'v'},
        {NULL, 0, NULL, 0}};
    int ch;
    while((ch = getopt_long(argc, argv, "O:T:d:hi:l:n:op:q:r:v:", longopts, NULL)) != -1)
    {
        switch(ch){
        case 'N':
            s->drawtext = 0;
            break;
        case 'O':
            dw_nullfree(s->outfile);
            s->outfile = strdup(optarg);
            break;
        case 'T':
            s->th = atof(optarg);
            s->dot_mode = DOT_MODE_TH;
            break;
        case 'd':
            s->dotfile = strdup(optarg);
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
            break;
        case 'i':
            s->image = strdup(optarg);
            break;
        case 'l':
            s->lsigma = atof(optarg);
            break;
        case 'n':
            s->ndots = atoi(optarg);
            s->dot_mode = DOT_MODE_NDOTS;
            break;
        case 'o':
            s->overwrite = 1;
            break;
        case 'p':
            s->plow = atof(optarg);
            break;
        case 'q':
            s->phigh = atof(optarg);
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
    if(s->image == NULL)
    {
        fprintf(stderr, "No --image specified\n");
        exit(EXIT_FAILURE);
    }


    if(s->outfile == NULL)
    {
        s->outfile = malloc(strlen(s->image) + 64);
        sprintf(s->outfile, "%s.dots.svg", s->image);
    }

    if(s->dotfile == NULL)
    {
        s->dotfile = malloc(strlen(s->image) + 64);
        sprintf(s->dotfile, "%s.dots.tsv", s->image);
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



cairo_surface_t * fim_to_cairo_surface(fim_t * I,
                                             float low, float high)
{
    cairo_surface_t * result;
    unsigned char * current_row;
    int stride;


    result = cairo_image_surface_create(CAIRO_FORMAT_RGB24, I->M, I->N);
    if (cairo_surface_status(result) != CAIRO_STATUS_SUCCESS)
        return result;

    cairo_surface_flush(result);
    current_row = cairo_image_surface_get_data(result);
    stride = cairo_image_surface_get_stride(result);
    for (int y = 0; y < I->N; y++) {
        uint32_t *row = (void *) current_row;
        for (int x = 0; x < I->M; x++) {
            float v =  ((I->V[x + I->M*y])-low) / (high-low);
            v*=255;
            v > 255? v = 255 : 0;
            v < 0 ? v = 0 : 0;
            uint8_t value = v;
            uint32_t r = value;
            uint32_t g = value;
            uint32_t b = value;
            row[x] = (r << 16) | (g << 8) | b;
        }

        current_row += stride;
    }
    cairo_surface_mark_dirty(result);
    return result;
}

void render(opts * s, const fim_t * I, ftab_t * T)
{
    printf("Creating max projection\n");
    float * m = fim_maxproj(I->V, I->M, I->N, I->P);
    fim_t * M = fim_image_from_array(m, I->M, I->N, 1);
    free(m);

    if(s->lsigma > 0)
    {
        fim_gsmooth(M->V, M->M, M->N, 1, s->lsigma);
    }

    fim_histogram_t * H = fim_histogram(M->V, M->M*M->N*M->P);
    float low =  fim_histogram_percentile(H, s->plow);
    float high = fim_histogram_percentile(H, s->phigh);
    printf("Using range [%f, %f]\n", low, high);
    fim_histogram_free(H);
    cairo_surface_t * mip = fim_to_cairo_surface(M, low, high);

    cairo_surface_t * surf = cairo_svg_surface_create (s->outfile,
                                                       M->M,
                                                       M->N);
    //cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_NEAREST);
//    cairo_pattern_set_filter(mip, CAIRO_FILTER_NEAREST);

    cairo_t * cr = cairo_create(surf);

    cairo_set_source_surface(cr, mip, 0, 0);
    /* This does not do anything :( */
    cairo_pattern_set_filter (cairo_get_source (cr), CAIRO_FILTER_NEAREST);
    cairo_paint(cr);

    cairo_set_line_width(cr, 1);
    cairo_set_source_rgb(cr, 1.0, 0.00, 0);

    int xcol = ftab_get_col(T, "x");
    int ycol = ftab_get_col(T, "y");
    int vcol = ftab_get_col(T, "value");
    int ucol = ftab_get_col(T, "use");

    if(s->dot_mode == DOT_MODE_AUTO)
    {
        if(ucol == -1)
        {
            printf("Warning: No column named 'use' will show all dots\n");
        }
    }

    float minvalue = 1e99;
    int nshow = 0;
    if(T == NULL)
    {
        if(s->verbose > 0)
        {
            printf("No dot table loaded\n");
        }
    }
    if(T != NULL)
    {
        if(s->verbose > 1)
        {
            printf("Drawing dots\n");
        }
        for(size_t kk = 0; kk<T->nrow; kk++)
        {
            //printf("%f %f\n", T->T[kk*T->ncol], T->T[kk*T->ncol+1]);
            int draw = 0;
            float x = T->T[kk*T->ncol + xcol];
            float y = T->T[kk*T->ncol + ycol];
            float use = 1;
            if(ucol >= 0)
            {
                use = T->T[kk*T->ncol + ucol];
            }
            float value = T->T[kk*T->ncol + vcol];


            if(s->dot_mode == DOT_MODE_AUTO)
            {
                use == 1 ? draw = 1: 0;
            }
            if(s->dot_mode == DOT_MODE_TH)
            {
                value > s->th ? draw = 1 : 0;
            }
            if(s->dot_mode == DOT_MODE_NDOTS)
            {
                nshow < s->ndots ? draw = 1 : 0;
            }
            if(draw)
            {
                value < minvalue ? minvalue = value : 0;
                cairo_arc(cr, x, y, 5, 0, 2 * M_PI);
                cairo_close_path(cr);
                cairo_stroke(cr);
                nshow++;
            }
        }
    }
    printf("Drew %d/%zu dots\n", nshow, T->nrow);

    if(s->drawtext){
        size_t bufflen = 1024;
        char * buff = malloc(bufflen);
        float ypos = 30;
        float dy = 30;
        float xpos = 20;
        cairo_set_source_rgb(cr, 1, 1, 1);
        cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL,
                               CAIRO_FONT_WEIGHT_NORMAL);
        cairo_set_font_size(cr, 20.0);

        cairo_move_to(cr, xpos, ypos);
        sprintf(buff, "Image: %s", s->image);
        cairo_show_text(cr, buff);
        cairo_stroke(cr);
        ypos += dy;

        if(s->dotfile != NULL)
        {
            cairo_move_to(cr, xpos, ypos);
            sprintf(buff, "Table: %s", s->dotfile);
            cairo_show_text(cr, buff);
            cairo_stroke(cr);
            ypos+=dy;
        }

        sprintf(buff, "Dots: %d/%zu shown", nshow, T->nrow);
        cairo_move_to(cr, xpos, ypos);
        cairo_show_text(cr, buff);
        cairo_stroke(cr);
        ypos += dy;

        sprintf(buff, "Dot threshold: >= %.1f", minvalue);
        if(s->dot_mode == DOT_MODE_AUTO)
        {
            strncat(buff, " (auto)", bufflen);
        }
        cairo_move_to(cr, xpos, ypos);
        cairo_show_text(cr, buff);
        cairo_stroke(cr);
        ypos += dy;

        sprintf(buff, "Dynamic range: [%.1f, %.1f]", low, high);
        cairo_move_to(cr, xpos, ypos);
        cairo_show_text(cr, buff);
        cairo_stroke(cr);
        ypos += dy;

        if(s->lsigma > 0)
        {
        sprintf(buff, "Gaussian filter, sigma=%.1f", s->lsigma);
        cairo_move_to(cr, xpos, ypos);
        cairo_show_text(cr, buff);
        cairo_stroke(cr);
        ypos += dy;
        }

        free(buff);
    }

    cairo_surface_destroy(surf);
    cairo_surface_destroy(mip);
    cairo_destroy(cr);
    free(M->V);
    free(M);
}




int dw_render(int argc, char ** argv)
{

    fim_tiff_init();
    opts * s = opts_new();

    argparsing(argc, argv, s);
    omp_set_num_threads(s->nthreads);
    if(s->verbose > 1)
    {
        opts_print(stdout, s);
    }


    char * inFile = s->image;
    char * outFile = s->outfile;

    if(!file_exist(inFile))
    {
        fprintf(stderr, "Can't open %s!\n", inFile);
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
    printf("Will write to %s\n", outFile);

    int64_t M = 0, N = 0, P = 0;
    float * A = fim_tiff_read(inFile, NULL, &M, &N, &P, s->verbose);
    fim_t * I = fim_image_from_array(A, M, N, P);
    free(A);

    ftab_t * T = NULL;
    if(s->dotfile != NULL)
    {
        T = ftab_from_tsv(s->dotfile);
    }

    {
        int valid_table = 1;
        if(ftab_get_col(T, "x") == -1)
        {
            fprintf(stderr, "No column named 'x'\n");
            valid_table = 0;
        }
        if(ftab_get_col(T, "y") == -1)
        {
            fprintf(stderr, "No column named 'x'\n");
            valid_table = 0;
        }
        if(ftab_get_col(T, "z") == -1)
        {
            fprintf(stderr, "No column named 'z'\n");
            valid_table = 0;
        }
        if(ftab_get_col(T, "value") == -1)
        {
            fprintf(stderr, "No column named 'value'\n");
            valid_table = 0;
        }
        if(valid_table == 0)
        {
            exit(EXIT_FAILURE);
        }
    }

    render(s, I, T);

    ftab_free(T);

    free(I->V);
    free(I);

    opts_free(s);
    return EXIT_SUCCESS;
}


#ifdef STANDALONE
int main(int argc, char ** argv)
{
    return dw_render(argc, argv);
}
#endif
