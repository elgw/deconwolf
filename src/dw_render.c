#include "dw_render.h"

typedef struct{
    int overwrite;
    int verbose;
    int optpos;
    char * image; /* Image to analyze */
    char * outfile; /* Where to write the tsv output */
    char * dotfile;
    int nthreads;
    float asigma; /* Axial sigma */
    float lsigma; /* Lateral sigma */
    int ndots; /* Number of dots to export */
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

static double clockdiff(struct timespec* end, struct timespec * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}



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

    s->lsigma = 1;
    s->asigma = 1;
    s->ndots = -1; /* Auto */

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
    nullfree(s->outfile);
    nullfree(s->image);

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
    fprintf(f, "Axial sigma: %.2f\n", s->asigma);

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
    printf("usage: %s [<options>] --image input.tif\n", argv[0]);
    printf("Options:\n");
    printf(" --lsigma s\n\t Lateral sigma\n");
    printf(" --asigma s\n\t Axial sigma\n");
    printf(" --overwrite\n\t Overwrite existing files\n");
    printf(" --help\n\t Show this message\n");
    printf(" --ndots n\n\t Number of dots to export\n");
    printf(" --dots file.tsv\n\t Specify TSV file to read\n");
    printf("\n");

}

static void argparsing(int argc, char ** argv, opts * s)
{
    struct option longopts[] = {

        {"out", required_argument, NULL, 'O'},
        {"asigma", required_argument, NULL, 'a'},
        {"dots",   required_argument, NULL, 'd'},
        {"help", no_argument, NULL, 'h'},
        {"image", required_argument, NULL, 'i'},
        {"lsigma", required_argument, NULL, 'l'},
        {"ndots",   required_argument, NULL, 'n'},
        {"overwrite", no_argument, NULL, 'o'},
        {"threads", required_argument, NULL, 't'},
        {"verbose", required_argument, NULL, 'v'},
        {NULL, 0, NULL, 0}};
    int ch;
    while((ch = getopt_long(argc, argv, "a:d:hi:l:n:op:r:v:", longopts, NULL)) != -1)
    {
        switch(ch){

        case 'O':
            nullfree(s->outfile);
            s->outfile = strdup(optarg);
            break;
        case 'a':
            s->asigma = atof(optarg);
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
            break;
        case 'o':
            s->overwrite = 1;
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



cairo_surface_t * fim_image_to_cairo_surface(fim_image_t * I)
{
    cairo_surface_t * result;
    unsigned char * current_row;
    int stride;

    float imax = fim_max(I->V, I->M*I->N);

    result = cairo_image_surface_create(CAIRO_FORMAT_RGB24, I->M, I->N);
    if (cairo_surface_status(result) != CAIRO_STATUS_SUCCESS)
        return result;

    cairo_surface_flush(result);
    current_row = cairo_image_surface_get_data(result);
    stride = cairo_image_surface_get_stride(result);
    for (int y = 0; y < I->M; y++) {
        uint32_t *row = (void *) current_row;
        for (int x = 0; x < I->N; x++) {
            float v =  1.5*(I->V[x + I->M*y]) / imax*255.0;
            v>255? v = 255 : 0;
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

void render(opts * s, fim_image_t * I, fim_table_t * T)
{
    printf("Creating max projection\n");
    float * M = fim_maxproj(I->V, I->M, I->N, I->P);
    I->V = M; I->P = 1;
    cairo_surface_t * mip = fim_image_to_cairo_surface(I);

    cairo_surface_t * surf = cairo_svg_surface_create (s->outfile,
                                                        I->M,
                                                         I->N);
    //cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_NEAREST);
//    cairo_pattern_set_filter(mip, CAIRO_FILTER_NEAREST);

    cairo_t * cr = cairo_create(surf);

    cairo_set_source_surface(cr, mip, 0, 0);
    /* This does not do anything :( */
    cairo_pattern_set_filter (cairo_get_source (cr), CAIRO_FILTER_NEAREST);
    cairo_paint(cr);

    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_select_font_face(cr, "Sans", CAIRO_FONT_SLANT_NORMAL,
                           CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 20.0);

    cairo_move_to(cr, 10.0, 50.0);
    cairo_show_text(cr, s->image);
    cairo_stroke(cr);

    cairo_set_line_width(cr, 1);
    cairo_set_source_rgb(cr, 0.0, 1.00, 0);


    if(T != NULL)
    {
        printf("Drawing dots\n");
        for(size_t kk = 0; kk<T->nrow; kk++)
        {
            //printf("%f %f\n", T->T[kk*T->ncol], T->T[kk*T->ncol+1]);
            float x = T->T[kk*T->ncol];
            float y = T->T[kk*T->ncol + 1];
            //cairo_translate(cr, x, y);
            cairo_arc(cr, x, y, 5, 0, 2 * M_PI);
            cairo_close_path(cr);
            cairo_stroke(cr);
        }
    }
    printf("Drew %zu dots\n", T->nrow);

    cairo_surface_destroy(surf);
    cairo_surface_destroy(mip);
    cairo_destroy(cr);

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


    int64_t M = 0, N = 0, P = 0;
    float * A = fim_tiff_read(inFile, NULL, &M, &N, &P, s->verbose);
    fim_image_t * I = fim_image_from_array(A, M, N, P);

    fim_table_t * T = NULL;
    if(s->dotfile != NULL)
    {
        T = fim_table_from_tsv(s->dotfile);
    }
    render(s, I, T);
    fim_table_free(T);

    free(I);
    free(A);

    opts_free(s);
    return EXIT_SUCCESS;
}


#ifdef STANDALONE
int main(int argc, char ** argv)
{
    return dw_render(argc, argv);
}
#endif
