/* (C) 2023 Erik L. G. Wernersson
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "dw_psf_sted.h"

typedef struct{
    int overwrite;
    int verbose;
    char * outfile;
    char * logfile;
    FILE * log;
    int nthreads;
    int M;
    int P;

    double Lfwhm_px;
    double Afwhm_px;

    double gamma;
    double sigma;
} opts;

static opts * opts_new();
static void opts_free(opts * s);
static void usage(__attribute__((unused)) int argc, char ** argv);
static void argparsing(int argc, char ** argv, opts * s);
static int file_exist(char * fname);


static double gauss1(double x, double s)
{
    return
        1.0 / (s * sqrt(2.0*M_PI))
        *exp(-0.5*pow(x/s, 2));
}

/* 2D Lorentzian centered at 0. */
static double lorentz2(double x, double y, double gamma)
{
    return
        1.0/(2.0*M_PI) *
        ( gamma /
          (
              pow(pow(x, 2) + pow(y, 2) + pow(gamma, 2), 3.0/2.0)
              )
            );
}

static double sigma_from_fwhm_gaussian(double fwhm)
{
    return fwhm / (2*sqrt(2*log(2)));
}

static double gamma_from_fwhm_lorentz2(double fwhm)
{
    return fwhm / (sqrt(pow(2, 2.0/3.0)-1)*2.0);
}

static opts * opts_new()
{
    opts * s = malloc(sizeof(opts));

    s->overwrite = 0;
    s->verbose = 1;
    s->outfile = NULL;
    s->logfile = NULL;
    s->log = NULL;
    s->nthreads = dw_get_threads();

    s->M = 81;
    s->P = 81;
    s->Lfwhm_px = 6;
    s->Afwhm_px = 5;

    s->sigma = 3;
    s->gamma = 3;
    return s;
}


static void opts_free(opts * s)
{
    dw_nullfree(s->outfile);
    dw_nullfree(s->logfile);
    if(s->log != NULL)
    {
        fclose(s->log);
    }
    free(s);
}

static void opts_print(FILE * f, opts * s)
{
    fprintf(f, "Deconwolf %s STED PSF generator.\n", deconwolf_version);
    fprintf(f, "overwrite = %d\n", s->overwrite);
    fprintf(f, "verbose = %d\n", s->verbose);
    fprintf(f, "nthreads = %d (ignored)\n", s->nthreads);
    fprintf(f, "outfile: %s\n", s->outfile);
    fprintf(f, "Out size: [%d x %d x %d] pixels\n", s->M, s->M, s->P);
    fprintf(f, "Lateral FWHM = %f px", s->Lfwhm_px);
    fprintf(f, " (gamma = %f)\n", s->gamma);
    fprintf(f, "Axial FWHM = %f px", s->Afwhm_px);
    fprintf(f, " (sigma = %f)\n", s->sigma);
    fprintf(f, "Lateral model: 2D Lorentzian\n");
    fprintf(f, "Axial model: Gaussian\n");
    const char * user = getenv("USER");
    if(user == NULL)
    {
        user = getenv("USERNAME");
    }

    if(user != NULL)
    {
        fprintf(f, "User: %s\n", getenv("USER"));
    }
    return;
}

static void usage(__attribute__((unused)) int argc, char ** argv)
{
    opts * s = opts_new();
    printf("Deconwolf %s STED PSF generator.\n", deconwolf_version);
    printf("usage: %s [<options>] output.tif\n", argv[0]);
    printf("Options:\n");
    printf(" --lateral l\n\t Lateral FWHM in pixels\n");
    printf(" --axial l\n\t Axial FWHM in pixels\n");
    printf(" --size N\n\t Set number of pixels along lateral dimensions\n");
    printf(" --nslice P\n\t Set number of pixels along axial dimension\n");
    printf("General:\n");
    printf(" --overwrite\n\t Overwrite existing files\n");
    printf(" --help\n\t Show this message\n");
    printf(" --verbose v\n\t Verbosity level\n");
    printf("\n");
    opts_free(s);
    return;
}

static void argparsing(int argc, char ** argv, opts * s)
{
    if(argc == 1)
    {
        printf("Deconwolf %s STED PSF generator.\n", deconwolf_version);
        printf("See `dw psf_sted --help` or `man dw`.\n");
        exit(EXIT_SUCCESS);
    }
    int hasL = 0;
    int hasA = 0;

    struct option longopts[] = {
        {"axial",     required_argument, NULL, 'A'},
        {"lateral",     required_argument, NULL, 'L'},
        {"help", no_argument, NULL, 'h'},
        {"overwrite", no_argument, NULL, 'o'},
        {"nslice", required_argument, NULL, 'p'},
        {"size", required_argument, NULL, 's'},
        {"threads", required_argument, NULL, 't'},
        {"verbose", required_argument, NULL, 'v'},
        {NULL, 0, NULL, 0}};
    int ch;
    while((ch = getopt_long(argc, argv,
                            "A:L:hop:s:t:b:", longopts, NULL)) != -1)
    {
        switch(ch){
        case 'A':
            s->Afwhm_px = atof(optarg);
            hasA = 1;
            break;
        case 'L':
            s->Lfwhm_px = atof(optarg);
            hasL = 1;
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
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
        case 'p':
            s->P = atoi(optarg);
            break;
        case 's':
            s->M = atoi(optarg);
            break;
        }
    }
    int ok = 1;
    if(hasL == 0)
    {
        fprintf(stderr, "Lateral FWHM [pixels] not specified (--lateral)\n");
        ok = 0;
    }
    if(hasA == 0)
    {
        fprintf(stderr, "Axial FWHM [pixels] not specified (--axial)\n");
        ok = 0;
    }

    s->gamma = gamma_from_fwhm_lorentz2(s->Lfwhm_px);
    s->sigma = sigma_from_fwhm_gaussian(s->Afwhm_px);


    if(ok == 0)
    {
        exit(EXIT_FAILURE);
    }

    /* We work with odd number of pixels */
    if(s->M % 2 == 0)
    {
        s->M++;
    }

    if(s->P % 2 == 0)
    {
        s->P++;
    }

    if(optind == argc)
    {
        s->outfile = malloc(1024);

        sprintf(s->outfile,
                "spsf_L%.1f_A%.1f.tif",
                s->Lfwhm_px, s->Afwhm_px);

        if(s->verbose > 0)
        {
            printf("Will write to %s\n", s->outfile);
        }
    } else {
        s->outfile = strdup(argv[optind]);
    }

    if(s->overwrite == 0)
    {
        if(file_exist(s->outfile))
        {
            printf("%s exists, leaving\n", s->outfile);
            printf("if you don't care, use --overwrite\n");
            exit(EXIT_FAILURE); /* TODO */
        }
    }

    s->logfile = malloc(strlen(s->outfile) + 32);
    sprintf(s->logfile, "%s.log.txt", s->outfile);
    s->log = fopen(s->logfile, "w");
    if(s->log == NULL)
    {
        fprintf(stderr, "Failed to open %s for writing\n", s->logfile);
        exit(EXIT_FAILURE);
    }

    for(int kk = 0; kk< argc; kk++)
    {
        fprintf(s->log, "%s", argv[kk]);
        if(kk+1 != argc)
        {
            fprintf(s->log, " ");
        } else {
            fprintf(s->log, "\n");
        }
    }

    opts_print(s->log, s);
    if(s->verbose > 1)
    {
        opts_print(stdout, s);
    }
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



static fim_t * gen_psf(opts * s)
{
    fim_t * PSF = fimt_zeros(s->M, s->M, s->P);
    const double gamma = s->gamma;
    const double sigma = s->sigma;

    assert(gamma > 0);
    assert(sigma > 0);

    int x0 = (s->M -1) / 2;
    const int M = s->M;

    for(int zz = 0; zz < s->P; zz++)
    {
        double zpos = zz-(s->P-1)/2;
        double Wz = gauss1(zpos, sigma);
        float * plane = PSF->V + zz*PSF->M*PSF->M;

        for(int xx = 0; xx < M; xx++)
        {
            for(int yy = 0; yy < M; yy++)
            {
                double x = xx-x0;
                double y = yy-x0;
                plane[yy + s->M*xx] = Wz*lorentz2(x, y, gamma);
            }
        }
    }

    return PSF;
}


static void unit_tests()
{
    /* Check the relation between FWHM and parameters */
    for(double fwhm = 1; fwhm < 10; fwhm+= 0.1)
    {
        double sigma = sigma_from_fwhm_gaussian(fwhm);
        double g0 = gauss1(0, sigma);
        double g1 = gauss1(fwhm/2.0, sigma);
        //printf("G %f %f %f\n", g0, g1, g0/g1);
        assert(fabs(g0/g1 -2) < 1e-5 );

        double gamma = gamma_from_fwhm_lorentz2(fwhm);
        double l0 = lorentz2(0, 0, gamma);
        double l1 = lorentz2(fwhm/2.0, 0, gamma);
        //printf("L %f %f %f\n", l0, l1, l0/l1);
        assert(fabs(l0/l1 -2) < 1e-5 );
    }
    return;
}

static void dw_psf_sted(opts * s)
{
    unit_tests();
    if(s->verbose > 0)
    {
        printf("Writing to %s\n", s->outfile);
    }
    fim_t * PSF = gen_psf(s);

    float sum = fimt_sum(PSF);
    fim_mult_scalar(PSF->V, fimt_nel(PSF), 1.0/sum);

    ttags * T = ttags_new();
    char * swstring = malloc(1024);
    sprintf(swstring, "deconwolf %s", deconwolf_version);
    ttags_set_software(T, swstring);
    ttags_set_imagesize(T, PSF->M, PSF->N, PSF->P);

    free(swstring);


    fim_tiff_write_float(s->outfile, PSF->V,
                         T,
                         PSF->M, PSF->N, PSF->P);
    ttags_free(&T);
    fim_free(PSF);
}

int dw_psf_sted_cli(int argc, char ** argv)
{
    /* Get default settings */
    opts * s = opts_new();
    /* Parse command line and open log file etc */
    argparsing(argc, argv, s);
    if(s->verbose > 2)
    {
        printf("Command lined parsed, continuing\n"); fflush(stdout);
    }
    /* Ready to go */
    dw_psf_sted(s);
    /* Clean up */
    opts_free(s);
    return 0;
}

#ifdef STANDALONE
int main(int argc, char ** argv)
{
    return dw_psf_sted_cli(argc, argv);
}
#endif
