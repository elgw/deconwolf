#include <assert.h>
#include <math.h>
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

typedef struct {
    int verbose;
    int overwrite;
    char * outfile;
    int optpos;
    float sigma;
} opts;

static opts *
opts_new(void)
{
    opts * s = calloc(1, sizeof(opts));
    assert(s != NULL);
    s->sigma = 100;
    return s;
}

static void
opts_free(opts * s)
{
    free(s->outfile);
    free(s);
    return;
}

static void usage(const char * progname)
{
    opts * defaults = opts_new();
    printf(
           "%s can be used to estimate vignetting.\n"
           "\n"
           "It calculates the average image intensity over several images, which can\n"
           "be either 2D or 3D. In case of 3D images they are averaged\n"
           "along the Z-axis to produce a 2D image. After averaging a Gaussian low pass\n"
           "filter is applied.\n"
           "The output is always a 2D tif image which is not normalized in any way.\n"
           "The approach can work well even for images with biological content, given that\n"
           "the number of images is high.\n",
           progname);

    printf("\n");
    printf("Usage:\n");
    printf("%s [OPTIONS] --out bg.tif file1.tif file2.tif ...\n", progname);
    printf("\n");
    printf("Options:\n");
    printf("  --overwrite\n\t"
           "overwrite existing files\n");
    printf("  --verbose n\n\t"
           "set verbosity level to n\n");
    printf("  --sigma s\n\t"
           "set the sigma value of the Gaussian smoothing kernel\n\t"
           "default value: %.1f [pixels]\n", defaults->sigma);
    opts_free(defaults);
    return;
}

static int
argparsing(int argc, char ** argv, opts * s)
{
    struct option longopts[] = {
        {"help", no_argument, NULL, 'h'},
        {"overwrite", no_argument, NULL, 'o'},
        {"verbose", required_argument, NULL, 'v'},
        {"sigma", required_argument, NULL, 's'},
        {"out", required_argument, NULL, 't'},
        {NULL, 0, NULL, 0}};
    int ch;
    while((ch = getopt_long(argc, argv, "hov:t:s:", longopts, NULL)) != -1)
    {
        switch(ch){
        case 'h':
            usage(argv[0]);
            exit(EXIT_SUCCESS);
            break;
        case 'o':
            s->overwrite = 1;
            break;
        case 'v':
            s->verbose = atoi(optarg);
            break;
        case 's':
            s->sigma = atof(optarg);
            break;
        case 't':
            free(s->outfile);
            s->outfile = strdup(optarg);
            break;
        }
    }
    s->optpos = optind;

    if(s->outfile == NULL)
    {
        printf("Please specify the output file with --out\n");
        return 1;
    }

    return 0;


}

int
dw_background(int argc, char ** argv)
{
    opts * s = opts_new();
    if(argparsing(argc, argv, s))
    {
        opts_free(s);
        exit(EXIT_FAILURE);
    }

    if(s->optpos == argc)
    {
        printf("No input files given\n");
        opts_free(s);
        exit(EXIT_FAILURE);
    }

    printf("Reading %s\n", argv[s->optpos]);
    fimo * I = fimo_tiff_read(argv[s->optpos]);
    fimo * bg = fimo_sumproj(I);
    fimo_free(I);

    for(int kk = s->optpos+1; kk < argc; kk++)
    {
        printf("Reading %s\n", argv[kk]);
        fimo * I = fimo_tiff_read(argv[kk]);
        if(I->M != bg->M || I->N != bg->M)
        {
            printf("%s dimensions %lu x %lu does not match previous: %lu x %lu\n"
                   "   ignoring this file\n", argv[kk], I->M, I->N, bg->M, bg->N);
            fimo_free(I);
            continue;
        }
        fimo * P = fimo_sumproj(I);
        fimo_free(I);
        fimo_add(bg, P);
        fimo_free(P);
    }

    printf("Gaussian filter, sigma = %.1f\n", s->sigma);
    fimo_gsmooth(bg, s->sigma);


    printf("Writing to %s\n", s->outfile);
    fimo_tiff_write(bg, s->outfile);
    fimo_free(bg);

    opts_free(s);
    return EXIT_SUCCESS;
}

#ifdef STANDALONE
int main(int argc, char ** argv)
{
    return dw_background(argc, argv);
}
#endif
