#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>

#include "fim.h"
#include "fim_tiff.h"
#include "dw_version.h"
#include "dw_util.h"

#include "dw_maxproj.h"

typedef enum {
    MODE_MAX,
    MODE_MAX_XYZ,
    MODE_SLICE,
    MODE_GM
} projection_type;

typedef struct{
    projection_type mode;
    int slice;
    int optpos; // Next argument not consumed by getargs
    int overwrite;
    int verbose;
} opts;

opts * opts_new()
{
    opts * s = calloc(1, sizeof(opts));
    assert(s != NULL);
    s->mode = MODE_MAX;
    s->verbose = 1;
    return s;
}

void opts_free(opts * s)
{
    free(s);
}

static void usage(__attribute__((unused)) int argc, char ** argv)
{
    printf("usage: %s [<options>] input1.tif input2.tif ... \n", argv[0]);
    printf("Modes:\n");
    printf("  --xyz \n\t"
           "A collage of max projections along x, y and z shown in a single image.\n");
    printf("  --slice N\n\t"
           "Extract slice N\n");
    printf("  --gm\n\t"
           "Extract the slice with the highest gradient magnitude\n");
    printf("Options:\n");
    printf("  --overwrite\n\t"
           "Overwrite existing files\n");
    printf("  --verbose v\n\t"
           "Set verbosity level\n");
    return;
}


static void argparsing(int argc, char ** argv, opts * s)
{
    struct option longopts[] = {
        {"help", no_argument, NULL, 'h'},
        {"slice", required_argument, NULL, 's'},
        {"overwrite", no_argument, NULL, 'o'},
        {"verbose", required_argument, NULL, 'v'},
        {"xyz", no_argument, NULL, 'x'},
        {"gm", no_argument, NULL, 'g'},
        {NULL, 0, NULL, 0}};
    int ch;
    while((ch = getopt_long(argc, argv, "gohs:v:x", longopts, NULL)) != -1)
    {
        switch(ch){
        case 'g':
            s->mode = MODE_GM;
            break;
        case 'o':
            s->overwrite = 1;
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
            break;
        case 's':
            s->slice = atoi(optarg);
            s->mode = MODE_SLICE;
            break;
        case 'v':
            s->verbose = atoi(optarg);
            break;
        case 'x':
            s->mode = MODE_MAX_XYZ;
            break;
        default:
            fprintf(stderr, "Unknown command line option");
            exit(EXIT_FAILURE);
            break;
        }
    }
    s->optpos = optind;
    return;
}


#ifdef STANDALONE
int main(int argc, char ** argv)
{
    return dw_tiff_max(argc, argv);
}
#endif

static char * get_outfile_name_for_max(const char * inFile)
{

    char * dname = dw_dirname(inFile);
    char * fname = dw_basename(inFile);
    char * outFile = malloc(strlen(inFile) + 20);
    assert(outFile != NULL);
    if(strlen(dname) > 0)
    {
        sprintf(outFile, "%s%cmax_%s", dname, FILESEP, fname);
    } else {
        sprintf(outFile, "max_%s", fname);
    }

    free(dname);
    free(fname);

    return outFile;
}

static void gen_gm(opts * s, const char * inFile, const char * outFile)
{
    fimo * I = fimo_tiff_read(inFile);
    if(I == NULL)
    {
        printf("Error reading %s\n", inFile);
        return;
    }
    float * gm = fim_focus_gm(I, 3);
    int slice = float_arg_max(gm, I->P);
    if(s->verbose > 0)
    {
        printf("Using slice %d/%d (index %d)\n", slice+1,
               (int) I->P, slice);
    }
    fimo * Z = fimo_get_plane(I, slice);
    fimo_free(I);
    fimo_tiff_write(Z, outFile);
    fimo_free(Z);
    return;
}

int
dw_tiff_max(int argc, char ** argv)
{

    opts * s = opts_new();

    if(argc < 2)
    {
        usage(argc, argv);
        exit(1);
    }

    argparsing(argc, argv, s);

    ftif_t * ftif = fim_tiff_new(stdout, s->verbose);

    char * inFile;

    for(int ff = s->optpos; ff<argc; ff++)
    {
        inFile = argv[ff];
        if(s->verbose > 1)
        {
            printf("Infile: %s\n", argv[ff]);
        }


        if(!dw_isfile(inFile))
        {
            printf("Can't open %s!\n", inFile);
            exit(EXIT_FAILURE);
        }

        fim_tiff_info info = {0};
        fim_tiff_get_info(ftif, inFile, &info);
        if(info.P <= 1)
        {
            printf("%s is 2D, skipping\n", inFile);
            continue;
        }

        /* Output file name */
        char *  outFile = NULL;
        switch(s->mode)
        {
        case MODE_MAX:
            outFile = get_outfile_name_for_max(inFile);
            break;
        case MODE_MAX_XYZ:
            outFile = get_outfile_name_for_max(inFile);
            break;
        case MODE_GM:
            outFile = get_outfile_name_for_max(inFile);
            break;
        case MODE_SLICE:
            outFile = malloc(strlen(inFile) + 20);
            sprintf(outFile, "s%04d_%s", s->slice, inFile);
            break;
        default:
            fprintf(stderr, "Unknown mode!\n");
            exit(EXIT_FAILURE);
        }

        if(s->verbose > 1)
        {
            fprintf(stdout, "Input file: %s\n", inFile);
            fprintf(stdout, "Output file: %s\n", outFile);
        }

        if(s->overwrite == 0 && dw_isfile(outFile))
        {
            printf("%s exists, skipping.\n", outFile);
            free(outFile);
            continue;
        }

        if(s->verbose > 0)
        {
            printf("%s -> %s\n", inFile, outFile);
        }

        switch(s->mode)
        {
        case MODE_MAX:
            if(fim_tiff_maxproj(ftif, inFile, outFile))
            {
                exit(EXIT_FAILURE);
            }
            break;
        case MODE_MAX_XYZ:
            fim_tiff_maxproj_XYZ(ftif, inFile, outFile);
            break;
        case MODE_GM:
            gen_gm(s, inFile, outFile);
            break;
        case MODE_SLICE:
            fim_tiff_extract_slice(ftif, inFile, outFile, s->slice);
            break;
        default:
            fprintf(stderr, "Unknown mode!\n");
            exit(EXIT_FAILURE);
            break;
        }

        free(outFile);
    }

    fim_tiff_destroy(ftif);
    opts_free(s);
    return EXIT_SUCCESS;
}
