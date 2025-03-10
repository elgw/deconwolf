#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

#include "fim_tiff.h"
#include "dw_util.h"

#include "dw_tiff_merge.h"

typedef struct {
    int overwrite;
    int verbose;
    int optind; /* Where the first non-consumed argument is stored */
    int test;
} tm_config_t;

tm_config_t * tm_config_new()
{
    tm_config_t * conf = malloc(sizeof(tm_config_t));
    assert(conf != NULL);
    conf->overwrite = 0;
    conf->verbose = 1;
    conf->test = 0;
    return conf;
}

void tm_config_free(tm_config_t * conf)
{
    free(conf);
}


void usage()
{
    printf("Usage for subcommand 'merge'\n");
    printf("dw merge [merge] output.tif input1.tif input2.tif ...\n");
    printf("Options:\n");
    printf("--help\n\t "
           "Show this help message and quit\n");
    printf("--verbose v\n\t "
           "Set verbosity level to v\n");
    printf("--overwrite\n\t "
           "Overwrite output image if it already exists\n");
    return;
}

static void argparsing(int argc, char ** argv, tm_config_t * conf)
{

    struct option longopts[] = {
        { "help",     no_argument,       NULL, 'h' },
        { "verbose",   required_argument, NULL, 'v' },
        { "overwrite",     no_argument,       NULL,  'o' },
        { "test",      no_argument, NULL, 't' },
        { NULL,           0,                 NULL,   0   }
    };

    int ch;
    while((ch = getopt_long(argc, argv,
                            "hv:ot",
                            longopts, NULL)) != -1)
    {
        switch(ch) {
        case 'h':
            usage();
            exit(EXIT_SUCCESS);
        case 'v':
            conf->verbose = atoi(optarg);
            break;
        case 'o':
            conf->overwrite = 1;
            break;
        case 't':
            conf->test = 1;
            break;
        default:
            fprintf(stderr, "Unknown argument\n");
            exit(EXIT_FAILURE);
        }
    }

    conf->optind = optind;
    return;
}

int dw_tiff_merge(int argc, char ** argv)
{
    tm_config_t * conf = tm_config_new();

    argparsing(argc, argv, conf);

    if(argc < 3)
    {
        usage();
        tm_config_free(conf);
        return EXIT_FAILURE;
    }

    char * outfile = argv[conf->optind];
    if(conf->overwrite == 0)
    {
        if(dw_isfile(outfile) == 1)
        {
            printf("%s already exists, doing nothing. Use --overwrite to overwrite existing files\n",
                   outfile);
            exit(EXIT_SUCCESS);
        }
    }

    fim_tiff_init();

    int64_t M = 0;
    int64_t N = 0;
    int64_t P = 0;

    int64_t M_out = 0;
    int64_t N_out = 0;
    int64_t P_out = 0;


    for(int kk = conf->optind+1 ; kk<argc; kk++)
    {

        if(kk > 1)
        {
            fim_tiff_get_size(argv[kk], &M, &N, &P);
            //printf(" %" PRId64 " x %" PRId64 " x %" PRId64 "", M, N, P);
            if(kk == optind+1)
            {
                M_out = M;
                N_out = N;
            }
            if(kk > optind+1)
            {
                if(M_out != M || N_out != N)
                {
                    fprintf(stderr, "Error while reading %s\n", argv[kk]);
                    fprintf(stderr, "Image sizes does not match. A previous image had size "
                            "%" PRIu64 " x %" PRIu64 " while %s has size %" PRIu64 " x %" PRIu64 "\n",
                            M_out, N_out, argv[kk], M, N);
                    return EXIT_FAILURE;
                }
            }
            P_out += P;
        }
    }

    printf("Output image size: %" PRIu64 " x %" PRIu64 " x %" PRIu64 "\n", M_out, N_out, P_out);

    float * im_out = malloc(M_out*N_out*P_out*sizeof(float));
    if(im_out == NULL)
    {
        fprintf(stderr, "Failed to allocate memory for the output image\n");
        exit(EXIT_FAILURE);
    }
    size_t plane = 0;
    for(int kk = conf->optind + 1; kk<argc; kk++)
    {
        printf("Reading %s\n", argv[kk]); fflush(stdout);
        if(!conf->test)
        {
            float * im = fim_tiff_read(argv[kk], NULL, &M, &N, &P, conf->verbose);
            assert(im != NULL);
            memcpy(im_out + plane*M_out*N_out,
                   im, M*N*P*sizeof(float));
            plane += P;
            free(im);
        }
    }
    printf("Writing to %s\n", outfile);
    if(! conf->test)
    {
        fim_tiff_write_float(outfile, im_out,
                             NULL, // tiff tags, TODO
                             M_out, N_out, P_out);
    }
    free(im_out);
    tm_config_free(conf);
    return EXIT_SUCCESS;
}
