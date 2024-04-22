#include <getopt.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include "fim_tiff.h"
#include "dw_util.h"
#include "sparse_preprocess.h"



typedef struct {
    float lambda;
    float lambda_s;
    int iter;
    int threads;
    int optpos;
    char * prefix;
    int overwrite;
    int verbose;
    int periodic;
    int directions;
} sparse_settings_t;

static void sparse_settings_write(sparse_settings_t * s, FILE * fid)
{
    fprintf(fid, "lambda=%e\n", s->lambda);
    fprintf(fid, "lambda_s=%e\n", s->lambda_s);
    fprintf(fid, "iter=%d\n", s->iter);
    fprintf(fid, "threads=%d\n", s->threads);
    fprintf(fid, "periodic=%d\n", s->periodic);
    fprintf(fid, "directions=%d\n", s->directions);
}

static sparse_settings_t * sparse_settings_new()
{
    sparse_settings_t * conf = calloc(1, sizeof(sparse_settings_t));
    assert(conf != NULL);
    conf->threads = -1;
    conf->iter = 60;
    conf->prefix = strdup("nf");
    conf->overwrite = 0;
    conf->verbose = 1;
    conf->lambda = 0.002;
    conf->lambda_s = 0.004;
    conf->periodic = 0;
    conf->directions = 3;
    return conf;
}

static void
sparse_settings_free(sparse_settings_t * conf)
{
    free(conf->prefix);
    free(conf);
}

static void usage(int argc, char ** argv)
{
    printf("For noise removal in images (before deconvolution)\n");
    printf("\n");
    printf("Usage:\n");
    sparse_settings_t * s = sparse_settings_new();

    if(argc > 0)
    {
        printf("%s [options] image.tif\n", argv[0]);
    }
    printf("\n");

    printf("Options:\n");
    printf("--lambda lambda\n\t"
           "set the lambda parameter (data fidelity, default=%e)\n", s->lambda);
    printf("--lambda_s lambda_s\n\t"
           "set the lambda_s parameter (sparsity, default=%e)\n", s->lambda_s);
    printf("--iter iter\n\t"
           "Number of iterations (default=%d)\n", s->iter);
    printf("--overwrite\n\t"
           "Overwrite existing output files\n");
    printf("--prefix p\n\t"
           "set the prefix of the output file (default=%s)\n", s->prefix);
    printf("--threads t\n\t"
           "Number of threads (default=auto)\n");
    printf("--verbose v\n\t"
           "Set the verbosity level (default=%d)\n", s->verbose);
    printf("--periodic\n\t"
           "Treat the image as periodic\n");
    printf("--nine\n\t"
           "use 9, not 3 directions for the Hessian\n");
    sparse_settings_free(s);
    return;
}


static sparse_settings_t * get_cli_options(int argc, char ** argv)
{

    sparse_settings_t * conf = sparse_settings_new();


    struct option longopts[] = {
        {"lambda", required_argument, NULL, '1'},
        {"lambda_s", required_argument, NULL, '2'},
        {"help", no_argument, NULL, 'h'},
        {"iter", required_argument, NULL, 'I'},
        {"nine", no_argument, NULL, 'n'},
        {"overwrite", no_argument, NULL, 'o'},
        {"periodic", no_argument, NULL, 'P'},
        {"prefix", required_argument, NULL, 'p'},
        {"threads", required_argument, NULL, 'T'},
        {"verbose", required_argument, NULL, 'v'},
        {NULL, 0, NULL, 0}};


    int opt;
    while( (opt = getopt_long(argc, argv, "1:2:hI:noPp:T:v:", longopts, NULL)) != -1)
    {
        switch(opt)
        {
        case '1':
            conf->lambda = atof(optarg);
            break;
        case '2':
            conf->lambda_s = atof(optarg);
            break;
        case 'h':
            usage(argc, argv);
            sparse_settings_free(conf);
            exit(EXIT_SUCCESS);
            break;
        case 'I':
            conf->iter = atoi(optarg);
            break;
        case 'n':
            conf->directions = 9;
            break;
        case 'o':
            conf->overwrite = 1;
            break;
        case 'P':
            conf->periodic = 1;
            break;
        case 'p':
            free(conf->prefix);
            conf->prefix = strdup(optarg);
            break;
        case 'T':
            conf->threads = atoi(optarg);
            break;
        case 'v':
            conf->verbose = atoi(optarg);
            break;
        }
    }
    if(conf->prefix == NULL)
    {
        fprintf(stderr, "Error: Empty prefix name\n");
        goto fail;
    }

    if(strlen(conf->prefix) == 0)
    {
        fprintf(stderr, "Error: The prefix must have a length > 0\n");
        goto fail;
    }

    if(conf->threads < 1)
    {
        conf->threads = dw_get_threads();
    }

    omp_set_num_threads(conf->threads);
    conf->optpos = optind;

    return conf;

 fail:
    sparse_settings_free(conf);
    return NULL;
}

int sparse_preprocess_cli(int argc, char ** argv)
{

    sparse_settings_t * conf = get_cli_options(argc, argv);
    if(conf == NULL)
    {
        goto fail;
    }

    fim_tiff_init();

    char * inFile;

    for(int ff = conf->optpos; ff<argc; ff++)
    {
        inFile = argv[ff];

        if(!dw_file_exist(inFile))
        {
            printf("Can't open %s!\n", inFile);
            continue;
        }

        assert(inFile != NULL);
        assert(conf->prefix != NULL);
        char * outfile = dw_prefix_file(inFile, conf->prefix);

        if((conf->overwrite == 0) && dw_file_exist(outfile))
        {
            printf("%s already exist, skipping\n", outfile);
            free(outfile);
            continue;
        }
        if(conf->verbose > 0)
        {
            printf("Processing %s -> %s(.log.txt)\n", inFile, outfile);
        }

        if(conf->verbose > 1)
        {
            printf("Reading %s\n", inFile);
        }
        ttags T;
        int64_t M, N, P;
        float * I = fim_tiff_read(inFile,
                                  &T,
                                  &M, &N, &P,
                                  conf->verbose);
        if(I == NULL)
        {
            fprintf(stderr, "Can't read %s. Does not look like a tif file\n", inFile);
            continue;
        }

        if(conf->verbose > 1)
        {
            printf("Opening log file\n");
        }
        /* TODO: Also dw version, and iterations status before closing */
        char * logfile = calloc(strlen(outfile) + 32, 1);
        assert(logfile != NULL);
        sprintf(logfile, "%s.log.txt", outfile);

        if(conf->verbose > 1)
        {
            printf("logfile = %s\n", logfile);
        }
        FILE * log = fopen(logfile, "w");
        if(log == NULL)
        {
            fprintf(stderr, "Unable to open %s\n", logfile);
            free(I);
            continue;
        }
        fprintf(log, "deconwolf: '%s'\n", deconwolf_version);
        fprintf(log, "CMD: ");
        for(int kk = 0; kk<argc; kk++)
        {
            fprintf(log, "%s ", argv[kk]);
        }
        fprintf(log, "\n");

        sparse_settings_write(conf, log);

        free(logfile);


        if(conf->verbose > 1)
        {
            printf("Iterative smoothing\n");
        }

        float * J = sparse_preprocess(I,
                                      M, N, P,
                                      conf->lambda, conf->lambda_s,
                                      conf->periodic, conf->directions,
                                      conf->iter, conf->verbose, log);

        free(I);
        fim_tiff_write_float(outfile, J,
                       &T, M, N, P);
        free(J);

        free(outfile);
        fclose(log);
    }

    if(conf->verbose > 1)
    {
        printf("Peak memory: %lu KB\n", get_peakMemoryKB());
    }

    sparse_settings_free(conf);

    return EXIT_SUCCESS;

 fail:
    return EXIT_FAILURE;
}

#ifdef STANDALONE
int main(int argc, char ** argv)
{
    return sparse_preprocess_cli(argc, argv);
}
#endif
