#include "dw_maxproj.h"

#define MODE_MAX 0
#define MODE_SLICE 1

typedef struct{
    int mode;
    int slice;
    int optpos; // Next argument not consumed by getargs
    int overwrite;
    int verbose;
} opts;

opts * opts_new()
{
    opts * s = malloc(sizeof(opts));

    s->mode = MODE_MAX;
    s->overwrite = 0;
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
    printf("Options:\n");
    printf(" --slice N\n\t Extract slice N\n");
    printf(" --overwrite\n\t Overwrite existing files\n");
    printf(" --verbose v\n\t Set verbosity level\n");
}


static void argparsing(int argc, char ** argv, opts * s)
{
    struct option longopts[] = {
                                 {"help", no_argument, NULL, 'h'},
                                 {"slice", required_argument, NULL, 's'},
                                 {"overwrite", no_argument, NULL, 'o'},
                                 {"verbose", required_argument, NULL, 'v'},
                                 {NULL, 0, NULL, 0}};
    int ch;
    while((ch = getopt_long(argc, argv, "ohs:v:", longopts, NULL)) != -1)
    {
        switch(ch){
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
        }
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


#ifdef STANDALONE
int main(int argc, char ** argv)
{
    return dw_tiff_max(argc, argv);
}
#endif


int dw_tiff_max(int argc, char ** argv)
{

    opts * s = opts_new();

  if(argc < 2)
  {
    usage(argc, argv);
    exit(1);
  }

  argparsing(argc, argv, s);

  fim_tiff_init();

  char * outFile;
  char * inFile;

  for(int ff = s->optpos; ff<argc; ff++)
  {
    inFile = argv[ff];

    if(!file_exist(inFile))
    {
        printf("Can't open %s!\n", inFile);
        exit(1);
    }

    outFile = malloc(strlen(inFile) + 20);
    if(s->mode == MODE_MAX)
    {
        char * _dname = strdup(inFile);
        char * dname = dirname(_dname);
        char * _fname = strdup(inFile);
        char * fname = basename(_fname);

        if(s->verbose > 1)
        {
            fprintf(stdout, "Input file: %s\n", inFile);
            fprintf(stdout, "'%s' /  '%s'\n", dname, fname);
        }

        sprintf(outFile, "%s/max_%s", dname, fname);
        free(_dname);
        free(_fname);
        if(s->verbose > 1)
        {
            fprintf(stdout, "Ouput file: %s\n", outFile);
        }

        if(s->overwrite == 0 && file_exist(outFile))
    {
      printf("%s exists, skipping.\n", outFile);
    } else {
      printf("%s -> %s\n", inFile, outFile);
      fim_tiff_maxproj(inFile, outFile);
    }}
    if(s->mode == MODE_SLICE)
    {
        sprintf(outFile, "s%04d_%s", s->slice, inFile);
        if(file_exist(outFile) && s->overwrite == 0)
        {
            printf("%s exists, skipping.\n", outFile);
        } else {
            printf("%s -> %s\n", inFile, outFile);
            fim_tiff_extract_slice(inFile, outFile, s->slice);
        }

    }
    free(outFile);
  }
  opts_free(s);
  return 0;
}
