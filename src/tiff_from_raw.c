#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include "fim_tiff.h"

static void usage(char ** argv)
{
    printf("Usage:\n");
    printf("$ %s input.raw M N P\n", argv[0]);
    printf("Will create input.raw.tiff\n");
    printf("Assumes that the raw data is of type float32_t\n");
    printf("The output image will be of type uint16_t\n");
}


static void check_nbytes(char * fname, size_t N)
{

    struct stat st;
    int status = stat(fname, &st);
    if(status != 0)
    {
        fprintf(stderr, "Failed to get file size of %s\n", fname);
        exit(1);
    }

    size_t nbytes = st.st_size;

    if(nbytes != N)
    {
        fprintf(stderr, "File has %zu bytes, expected %zu\n", nbytes, N);
        exit(1);
    }

    return;
}

int main(int argc, char ** argv)
{
    if(argc < 5)
    {
        usage(argv);
        exit(EXIT_FAILURE);
    }

    char * input = argv[1];
    char * output = malloc(strlen(input)+10);
    assert(output != NULL);
    sprintf(output, "%s.tif", input);
    size_t M = atol(argv[2]);
    size_t N = atol(argv[3]);
    size_t P = atol(argv[4]);
    size_t nel = M*N*P;


    check_nbytes(input, nel*sizeof(float));

    return fim_tiff_from_raw(output, M, N, P, input);
}
