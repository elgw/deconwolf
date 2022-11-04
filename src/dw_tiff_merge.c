#include "dw_tiff_merge.h"

void usage()
{
    printf("Usage\n");
    printf("dw merge output.tif input1.tif input2.tif ...\n");
    return;
}

int dw_tiff_merge(int argc, char ** argv)
{
    int verbose = 1;

    if(argc < 2)
    {
        return EXIT_FAILURE;
    }
    fim_tiff_init();

    int64_t M = 0;
    int64_t N = 0;
    int64_t P = 0;

    int64_t M_out = 0;
    int64_t N_out = 0;
    int64_t P_out = 0;

    fim_tiff_init();
    for(int kk = 0 ; kk<argc; kk++)
    {
        printf("argv[%d] = %s", kk, argv[kk]);
        if(kk > 1)
        {
            fim_tiff_get_size(argv[kk], &M, &N, &P);
            printf(" %lld x %lld x %lld", M, N, P);
            if(kk == 2)
            {
                M_out = M;
                N_out = N;
            }
            if(kk > 2)
            {
                if(M_out != M || N_out != N)
                {
                    fprintf(stderr, "Image sizes does not match. A previous image had size"
                            "%llu x %llu while %s has size %llu x %llu\n",
                            M_out, N_out, argv[kk], M, N);
                    return EXIT_FAILURE;
                }
            }
            P_out += P;
        }
        printf("\n");
    }

    printf("Output image size: %llu x %llu x %llu\n", M_out, N_out, P_out);

    float * im_out = malloc(M_out*N_out*P_out*sizeof(float));
    if(im_out == NULL)
    {
        fprintf(stderr, "Failed to allocate memory for the output image\n");
        return EXIT_FAILURE;
    }
    size_t plane = 0;
    for(int kk = 2; kk<argc; kk++)
    {
        ttags T;
        printf("Reading %s\n", argv[kk]); fflush(stdout);
        float * im = fim_tiff_read(argv[kk], &T, &M, &N, &P, verbose);
        assert(im != NULL);
        memcpy(im_out + plane*M_out*N_out,
               im, M*N*P*sizeof(float));
        plane += P;
        free(im);
    }
    printf("Writing to %s\n", argv[1]); fflush(stdout);
    fim_tiff_write_float(argv[1], im_out,
                         NULL, // tiff tags, TODO
                         M_out, N_out, P_out);
    free(im_out);
    return EXIT_SUCCESS;
}
