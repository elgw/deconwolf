/*    Copyright (C) 2020 Erik L. G. Wernersson
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

/* Command line interface for deconwolf */

/* Extra modules can be enabled by un-commenting in
 * the header file. */
#include "dw.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>


#include "dw_util.h"
#include "dw_maxproj.h"
#include "dw_tiff_merge.h"
#include "dw_imshift.h"
#include "dw_version.h"
#ifndef WINDOWS
#include "dw_nuclei.h"
#endif
#include "dw_background.h"

/* Uncomment to include, requires linking with libpng and libz
 * can be build separately by the makefile in the src folder
 */
// #include "dw_otsu.h"
#ifndef WINDOWS
#include "dw_dots.h"
#include "dw_psf.h"
#include "dw_psf_sted.h"
#include "dw_align_dots.h"
#endif

#include "fim.h"
#include "fim_tiff.h"
#include "fft.h"
#include "tiling.h"
#include "sparse_preprocess_cli.h"

static int
deconwolf_cli_help(int argc, char ** argv)
{
    printf("usage: %s [--version] <command> [<options>]\n", argv[0]);
    printf("\n");
    printf("main commands\n");
    printf("   deconvolve   deconvolve 3D images\n");
    printf("   maxproj      Z-projections of 3D images\n");
#ifdef dw_module_dots
    printf("   dots         detect dots with sub pixel precision\n");
#endif
    printf("\n");
    printf("also available\n");
    printf("   merge        merge individual slices to volume\n");
#ifdef dw_module_psf
    printf("   psf          generate PSFs for Widefield and Confocal\n");
#endif
#ifdef dw_module_psf_sted
    printf("   psf-STED     PSFs for 3D STED\n");
#endif
#ifdef dw_module_nuclei
    printf("   nuclei       pixel classifier\n");
#endif
#ifdef dw_module_background
    printf("   background   vignetting/background estimation\n");
#endif
    printf("   align-dots   Estimate alignment between dot\n");
    printf("   imshift      Shift/translate tif images\n");
    printf("   tif2npy      convert a tif file to a Numpy .npy file\n");
    printf("   npy2tif      convert a Numpy .npy file to a tif file\n");
    printf("\n");
    printf("Web page: https://www.github.com/elgw/deconwolf/\n");
    return EXIT_SUCCESS;
}

static int
deconwolf_cli_version(int argc, char ** argv)
{
    printf("deconwolf version %s\n", deconwolf_version);
    return EXIT_SUCCESS;
}

static int
npy2tif(int argc, char ** argv)
{
    if(argc > 1 && strcmp(argv[1], "--help") == 0)
    {
        printf("Usage:\n");
        printf("%s input.npy output.tif\n", argv[0]);
        exit(EXIT_SUCCESS);
    }
    if(argc < 3)
    {
        printf("Usage:\n");
        printf("%s input.npy output.tif\n", argv[0]);
        return EXIT_FAILURE;
    }

    npio_t * npy = npio_load(argv[1]);
    if(npy == NULL)
    {
        printf("Failed to open %s as a npy file\n", argv[1]);
        return EXIT_FAILURE;
    }

    int64_t M = 1;
    int64_t N = 1;
    int64_t P = 1;

    if(npy->ndim > 3 || npy->ndim < 2)
    {
        printf("This function does only work with 2D and 3D arrays\n");
        npio_free(npy);
        return EXIT_FAILURE;
    }

    switch(npy->ndim)
    {
    case 3:
        M = npy->shape[2];
        N = npy->shape[1];
        P = npy->shape[0];
        break;
    case 2:
        M = npy->shape[1];
        N = npy->shape[0];
        P = 1;
        break;
    default:
        fprintf(stderr, "Unsupported number of dimensions: %d\n", npy->ndim);
        npio_free(npy);
        return EXIT_FAILURE;
    }


    size_t nel = M*N*P;
    float * V = fim_malloc(nel*sizeof(float));

    if(npy->dtype == NPIO_F32)
    {
        memcpy(V, npy->data, nel*sizeof(float));
        goto success;
    }

    if(npy->dtype == NPIO_U8)
    {
        uint8_t * IN = (uint8_t * ) npy->data;
        for(size_t kk = 0; kk < nel; kk++)
        {
            V[kk] = (float) IN[kk];
        }
        goto success;
    }

    if(npy->dtype == NPIO_U16)
    {
        uint16_t * IN = (uint16_t * ) npy->data;
        for(size_t kk = 0; kk < nel; kk++)
        {
            V[kk] = (float) IN[kk];
        }
        goto success;
    }

    if(npy->dtype == NPIO_I32)
    {
        int32_t * IN = (int32_t * ) npy->data;
        for(size_t kk = 0; kk < nel; kk++)
        {
            V[kk] = (float) IN[kk];
        }
        goto success;
    }

    printf("Unable to convert the following npy file to float\n");
    npio_print(stdout, npy);
    npio_free(npy);
    fim_free(V);
    return EXIT_FAILURE;

success:
    ftif_t * ftif = fim_tiff_new(stdout, 1);
    fim_tiff_write_float(ftif,
                         argv[2], V, NULL,
                         M, N, P);
    fim_tiff_destroy(ftif);
    ftif = NULL;

    fim_free(V);
    npio_free(npy);
    return EXIT_SUCCESS;
}


static int
tif2npy(int argc, char ** argv)
{
    if(argc > 1 && strcmp(argv[1], "--help") == 0)
    {
        printf("Usage:\n");
        printf("%s input.tif output.npy\n", argv[0]);
        exit(EXIT_SUCCESS);
    }
    if(argc < 3)
    {
        printf("Usage:\n");
        printf("%s input.tif output.npy\n", argv[0]);
        return EXIT_FAILURE;
    }
    ftif_t * ftif = fim_tiff_new(stdout, 1);
    int64_t M, N, P;
    float * I = fim_tiff_read(ftif, argv[1], NULL,
                              &M, &N, &P);
    fim_tiff_destroy(ftif);
    ftif = NULL;
    if(I == NULL)
    {
        printf("Failed to open %s as a tif file\n", argv[1]);
        return EXIT_FAILURE;
    }
    int shape[3] = {P, N, M};
    if(npio_write(argv[2], 3, shape, I, NPIO_F32, NPIO_F32) <= 0)
    {
        printf("Failed to write to %s\n", argv[2]);
        fim_free(I);
        return EXIT_FAILURE;
    }

    fim_free(I);
    return EXIT_SUCCESS;
}


int main(int argc, char ** argv)
{
    /* Most of the commands can be compiled into separate programs.

       If the command is not recognized, the 'deconvolve' command will
       be assumed.
    */

    if(argc <= 1)
    {
        return deconwolf_cli_help(argc, argv);
    }

    assert(argc > 1);

    if(strcmp(argv[1], "--version") == 0)
    {
        return deconwolf_cli_version(argc, argv);
    }

    if(strcmp(argv[1], "deconvolve") == 0)
    {
        return dw_deconvolve_cli(argc-1, argv+1);
    }
    if(strcmp(argv[1], "maxproj") == 0)
    {
        return dw_tiff_max(argc-1, argv+1);
    }
    if(strcmp(argv[1], "imshift") == 0)
    {
        return dw_imshift(argc-1, argv+1);
    }
    if(strcmp(argv[1], "merge") == 0)
    {
        return dw_tiff_merge(argc-1, argv+1);
    }

    if(strcmp(argv[1], "nuclei") == 0)
    {
#ifdef dw_module_nuclei
        return dw_nuclei(argc-1, argv+1);
#else
        fprintf(stderr, "dw was built without the 'nuclei' module\n");
        exit(EXIT_FAILURE);
#endif
    }

    if(strcmp(argv[1], "dots") == 0)
    {
#ifdef dw_module_dots
        return dw_dots(argc-1, argv+1);
#else
        fprintf(stderr, "dw was built without the 'dots' module\n");
        exit(EXIT_FAILURE);
#endif
    }

    if(strcmp(argv[1], "psf") == 0)
    {
#ifdef dw_module_psf
        return dw_psf_cli(argc-1, argv+1);
#else
        fprintf(stderr, "dw was built without the 'psf' module\n");
        exit(EXIT_FAILURE);
#endif
    }

    if(strcmp(argv[1], "psf-STED") == 0)
    {
#ifdef dw_module_psf_sted
        return dw_psf_sted_cli(argc-1, argv+1);
#else
        fprintf(stderr, "dw was built without the 'psf-STED' module\n");
        exit(EXIT_FAILURE);
#endif
    }

    if(strcmp(argv[1], "noise1") == 0)
    {
        return sparse_preprocess_cli(argc-1, argv+1);
    }

    if( strcmp(argv[1], "align-dots") == 0)
    {
#ifdef dw_module_align_dots
        return dw_align_dots(argc-1, argv+1);
#else
        fprintf(stderr, "dw was not built with the 'align-dots' module\n");
        exit(EXIT_FAILURE);
#endif
    }

    if( strcmp(argv[1], "background") == 0)
    {
#ifdef dw_module_background
        return dw_background(argc-1, argv+1);
#else
        fprintf(stderr, "dw was not built with the 'background' module\n");
        exit(EXIT_FAILURE);
#endif
    }

    if( strcmp(argv[1], "tif2npy") == 0)
    {
        return tif2npy(argc-1, argv+1);
    }

    if( strcmp(argv[1], "npy2tif") == 0)
    {
        return npy2tif(argc-1, argv+1);
    }

    // Fallback to 'deconvolve' if no command was specified
    return dw_deconvolve_cli(argc, argv);
}
