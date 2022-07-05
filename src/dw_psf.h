#ifndef __dw_psf_h__
#define __dw_psf_h__


#include <getopt.h>
#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>


#include "fim.h"
#include "ftab.h"
#include "fim_tiff.h"
#include "dw_version.h"
#include "dw_util.h"

/* For the PSF model see https://doi.org/10.1364/OL.28.000801
 * For the confocal model, see  https://doi.org/10.1111/j.1365-2818.1989.tb00577.x
*/

int dw_psf_cli(int argc, char ** argv);

#endif
