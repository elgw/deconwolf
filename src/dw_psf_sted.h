#pragma once

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

/* Provides a PSF model for STED images taking FWHM as input
 * parameter */

/* Command line interface */
int dw_psf_sted_cli(int argc, char ** argv);
