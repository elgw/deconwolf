#pragma once

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <libgen.h>

#include "fim_tiff.h"
#include "dw_util.h"

/* Create a tiff file argv[1] by merging the files in argv[2]
 * ... argv[argc-1] */

int dw_tiff_merge(int argc, char ** argv);
