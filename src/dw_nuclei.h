#ifndef _dw_nuclei_h_
#define _dw_nuclei_h_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <libgen.h>

#include <png.h>
#include <zlib.h>

#include "fim.h"
#include "fim_tiff.h"
#include "dw_util.h"
#include "dw_version.h"
#include "random_forest/prf_forest.h"

int dw_nuclei(int argc, char ** argv);

#endif