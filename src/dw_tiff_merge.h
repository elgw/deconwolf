#ifndef __dw_tiff_merge_h__
#define __dw_tiff_merge_h__

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "fim_tiff.h"


/* Create a tiff file argv[1] by merging the files in argv[2]
 * ... argv[argc-1] */

int dw_tiff_merge(int argc, char ** argv);

#endif
