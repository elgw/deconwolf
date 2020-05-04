#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <tiffio.h>

int writetif(char * fName, float * V, 
    int M, int N, int P);

float * readtif_asFloat(char * fName, 
    int * M0, int * N0, int * P0);

void floatimage_normalize(float * restrict, const size_t);
