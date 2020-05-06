#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <tiffio.h>

int writetif(char * fName, float * V, 
    int M, int N, int P);

// Read a 3D tif stack as a float array
float * readtif_asFloat(char * fName, 
    int * M0, int * N0, int * P0, int verbosity);

void floatimage_normalize(float * restrict, const size_t);
void floatimage_show_stats(float * I, size_t N, size_t M, size_t P);

