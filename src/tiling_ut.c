#include <stdio.h>
#include <stdlib.h>
#include "tiling.h"
#include "tiffio.h"

int main(int argc, char ** argv)
{

  int M = 1024, N = 1024, P = 1;
  int overlap = 2;
  int maxSize = 400;

  tiling * T = tiling_create(M,N,P, maxSize, overlap);
  tiling_show(T);

  float w = tiling_getWeights(T, 512, 512, 0);
  printf("w = %f\n", w);

  float * W = malloc(M*N*sizeof(float));
  for(int aa = 0; aa<M; aa++)
  {
    for(int bb = 0; bb<N; bb++)
    {
      W[aa + bb*M] = tiling_getWeights(T, aa, bb, 0);
    }
  }

  char outFile[] = "tiling_ut_weights.tif";
  printf("Writing weights to %s\n", outFile);
  writetif(outFile, W, M, N, P);
  free(W);

  tiling_free(T);
  free(T);


}
