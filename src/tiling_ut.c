#include <stdio.h>
#include <stdlib.h>
#include "tiling.h"
#include "tiffio.h"

int main(int argc, char ** argv)
{

  int M = 1024, N = 1024, P = 1;
  int overlap = 2;
  int maxSize = 400;

  if(argc == 6)
  {
    M = atol(argv[1]);
    N = atol(argv[2]);
    P = atol(argv[3]);
    maxSize = atol(argv[4]);
    overlap = atol(argv[5]);
  }

  tiling * T = tiling_create(M,N,P, maxSize, overlap);
  tiling_show(T);

  printf("Generating weight map\n"); fflush(stdout);
  float * W = malloc(M*N*sizeof(float));
  for(int aa = 0; aa<M; aa++)
  {
    for(int bb = 0; bb<N; bb++)
    {
      float w = tiling_getWeights(T, aa, bb, 0);
      assert(w>0);
      W[aa + bb*M] = w;
    }
  }

  char outFile[] = "tiling_ut_weights.tif";
  printf("Writing weights to %s\n", outFile); fflush(stdout);
  writetif(outFile, W, M, N, P);
  free(W);

  tiling_free(T);
  free(T);


}
