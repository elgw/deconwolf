#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tiling.h"
#include "fim.h"
#include "fim_tiff.h"

void test_getWeight1d()
{
  float a = 1, b = 6, c = 8, d = 12;
  printf("test_getWeight1d\n");
  printf("[%.0f %.0f %.0f %.0f]\n", a, b, c, d);
  for(int kk = 0; kk<14; kk++)
  {
    printf("%d %.2f\n", kk, getWeight1d(a, b, c, d, kk));
  }
  printf("\n");
}


void test_weights(int M, int N, int P, int maxSize, int overlap)
  {
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
  writetif(outFile, W, M, N, 1);

  // Assuming that is it same behaviour for all P
 for(size_t kk = 0; kk<M*N; kk++)
  {
    assert(W[kk] > 0);
  }

  free(W);

  tiling_free(T);
  free(T);
  }

void test_copy_paste(int M, int N, int P, int maxSize, int overlap)
  /* Extract tiles and put them into a new image, one by one */
  {
    printf("-> test_copy_paste\n"); fflush(stdout);

  tiling * T = tiling_create(M,N,P, maxSize, overlap);
//  tiling_show(T);

  float * source = malloc(M*N*P*sizeof(float));
  for(size_t kk = 0; kk<M*N*P; kk++)
  { 
    source[kk] = 7; //rand()/RAND_MAX; 
  }

  float * target = malloc(M*N*P*sizeof(float));
  char * fname = malloc(100);

  for(int tt = 0; tt<T->nTiles; tt++)
  {
    float * Ttile = tiling_get_tile(T, tt, source);
    tile * tl = T->tiles[tt];
    sprintf(fname, "tile%d.tif", tt);
    writetif(fname, Ttile, tl->xsize[0], tl->xsize[1], tl->xsize[2]);

    // P has size T->tiles[kk]->xsize;
    tiling_put_tile(T, tt, target, Ttile);
    free(Ttile);
  }
free(fname);
  writetif("joined.tif", target, M, N, P);

  for(size_t kk = 0; kk<M*N*P; kk++)
  { 
    assert( fabs(source[kk]-target[kk]) < 1e-4); 
  }

  free(target);
  free(source);
  tiling_free(T);
  free(T);
  }

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
  } else {
    printf("Please use:\n$ %s M N P maxSize overlap\n", argv[0]);
    exit(1);
  }

  test_getWeight1d();
  test_weights(M, N, P, maxSize, overlap);
  test_copy_paste(M, N, P, maxSize, overlap);

}
