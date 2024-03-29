#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tiling.h"
#include "fim.h"
#include "fim_tiff.h"

void test_getWeight1d()
{
  float a = 1, b = 6, c = 8, d = 12;
  printf(" -> test_getWeight1d\n");
  printf("Assuming that the tile goes from %.0f to %.0f\n", b, c);
  printf("See the values over the range from %.0f to %.0f\n", a, d);

  for(int kk = 0; kk<14; kk++)
  {
    float w = getWeight1d(a, b, c, d, kk);
    printf("%d %.2f\n", kk, w);
    if(kk >= b && kk <= c)
    {
      assert(w == 1);
    } else {
      assert(w < 1);
    }
  }
  printf("\n");
}


void test_weights(int M, int N, int P, int maxSize, int overlap)
  {
printf(" -> test_weights\n");
  tiling * T = tiling_create(M,N,P, maxSize, overlap);
  tiling_show(T);

  printf("Generating weight map ... \n"); fflush(stdout);
  float * W = malloc(M*N*sizeof(float));
  assert(W != NULL);
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
  printf("Writing weights from all tiles to %s\n", outFile); fflush(stdout);

  // Assuming that is it same behaviour for all P
  for(size_t kk = 0; kk< (size_t) M*N; kk++)
  {
    assert(W[kk] > 0);
  }


  fim_tiff_write(outFile, W, NULL, M, N, 1);

  for(int aa = 0; aa<M; aa++)
  {
    for(int bb = 0; bb<N; bb++)
    {
      float w = tile_getWeight(T->tiles[1], aa, bb, 0);
      assert(w>=0);
      W[aa + bb*M] = w;
    }
  }

  char outFileTile[] = "tiling_ut_weights_tile1.tif";
  printf("Writing weights from the first tile to %s\n", outFileTile); fflush(stdout);
  fim_tiff_write(outFileTile, W, NULL, M, N, 1);


  // Assuming that is it same behaviour for all P
  for(size_t kk = 0; kk< (size_t) M*N; kk++)
  {
    assert(W[kk] >= 0);
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
  size_t MNP = M*N*P;

  float * source = malloc(MNP*sizeof(float));
  assert(source != NULL);
  for(size_t kk = 0; kk < MNP; kk++)
  {
    source[kk] = 7; //rand()/RAND_MAX;
  }

  float * target = calloc(MNP, sizeof(float));
  assert(target != NULL);
  char * fname = malloc(100);
  assert(fname != NULL);

  for(int tt = 0; tt < T->nTiles; tt++)
  {
    float * Ttile = tiling_get_tile(T, tt, source);
    tile * tl = T->tiles[tt];
    size_t mnp = tl->xsize[0]*tl->xsize[1]*tl->xsize[2];
    for(size_t uu = 0; uu<mnp; uu++)
    {
        assert(Ttile[uu] == 7);
    }

    sprintf(fname, "tile%d.tif", tt);
    fim_tiff_write(fname, Ttile, NULL,
                   tl->xsize[0], tl->xsize[1], tl->xsize[2]);

    // P has size T->tiles[kk]->xsize;
    tiling_put_tile(T, tt, target, Ttile);
    free(Ttile);
  }

  free(fname);
  fim_tiff_write("joined.tif", target, NULL, M, N, P);

  if(MNP < 500)
  {
      for(int kk = 0; kk< M; kk++)
      {
          for(int ll = 0; ll < N; ll++)
          {
              printf("%.0f/%.0f ", source[kk + ll*M], target[kk + ll*M]);
          }
          printf("\n");
      }
  }

  for(size_t kk = 0; kk < MNP; kk++)
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
