#ifndef tiling_h
#define tiling_h

typedef struct{
  int * size; // M, N, P
  int * xsize; // M, N, P with padding
  int * pos; // Position in original image (M0, M1, N0, N1, P0, P1)
  int * xpos; // Extended position, including overlap
} tile;

typedef struct{
  int M, N, P;
  int nTiles;
  tile ** tiles;
  int maxSize;
  int overlap;
} tiling;

tiling * tiling_create(int M, int N, int P, int maxSize, int overlap);
void tiling_show(tiling * T);
void tiling_free(tiling * T);
float tiling_getWeights(tiling * T, int m, int n, int p);
/* Extract tile #t from V */
float * tiling_get_tile(tiling * T, int t, float * V);
/* Put back data extracted by tiling_get_tile
 * S extracted data from tile t
 * V target image, dimensions given by T->M, N, P
 * */
void tiling_put_tile(tiling * T, int t, float * V, float * S);


tile * tile_create();
void tile_free(tile *);
void tile_show(tile *);
float tile_getWeight(tile *, int m, int n, int p, float pad);

#endif
