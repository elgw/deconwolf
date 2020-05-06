#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "tiling.h"


int * tiling_getDivision(int M, int m, int * nDiv)
{

  // How many sections in M?
  float ns = (float) ceil(M/m);
  // Size of the sections?
  float width = M/ns;
  printf("nb: %.0f, w: %.0f\n", ns, width);
  int * divs = malloc(2*ns*sizeof(int));
  divs[0] = 0;
  divs[(int)( 2*ns-1)] = M-1;

  // Settle the end and start of all regions between the end points
  for(int dd = 1 ; dd < ns; dd++)
  {
    divs[dd*2-1] = (int) ceil(width*dd);
    divs[dd*2] = divs[dd*2-1]+1;  
  }
  nDiv[0] = (int) ns;
  return divs;
}

int imin(int a, int b)
{
  if(a < b){ return(a); } else { return(b); } ;
}

int imax(int a, int b)
{
  if(a > b){ return(a); } else { return(b); } ;
}

tiling * tiling_create(int M, int N, int P, int maxSize, int overlap)
{

  int nM = 0;
  int * divM = tiling_getDivision(M, maxSize, &nM);
  int nN = 0;
  int * divN = tiling_getDivision(N, maxSize, &nN);

  printf("Dividing %d into:\n", M);
  for(int kk = 0; kk<nM; kk++)
  {
    printf("[%d, %d]\n", divM[2*kk], divM[2*kk+1]);
  }

  // Create tiles
  tiling * T = malloc(sizeof(tiling));
  T->maxSize = maxSize;
  T->overlap = overlap;
  T->nTiles = nM*nN;
  T->tiles = malloc(T->nTiles * sizeof(tile*));
  T->M = M;
  T->N = N;
  T->P = P;
  T->maxSize = maxSize;
  T->overlap = overlap;

  int bb = 0;
  for(int mm = 0; mm<nM; mm++)
  {
    for(int nn = 0; nn<nN; nn++)
    {
      T->tiles[bb] = tile_create();
      tile * t = T->tiles[bb]; 
      t->size[0] = divM[mm*2+1] - divM[mm*2];
      t->size[1] = divN[nn*2+1] - divN[nn*2];
      t->size[2] = P;
      t->pos[0]=divM[mm*2]; t->pos[1]=divM[mm*2+1];
      t->pos[2]=divN[nn*2]; t->pos[3]=divN[nn*2+1];
      t->pos[4]=0; t->pos[5] = P-1;

      t->xpos[0] = imax(0, t->pos[0]-overlap);
      t->xpos[1] = imin(t->pos[1]+overlap, M-1);
      t->xpos[2] = imax(0, t->pos[2]-overlap);
      t->xpos[3] = imin(t->pos[3]+overlap, N-1);
      t->xpos[4] = t->pos[4]; t->xpos[5] = t->pos[5];

      bb++;
    }
  }
  free(divM);
  free(divN);
  return T;
}

void tiling_show(tiling * T)
{
  printf("Tiling with %d tiles\n", T->nTiles);
  printf("Generated for [%d x %d x %d], maxSize: %d, overlap: %d\n",
      T->M, T->N, T->P, T->maxSize, T->overlap);
  for(int kk = 0; kk<T->nTiles; kk++)
  {
    tile_show(T->tiles[kk]);
  }
}

void tiling_free(tiling * T)
{
  for(int kk = 0; kk<T->nTiles; kk++)
  {
    tile_free(T->tiles[kk]);
    free(T->tiles[kk]);
  }
  free(T->tiles);
}

void tile_show(tile * T)
{
  printf("tile: [%d x %d x %d]\n", 
      T->size[0], T->size[1], T->size[2]);
  printf("pos: %d--%d, %d--%d, %d--%d\n", 
      T->pos[0], T->pos[1], T->pos[2],
      T->pos[3], T->pos[4], T->pos[5]);
  printf("xpos: %d--%d, %d--%d, %d--%d\n", 
      T->xpos[0], T->xpos[1], T->xpos[2],
      T->xpos[3], T->xpos[4], T->xpos[5]);
}


float getWeight1d(float a, float b, float c, float d, int x)
{
  // f(x) = 0, x<a, or x>d
  //        1    b < x < c
  //        and linear ramp between a,b and c,d
  if(x<=a || x>=d)
  {
    return 0;
  }
  if(x>= b && x<=c)
  {
    return 1;  
  }
  if(x>=a && x <=b)
  {
    return (x-a)/(b-a);
  }
if(x>=c && x <=d)
  {
    return 1 -(x-c)/(d-c);
  }
assert(0);
}

float tile_getWeight(tile * t, int m, int n, int p, float pad)
{
  assert(p>=t->xpos[4]);
  assert(p<=t->xpos[5]);
  float wm = getWeight1d(t->xpos[0], t->pos[0], t->pos[1], t->xpos[1], m);
  float wn = getWeight1d(t->xpos[2], t->pos[2], t->pos[3], t->xpos[3], n);
//  float w = sqrt( pow((wm-pad), 2) + pow((wn-pad),2));
  float w = fmin(wm,wn);
//  printf("wm: %f, wn: %f, w: %f\n", wm, wn, w);
  assert(w>= 0);
  assert(w<=1);
  return w;
}

float tiling_getWeights(tiling * T, int M, int N, int P)
{
  float w = 0;
  const float pad = (float) T->overlap;
  for(int tt = 0; tt<T->nTiles; tt++)
  {
    w+=tile_getWeight(T->tiles[tt], M, N, P, pad);
  }
  return w;
}
tile * tile_create()
{
  tile * t = malloc(sizeof(tile));
  t->size = malloc(3*sizeof(int));
  t->pos = malloc(6*sizeof(int));
  t->xpos = malloc(6*sizeof(int));
  return t;
}

void tile_free(tile * t)
{
  free(t->size);
  free(t->pos);
  free(t->xpos);
}
