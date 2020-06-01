/*    Copyright (C) 2020 Erik L. G. Wernersson
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <fftw3.h>
#include "tiling.h"
#include "fim_tiff.h"

int * tiling_getDivision(const int M, const int m, int * nDiv)
{

  // How many sections in M?
  float ns = (float) ceil((double) M/ (double) m);
  //  printf("%.0f sections\n", ns);
  //  exit(1);
  // Size of the sections?
  float width = M/ns;
#ifndef NDEBUG
  printf("nb: %.0f, w: %.0f\n", ns, width);
#endif
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

tiling * tiling_create(const int M, const int N, const int P, const int maxSize, const int overlap)
{

  int nM = 0;
  int * divM = tiling_getDivision(M, maxSize, &nM);
  int nN = 0;
  int * divN = tiling_getDivision(N, maxSize, &nN);

#ifndef NDEBUG
  printf("Dividing %d into:\n", M);
  for(int kk = 0; kk<nM; kk++)
  {
    printf("[%d, %d]\n", divM[2*kk], divM[2*kk+1]);
  }
#endif

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
      t->size[0] = divM[mm*2+1] - divM[mm*2] + 1;
      t->size[1] = divN[nn*2+1] - divN[nn*2] + 1;
      t->size[2] = P;

      t->pos[0]=divM[mm*2]; t->pos[1]=divM[mm*2+1];
      t->pos[2]=divN[nn*2]; t->pos[3]=divN[nn*2+1];
      t->pos[4]=0; t->pos[5] = P-1;

      t->xpos[0] = imax(0, t->pos[0]-overlap);
      t->xpos[1] = imin(t->pos[1]+overlap, M-1);
      t->xpos[2] = imax(0, t->pos[2]-overlap);
      t->xpos[3] = imin(t->pos[3]+overlap, N-1);
      t->xpos[4] = t->pos[4]; t->xpos[5] = t->pos[5];

      t->xsize[0] = t->xpos[1] - t->xpos[0] + 1;
      t->xsize[1] = t->xpos[3] - t->xpos[2] + 1;
      t->xsize[2] = t->xpos[5] - t->xpos[4] + 1;

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
  printf("-> tile_show\n");
  printf("size: [%d x %d x %d]\n", 
      T->size[0], T->size[1], T->size[2]);
  printf("xsize: [%d x %d x %d]\n", 
      T->xsize[0], T->xsize[1], T->xsize[2]);
  printf("pos: %d--%d, %d--%d, %d--%d\n", 
      T->pos[0], T->pos[1], T->pos[2],
      T->pos[3], T->pos[4], T->pos[5]);
  printf("xpos: %d--%d, %d--%d, %d--%d\n", 
      T->xpos[0], T->xpos[1], T->xpos[2],
      T->xpos[3], T->xpos[4], T->xpos[5]);
}


float getWeight1d(const float a, const float b, const float c, const float d, const int x)
{
  assert(a<=b);
  assert(b<c); // a tile has to have a size
  assert(c<=d);
  // f(x) = 0, x<a, or x>d
  //        1    b < x < c
  //        and linear ramp between a,b and c,d
  if(x<a || x>d)
  {
    return 0;
  }
  if(x>= b && x<=c)
  {
    return 1;  
  }
  if(x>=a && x <=b)
  {
    return (x-a+1)/(b-a+1);
  }
  if(x>=c && x <=d)
  {
    return 1 -(x-c)/(d-c+1);
  }
  assert(0);
  return -1;
}

float tile_getWeight(tile * t, 
    const int m, const int n, const int p, 
    const float pad)
  /* Calculate the weight for tile t at position (m,n,p) */
{
  assert(p>=t->xpos[4]);
  assert(p<=t->xpos[5]);
  float wm = getWeight1d(t->xpos[0], t->pos[0], t->pos[1], t->xpos[1], m);
  float wn = getWeight1d(t->xpos[2], t->pos[2], t->pos[3], t->xpos[3], n);
  //  float w = sqrt( pow((wm-pad), 2) + pow((wn-pad),2));
  assert(wm>=0); assert(wn>=0); assert(wm<=1); assert(wn<=1);
  float w = fminf(wm,wn);
  //  printf("wm: %f, wn: %f, w: %f\n", wm, wn, w);
  assert(w>= 0);
  assert(w<=1);
  return w;
}

float tiling_getWeights(tiling * T, const int M, const int N, const int P)
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
  t->xsize = malloc(3*sizeof(int));
  t->pos = malloc(6*sizeof(int));
  t->xpos = malloc(6*sizeof(int));
  return t;
}

void tile_free(tile * t)
{
  free(t->size);
  free(t->xsize);
  free(t->pos);
  free(t->xpos);
}

float * tiling_get_tile_tiff(tiling * T, const int tid, const char * fName)
{
  tile * t = T->tiles[tid];
  int verbosity = 1;
  int M = 0; int N = 0; int P = 0; // Will be set to the image size
  float * R = fim_tiff_read_sub(fName, &M, &N, &P, verbosity, 
      1, 
      t->xpos[0], t->xpos[2], t->xpos[4], // Start pos
      t->xsize[0], t->xsize[1], t->xsize[2]); // size

  printf("%d %d %d\n", t->xsize[0], t->xsize[1], t->xsize[2]);

  if(0)
  {
    printf("Writing to tile.tif\n");
    fim_tiff_write("tile.tif", R, t->xsize[0], t->xsize[1], t->xsize[2]);
    printf("ok\n"); getchar();
  }
  return R;
}

float * tiling_get_tile(tiling * T, const int tid, const float * restrict V)
  /* Extract tile number tid from the image V */
{
  tile * t = T->tiles[tid];
  int M = T->M; int N = T->N; 
#ifndef NDEBUG
  int P = T->P;
  tile_show(t);
#endif
  int m = t->xsize[0];
  int n = t->xsize[1];
  int p = t->xsize[2];
  float * R = fftwf_malloc(m*n*p*sizeof(float));
  for(int cc = t->xpos[4]; cc <= t->xpos[5]; cc++)
  {
    for(int bb = t->xpos[2]; bb <= t->xpos[3]; bb++)
    {
      for(int aa = t->xpos[0]; aa <= t->xpos[1]; aa++)
      {
        // printf("aa:%d, bb:%d, cc:%d\n", aa, bb, cc);
        size_t Vidx = aa + bb*M + cc*M*N;
        assert(Vidx < M*N*P);
        // New coordinates are offset ...
        size_t Ridx = (aa - t->xpos[0]) + 
          (bb - t->xpos[2])*m + 
          (cc - t->xpos[4])*m*n;
        assert(Ridx < m*n*p);
        R[Ridx] = V[Vidx];
      }
    }
  }
  return R;
}

// Write tile directly to raw float file
void tiling_put_tile_raw(tiling * T, int tid, const char * fname, float * restrict S)  
{
  /* 
   * Assumes that the raw file is already created and big enough 
   * To reduce the number of fseeks and writes, one column at a time is read/written
   *
   * TIFF files does not support altering of the contents,
   * that is why I settled for this solution.
   * */

  tile * t = T->tiles[tid];
  int M = T->M; int N = T->N; 
  int m = t->xsize[0];
  int n = t->xsize[1];

  printf("Opening %s for r/w\n", fname); fflush(stdout);
  FILE * fid = fopen(fname, "r+");
  if(fid == NULL)
  {
    printf("ERROR: Failed to open %s\n", fname);
    exit(1);
  }

  size_t buf_size = M*sizeof(float);
  float * buf = malloc(buf_size);

  for(int cc = t->xpos[4]; cc <= t->xpos[5]; cc++)
  {
    for(int bb = t->xpos[2]; bb <= t->xpos[3]; bb++)
    {
      size_t colpos = (bb*M + cc*M*N)*sizeof(float);
      //      printf("colpos: %zu\n", colpos);
      //fsetpos(fid, &colpos);
      fseek(fid, colpos, SEEK_SET);
      fread(buf, buf_size, 1, fid);
      size_t buf_pos = t->xpos[0];
      for(int aa = t->xpos[0]; aa <= t->xpos[1]; aa++)
      {
        // Index in the tile ...
        size_t Sidx = (aa - t->xpos[0]) + 
          (bb - t->xpos[2])*m + 
          (cc - t->xpos[4])*m*n;
        float w = tile_getWeight(t, aa, bb, cc, T->overlap);
        w/= tiling_getWeights(T, aa, bb, cc);
        buf[buf_pos++] += w*(float) S[Sidx];
      }

      fseek(fid, colpos, SEEK_SET);
      //fsetpos(fid, &colpos);
      fwrite(buf, buf_size, 1, fid);
    }
  }
  fclose(fid);
  free(buf);

}

void tiling_put_tile(tiling * T, int tid, float * restrict V, float * restrict S)
{
  tile * t = T->tiles[tid];
  int M = T->M; int N = T->N; 
  int m = t->xsize[0];
  int n = t->xsize[1];
  for(int cc = t->xpos[4]; cc <= t->xpos[5]; cc++)
  {
    for(int bb = t->xpos[2]; bb <= t->xpos[3]; bb++)
    {
      for(int aa = t->xpos[0]; aa <= t->xpos[1]; aa++)
      {
        size_t Vidx = aa + bb*M + cc*M*N;
        size_t Sidx = (aa-t->xpos[0]) + 
          (bb-t->xpos[2])*m + 
          (cc-t->xpos[4])*m*n;
        float w = tile_getWeight(t, aa, bb, cc, T->overlap);
        w/= tiling_getWeights(T, aa, bb, cc);
        V[Vidx] += w*S[Sidx];
      }
    }
  }

}
