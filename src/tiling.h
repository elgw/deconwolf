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
float * tiling_get_tile(tiling * T, int t, const float * restrict V);

/* Extract tile #t from tiff file */
float * tiling_get_tile_tiff(tiling * T, int t, const char * fName);

/* Extract tile #t from raw float file */
float * tiling_get_tile_raw(tiling * T, int t, const char * fName);


/* Put back data extracted by tiling_get_tile
 * S extracted data from tile t
 * V target image, dimensions given by T->M, N, P
 * */
void tiling_put_tile(tiling * T, int t, float * V, float * S);

/* Put back data extracted by tiling_get_tile
 * S extracted data from tile t
 * V target image, dimensions given by T->M, N, P
 * */
void tiling_put_tile_raw(tiling * T, int t, const char * fName, float * S);

tile * tile_create();
void tile_free(tile *);
void tile_show(tile *);
float tile_getWeight(tile *, int m, int n, int p, float pad);
float getWeight1d(float a, float b, float c, float d, int x);

#endif
