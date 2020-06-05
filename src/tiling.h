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
  int64_t * size; // M, N, P
  int64_t * xsize; // M, N, P with padding
  int64_t * pos; // Position in original image (M0, M1, N0, N1, P0, P1)
  int64_t * xpos; // Extended position, including overlap
} tile;

typedef struct{
  int64_t M, N, P;
  int nTiles;
  tile ** tiles;
  int maxSize;
  int overlap;
} tiling;

tiling * tiling_create(int64_t M, int64_t N, int64_t P, int64_t maxSize, int64_t overlap);
void tiling_show(tiling * T);
void tiling_free(tiling * T);
float tiling_getWeights(tiling * T, int64_t m, int64_t n, int64_t p);

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
float tile_getWeight(tile *, int64_t m, int64_t n, int64_t p, float pad);
float getWeight1d(float a, float b, float c, float d, int64_t x);

#endif
