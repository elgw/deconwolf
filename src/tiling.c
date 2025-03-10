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

#include "tiling.h"
#include "fim_tiff.h"
#include "dw_util.h"
#include "fim.h"

int64_t * tiling_getDivision(const int64_t M, const int64_t m, int64_t * nDiv)
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
    assert(ns > 0);
    int64_t * divs = calloc(2*ns, sizeof(int64_t));
    assert(divs != NULL);
    // Start point for the first tile
    divs[0] = 0;
    // End point for the last tile
    divs[(int)( 2*ns-1)] = M-1;

    // Settle the end and start of all regions between the end points
    for(int64_t dd = 1 ; dd < ns; dd++)
    {
        divs[dd*2-1] = (int) ceil(width*dd)-1;
        divs[dd*2] = divs[dd*2-1]+1;
    }
    nDiv[0] = (int) ns;

#ifndef NDEBUG
    for(int64_t dd = 1; dd < ns; dd++)
    {
        assert( (divs[dd*2 - 1] + 1) == divs[dd*2] );
    }
#endif

    return divs;
}

int64_t imin(int64_t a, int64_t b)
{
    if(a < b){ return(a); } else { return(b); } ;
}

int64_t imax(int64_t a, int64_t b)
{
    if(a > b){ return(a); } else { return(b); } ;
}

tiling * tiling_create(const int64_t M, const int64_t N, const int64_t P, const int64_t maxSize, const int64_t overlap)
{

    int64_t nM = 0;
    int64_t * divM = tiling_getDivision(M, maxSize, &nM);
    int64_t nN = 0;
    int64_t * divN = tiling_getDivision(N, maxSize, &nN);

#ifndef NDEBUG
    printf("Dividing %" PRId64 " into:\n", M);
    for(int64_t kk = 0; kk<nM; kk++)
    {
        printf("[%" PRId64 ", %" PRId64 "]\n", divM[2*kk], divM[2*kk+1]);
    }
#endif

    // Create tiles
    tiling * T = malloc(sizeof(tiling));
    T->maxSize = maxSize;
    T->overlap = overlap;
    T->nTiles = nM*nN;
    T->tiles = calloc(T->nTiles, sizeof(tile*));
    assert(T->tiles != NULL);
    T->M = M;
    T->N = N;
    T->P = P;
    T->maxSize = maxSize;
    T->overlap = overlap;

    int64_t bb = 0;
    for(int64_t mm = 0; mm<nM; mm++)
    {
        for(int64_t nn = 0; nn<nN; nn++)
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
    printf("Generated for [%" PRId64 " x %" PRId64 " x %" PRId64 "], maxSize: %d, overlap: %d\n",
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
    printf("size: [%" PRId64 " x %" PRId64 " x %" PRId64 "]\n",
           T->size[0], T->size[1], T->size[2]);
    printf("xsize: [%" PRId64 " x %" PRId64 " x %" PRId64 "]\n",
           T->xsize[0], T->xsize[1], T->xsize[2]);
    printf("pos: %" PRId64 "--%" PRId64 ", %" PRId64 "--%" PRId64 ", %" PRId64 "--%" PRId64 "\n",
           T->pos[0], T->pos[1], T->pos[2],
           T->pos[3], T->pos[4], T->pos[5]);
    printf("xpos: %" PRId64 "--%" PRId64 ", %" PRId64 "--%" PRId64 ", %" PRId64 "--%" PRId64 "\n",
           T->xpos[0], T->xpos[1], T->xpos[2],
           T->xpos[3], T->xpos[4], T->xpos[5]);
}


float getWeight1d(const float a, const float b, const float c, const float d, const int64_t x)
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
                     const int64_t m, const int64_t n, const int64_t p)
/* Calculate the weight for tile t at position (m,n,p) */
{
    ((void) p);
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

float tiling_getWeights(tiling * T, const int64_t M, const int64_t N, const int64_t P)
{
    float w = 0;
    for(int tt = 0; tt < T->nTiles; tt++)
    {
        w+=tile_getWeight(T->tiles[tt], M, N, P);
    }
    return w;
}

tile * tile_create()
{
    tile * t = malloc(sizeof(tile));
    t->size = malloc(3*sizeof(int64_t));
    t->xsize = malloc(3*sizeof(int64_t));
    t->pos = malloc(6*sizeof(int64_t));
    t->xpos = malloc(6*sizeof(int64_t));
    return t;
}

void tile_free(tile * t)
{
    free(t->size);
    free(t->xsize);
    free(t->pos);
    free(t->xpos);
}

float * tiling_get_tile_raw(tiling * T, const int tid, const char * fName)
{
//  printf("Reading tile from %s\n", fName);
    tile * t = T->tiles[tid];
    FILE * fid = fopen(fName, "rb");
    if(fid == NULL)
    {
        printf("ERROR: Can't read %s\n", fName);
        exit(1);
    }

    size_t m = t->xsize[0];
    size_t n = t->xsize[1];
    size_t p = t->xsize[2];

    size_t npixels = m*n*p;
//  printf("To populate %zu pixels (%zu x %zu x %zu)\n", npixels, m, n, p);
    float * R = fim_malloc(npixels*sizeof(float));
    if(R == NULL)
    {
        printf("ERROR: memory allocation failed\n");
        exit(-1);
    }
    memset(R, 0, npixels*sizeof(float));

    for(size_t pp = t->xpos[4]; pp <= (size_t) t->xpos[5]; pp++)
    {
        for(size_t nn = t->xpos[2]; nn <= (size_t) t->xpos[3]; nn++)
        {
            // seek position in big raw file
            size_t spos = t->xpos[0] + nn*T->M + pp*T->M*T->N;
            // write position in tile
            size_t wpos = (nn - t->xpos[2])*t->xsize[0] +
                (pp - t->xpos[4])*t->xsize[0]*t->xsize[1]; // in tile
//      printf("spos: %zu, wpos: %zu\n", spos, wpos);
            assert(wpos+m <= npixels);
            assert(wpos <= spos);
            dw_fseek(fid, spos*sizeof(float), SEEK_SET);
            errno = 0;
            size_t nread = fread(R+wpos, sizeof(float), m, fid);
            if(nread != m)
            {
                perror("fread error");
                fprintf(stderr,
                        "Error reading from %s. Read %zu elements, expected %zu\n"
                        "In %s %d\n",
                        fName, nread, m,
                        __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
        }
    }

    fclose(fid);
    return R;
}

float * tiling_get_tile_tiff(tiling * T, const int tid, const char * fName)
{
    tile * t = T->tiles[tid];
    int verbosity = 1;
    int64_t M = 0; int64_t N = 0; int64_t P = 0; // Will be set to the image size
    float * R = fim_tiff_read_sub(fName, NULL, &M, &N, &P, verbosity,
                                  1,
                                  t->xpos[0], t->xpos[2], t->xpos[4], // Start pos
                                  t->xsize[0], t->xsize[1], t->xsize[2]); // size

    printf("%" PRId64 " %" PRId64 " %" PRId64 "\n", t->xsize[0], t->xsize[1], t->xsize[2]);

    if(0)
    {
        printf("Writing to tile.tif\n");
        fim_tiff_write("tile.tif", R, NULL, t->xsize[0], t->xsize[1], t->xsize[2]);
        printf("ok\n"); getchar();
    }
    return R;
}

float * tiling_get_tile(tiling * T, const int tid, const float * restrict V)
/* Extract tile number tid from the image V */
{
    tile * t = T->tiles[tid];
    int64_t M = T->M; int64_t N = T->N;
#ifndef NDEBUG
    int64_t P = T->P;
    tile_show(t);
#endif
    int64_t m = t->xsize[0];
    int64_t n = t->xsize[1];
    int64_t p = t->xsize[2];
    float * R = fim_malloc(m*n*p*sizeof(float));
    for(int64_t cc = t->xpos[4]; cc <= t->xpos[5]; cc++)
    {
        for(int64_t bb = t->xpos[2]; bb <= t->xpos[3]; bb++)
        {
            for(int64_t aa = t->xpos[0]; aa <= t->xpos[1]; aa++)
            {
                // printf("aa:%d, bb:%d, cc:%d\n", aa, bb, cc);
                size_t Vidx = aa + bb*M + cc*M*N;
                assert(Vidx < (size_t) M*N*P);
                // New coordinates are offset ...
                size_t Ridx = (aa - t->xpos[0]) +
                    (bb - t->xpos[2])*m +
                    (cc - t->xpos[4])*m*n;
                assert(Ridx < (size_t) m*n*p);
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
    int64_t M = T->M; int64_t N = T->N;
    int64_t m = t->xsize[0];
    int64_t n = t->xsize[1];

//  printf("Opening %s for r/w\n", fname); fflush(stdout);
    FILE * fid = fopen(fname, "rb+");
    if(fid == NULL)
    {
        printf("ERROR: Failed to open %s\n", fname);
        exit(1);
    }

    size_t buf_size = M*sizeof(float);
    float * buf = calloc(buf_size, 1);
    assert(buf != NULL);

    for(int64_t cc = t->xpos[4]; cc <= t->xpos[5]; cc++)
    {
        for(int64_t bb = t->xpos[2]; bb <= t->xpos[3]; bb++)
        {
            size_t colpos = (bb*M + cc*M*N)*sizeof(float);
            //      printf("colpos: %zu\n", colpos);
            //fsetpos(fid, &colpos);
            dw_fseek(fid, colpos, SEEK_SET);
            size_t nread = fread(buf, buf_size, 1, fid);
            (void) nread;
            size_t buf_pos = t->xpos[0];
            for(int64_t aa = t->xpos[0]; aa <= t->xpos[1]; aa++)
            {
                // Index in the tile ...
                size_t Sidx = (aa - t->xpos[0]) +
                    (bb - t->xpos[2])*m +
                    (cc - t->xpos[4])*m*n;
                float w = tile_getWeight(t, aa, bb, cc);
                w/= tiling_getWeights(T, aa, bb, cc);
                buf[buf_pos++] += w*(float) S[Sidx];
            }

            dw_fseek(fid, colpos, SEEK_SET);
            //fsetpos(fid, &colpos);
            fwrite(buf, buf_size, 1, fid);
        }
    }
    fclose(fid);
    free(buf);

}

void tiling_put_tile(tiling * T, int tid, float * restrict V, float * restrict S)
{
    /* Using the Tiling T, and the tile tid, write the contents of the tile S
       into the target volume V
       Typically V should to be initialized to 0 first.
    */

    tile * t = T->tiles[tid];
    int64_t M = T->M; int64_t N = T->N;
    int64_t m = t->xsize[0];
    int64_t n = t->xsize[1];
    for(int64_t cc = t->xpos[4]; cc <= t->xpos[5]; cc++)
    {
        for(int64_t bb = t->xpos[2]; bb <= t->xpos[3]; bb++)
        {
            for(int64_t aa = t->xpos[0]; aa <= t->xpos[1]; aa++)
            {
                size_t Vidx = aa + bb*M + cc*M*N;
                size_t Sidx = (aa-t->xpos[0]) +
                    (bb-t->xpos[2])*m +
                    (cc-t->xpos[4])*m*n;
                float w = tile_getWeight(t, aa, bb, cc);
                w/= tiling_getWeights(T, aa, bb, cc);
                V[Vidx] += w*S[Sidx];
            }
        }
    }

}
