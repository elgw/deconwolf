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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <tiffio.h>
#include <unistd.h>
#include "fim_tiff.h"
#include "fim.h"

/* see man 3 tifflib
 *
 * From: https://www.cs.rochester.edu/u/nelson/courses/vision/resources/tiff/libtiff.html#Errors
 * Finally, note that the last strip of data in an image may have fewer rows in it than specified by the 
 * RowsPerStrip tag. A reader should not assume that each decoded strip contains a full set of rows in it. 
 */


void fim_tiff_ut()
{
  printf("-> fim_tiff_ut (write and read back a tif file)\n");
  /* Create and write a 3D image to disk,
   * read it back and check that it is the same thing*/

  char fname[] = "_deconwolf_temporary_XXXXXX";
  int fd = mkstemp(fname);
  if(fd == -1)
  {
    printf("Could not get a good temporary file name, skipping test\n");
    return;
  }
  close(fd);

  //  printf("%s\n", fname);
  //  getchar();
  int M = 1024, N = 2048, P = 2;
  float * im = fim_zeros(M*N*P);

  size_t pos1 = 1111;
  size_t pos2 = 2222;
  size_t pos3 = 3333;

  im[pos1] = 1;
  im[pos2] = 2;
  im[pos3] = 3;

  writetif(fname, im, M, N, P);

  int M2 = 0, N2 = 0, P2 = 0;
  float * im2 = readtif_asFloat(fname, &M2, &N2, &P2, 0);
  if(im2 == NULL)
  {
    printf("Could not read back the image\n");
    free(im);
    remove(fname);
    return;
  }

  if(M != M2 || N != N2 || P != P2)
  {
    printf("Dimensions does not match!\n");
    printf("Wrote: [%d x %d x %d], Read: [%d x %d x %d]\n",
        M, N, P, M2, N2, P2);
    free(im);
    free(im2);
    remove(fname);
    return;
  }

  if( fabs(im[pos1]/im[pos2] - im2[pos1]/im2[pos2])>1e-5)
  {
    printf("Error: images values does not match\n");
  }

  free(im);
  free(im2);
  remove(fname);
}


void floatimage_normalize(float * restrict I, const size_t N)
{
  // Scale image to span the whole 16 bit range.
  float imax = 1e7; float imin = 1e7;
  int ok = 1;
  for(size_t kk=0; kk<N; kk++)
  {
    if(~isfinite(I[kk]))
    {
      I[kk] = 0;
      ok = 0;
    }
    I[kk] > imax ? imax = I[kk] : 0 ;
    I[kk] < imin ? imin = I[kk] : 0 ;
  }

  if(!ok)
  {
    printf("floatimage_normalize got non-normal numbers\n");
  }

  printf("floatimage_normalize imin: %f imax: %f\n", imin, imax);

  if(imax>0)
  {
    for(size_t kk=0; kk<N; kk++)
      I[kk]*=(pow(2,16)-1)/imax;
  }

}

void floatimage_show_stats(float * I, size_t N, size_t M, size_t P)
{
  float isum = 0;
  float imin = INFINITY;
  float imax = -INFINITY;
  for(size_t kk = 0; kk<M*N*P; kk++)
  {
    I[kk] > imax ? imax = I[kk] : 0 ;
    I[kk] < imin ? imin = I[kk] : 0 ;
    isum += I[kk];
  }
  printf("min: %f max: %f mean: %f\n", imin, imax, isum/(M*N*P));
}

void readUint8(TIFF * tfile, float * V, 
    const uint32_t ssize, 
    const uint32_t ndirs,
    const uint32_t nstrips,
    const uint32_t perDirectory
    )
{
  // Number of elements per strip
  size_t nes = ssize/sizeof(uint8_t);
  uint8_t * buf = _TIFFmalloc(ssize);
  assert(buf != NULL);

  for(int dd=0; dd<ndirs; dd++) {
    TIFFSetDirectory(tfile, dd);
    for(int kk=0; kk<nstrips; kk++) {
      int strip = kk;
      tsize_t read = TIFFReadEncodedStrip(tfile, strip, buf, (tsize_t) - 1);
      assert(read>0);

      for(int ii = 0; ii<read/sizeof(uint8_t); ii++) {
        size_t pos = ii+kk*nes + dd*perDirectory; 
        V[pos] = (float) buf[ii];
      }
    } 
  }
  _TIFFfree(buf);
}

void readUint16(TIFF * tfile, float * V, 
    const uint32_t ssize, 
    const uint32_t ndirs,
    const uint32_t nstrips,
    const uint32_t perDirectory
    )
{
  // Number of elements per strip
  size_t nes = ssize/sizeof(uint16_t);
  uint16_t * buf = _TIFFmalloc(ssize);
  assert(buf != NULL);

  for(int dd=0; dd<ndirs; dd++) {
    TIFFSetDirectory(tfile, dd);
    for(int kk=0; kk<nstrips; kk++) {
      int strip = kk;
      tsize_t read = TIFFReadEncodedStrip(tfile, strip, buf, (tsize_t) - 1);
      assert(read>0);

      for(int ii = 0; ii<read/sizeof(uint16_t); ii++) {
        size_t pos = ii+kk*nes + dd*perDirectory; 
        V[pos] = (float) buf[ii];
      }
    } 
  }
  _TIFFfree(buf);
}

void readFloat(TIFF * tfile, float * V,
    uint32_t ssize, 
    uint32_t ndirs,
    uint32_t nstrips,
    uint32_t perDirectory
    )
{
  // Number of elements per strip
  size_t nes = ssize/sizeof(float);
  float * buf = _TIFFmalloc(ssize);
  assert(buf != NULL);

  for(int dd=0; dd<ndirs; dd++) {
    TIFFSetDirectory(tfile, dd);
    for(int kk=0; kk<nstrips; kk++) {
      int strip = kk;
      tsize_t read = TIFFReadEncodedStrip(tfile, strip, buf, (tsize_t)-1);
      assert(read > 0);

      for(int ii = 0; ii<read/sizeof(float); ii++) {
        size_t pos = ii+kk*nes + dd*perDirectory; 
        V[pos] = buf[ii];
      }
    } 
  }
  _TIFFfree(buf);
}


int writetif(char * fName, float * V, 
    int N, int M, int P)
{
  float imax = -INFINITY;
  for(size_t kk = 0; kk<M*N*P; kk++)
  {
    if(isfinite(V[kk]))
    {
      if(V[kk] > imax)
      {
        imax = V[kk];
      }
    }
  }
  float scaling = 1.0/imax*(pow(2,16)-1.0);
//  printf("scaling: %f\n", scaling);

  size_t bytesPerSample = sizeof(uint16_t);
  TIFF* out = TIFFOpen(fName, "w");
  assert(out != NULL);

  size_t linbytes = (M+N)*bytesPerSample;
  uint16_t * buf = _TIFFmalloc(linbytes);

  for(size_t dd = 0; dd<P; dd++)
  {

    TIFFSetField(out, TIFFTAG_IMAGEWIDTH, N);  // set the width of the image
    TIFFSetField(out, TIFFTAG_IMAGELENGTH, M);    // set the height of the image
    TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1);   // set number of channels per pixel
    TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8*bytesPerSample);    // set the size of the channels
    TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image.
    //   Some other essential fields to set that you do not have to understand for now.
    TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField(out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);

    /* We are writing single page of the multipage file */
    TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
    /* Set the page number */
    TIFFSetField(out, TIFFTAG_PAGENUMBER, dd, P); 


    for(size_t kk = 0; kk<M; kk++)
    {
      for(size_t ll = 0; ll<N; ll++)
      {
        float value = V[M*N*dd + kk*N + ll]*scaling;
        if(!isfinite(value))
        { value = 0; }

        buf[ll] = (uint16_t) round(value);

      }

      
          int ok = TIFFWriteScanline(out, // TIFF
            buf, 
            kk, // row
            0); //sample
        if(ok != 1)
        {
          printf("TIFFWriteScanline failed\n");
        }
    }

    TIFFWriteDirectory(out);
  }

  _TIFFfree(buf);

  TIFFClose(out);
  return 0;
}

float * readtif_asFloat(char * fName, 
    int * N0, int * M0, int * P0, int verbosity)
{
  /* Reads the content of the tif file with fName
   * Puts the images size in M0, N0, P0
   * Fortran or C order?
   * Should be able to read also floating point images etc...
   */

  TIFF * tfile = TIFFOpen(fName, "r");

  if(tfile == NULL) {
    return NULL;
  }

  // Tags: ImageWidth and ImageLength
  uint32_t M = 0, N = 0, BPS = 0, SF = 0, CMP=0;
  TIFFGetField(tfile, TIFFTAG_IMAGELENGTH, &M);
  TIFFGetField(tfile, TIFFTAG_IMAGEWIDTH, &N);
  TIFFGetField(tfile, TIFFTAG_BITSPERSAMPLE, &BPS);
  int gotSF = TIFFGetField(tfile, TIFFTAG_SAMPLEFORMAT, &SF);
  int gotCMP = TIFFGetField(tfile, TIFFTAG_COMPRESSION, &CMP);

  if(gotCMP)
  {
    if(CMP != 1)
    {
      printf("TIFFTAG_COMPRESSION=%u is not supported\n", CMP);
      return NULL;
    }
  }

  int isUint = 0;
  int isFloat = 0;

  if(gotSF)
  {
    if( SF == SAMPLEFORMAT_UINT)
    {
      isUint = 1;
    }
    if( SF == SAMPLEFORMAT_IEEEFP)
    {
      isFloat = 1;
    }
  } 
  else {
    printf("Warning: TIFFTAG_SAMPLEFORMAT not specified, assuming uint but that could be wrong!\n");
    isUint = 1;
  }

  if(!(isUint || isFloat))
  {
    printf("Only unsigned integer images are supported\n");
    return NULL;
  }

  if(isUint && !(BPS == 16 || BPS == 8))
  {
    printf("Unsigned %d-bit images are not supported, only 8 and 16.\n", BPS);
    return NULL;
  }

  if(isFloat && (BPS != 32))
  {
    printf("For floating point images, only 32-bit samples are supported %u\n", BPS);
  }

  M0[0] = (size_t) M;
  N0[0] = (size_t) N;

  uint32_t B = 0;
  int gotB = TIFFGetField(tfile, TIFFTAG_IMAGEDEPTH, &B);

  tmsize_t ssize = TIFFStripSize(tfile); // Seems to be in bytes
  uint32_t nstrips = TIFFNumberOfStrips(tfile);
  uint32_t ndirs = TIFFNumberOfDirectories(tfile);
  P0[0] = (size_t) ndirs;

  if(verbosity > 1)
  {
    if(gotB){
      printf(" TIFFTAG_IMAGEDEPTH: %u\n", B);}
    printf(" TIFFTAG_BITSPERSAMPLE: %u\n", BPS);
    printf(" size: %zu x %zu, %zu bits\n", (size_t) M, (size_t) N, (size_t) BPS);
    printf(" # strips: %zu \n", (size_t) nstrips);
    printf(" strip size (ssize): %zu bytes \n", (size_t) ssize);
    printf(" # dirs (slices): %zu\n", (size_t) ndirs);
  }

  if(TIFFIsTiled(tfile))
  {
    printf("Tiled files are not supported\n");
    TIFFClose(tfile);
    return NULL;
  }

  // assert(M*N*BPS/8 == nstrips * ssize);

  //  size_t nel = nstrips * ssize * ndirs; 
  size_t nel = M*N*ndirs;
  float * V = malloc(nel*sizeof(float));

  if(isFloat)
  {
    if(verbosity > 1)
    {
      printf("ReadFloat ...\n");
    }
    readFloat(tfile, V, ssize, ndirs, nstrips, M*N);
  }
  if(isUint)
  {
    if(verbosity > 1)
    {
      printf("ReadUint ...\n");
    }
    if(BPS == 16)
    {
      readUint16(tfile, V, ssize, ndirs, nstrips, M*N);
    }
    if(BPS == 8)
    {
      readUint8(tfile, V, ssize, ndirs, nstrips, M*N);
    }
  }

  TIFFClose(tfile);

  return V;
}

#ifdef unittest
int main(int argc, char ** argv) 
{

  char * outname = NULL;
  char * inname = NULL;

  if(argc == 1)
  {
    printf("Use %s file.tif\n", argv[0]);
    exit(1);
  }
  inname = argv[1];

  if(argc>2)
  {
    outname = argv[2];
  } else {
    outname = malloc(100*sizeof(char));
    sprintf(outname, "foo.tif");
  }
  printf("Will read from %s and write to %s.\n", inname, outname);

  size_t M = 0, N = 0, P = 0;

  float * I = (float *) readtif_asFloat(inname, &M, &N, &P);

  if(I == NULL)
  {
    printf("Failed to read %s\n", inname);
    exit(1);
  }

  floatimage_show_stats(I, M, N, P);

  floatimage_normalize(I, M*N*P);
  writetif(outname, I, M, N, P);

  free(I);

  return 0;
}
#endif
