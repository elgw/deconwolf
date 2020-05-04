#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <tiffio.h>

/* see man 3 tifflib
 *
 * From: https://www.cs.rochester.edu/u/nelson/courses/vision/resources/tiff/libtiff.html#Errors
 * Finally, note that the last strip of data in an image may have fewer rows in it than specified by the 
 * RowsPerStrip tag. A reader should not assume that each decoded strip contains a full set of rows in it. 
 */

void floatimage_normalize(float * restrict I, const size_t N)
{
  // Scale image to span the whole 16 bit range.
  float imax = 0;
  for(size_t kk=0; kk<N; kk++)
    I[kk] > imax ? imax = I[kk] : 0 ;
  //printf("imax: %f\n", imax);
  for(size_t kk=0; kk<N; kk++)
    I[kk]*=(pow(2,16)-1)/imax;
}

void floatimage_show_stats(float * I, size_t M, size_t N, size_t P)
{
  float isum = 0;
  float imin = 100000;
  float imax = 0;
  for(size_t kk = 0; kk<M*N*P; kk++)
  {
    I[kk] > imax ? imax = I[kk] : 0 ;
    I[kk] < imin ? imin = I[kk] : 0 ;
    isum += I[kk];
  }
  printf("min: %f max: %f mean: %f\n", imin, imax, isum/(M*N*P));
}

void readUint(TIFF * tfile, float * V, 
    uint32_t ssize, 
    uint32_t ndirs,
    uint32_t nstrips,
    uint32_t perDirectory
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
    int M, int N, int P)
{

  size_t bytesPerSample = sizeof(uint16_t);
  TIFF* out = TIFFOpen(fName, "w");
  assert(out != NULL);

  size_t linbytes = M*bytesPerSample;
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


    for(size_t kk = 0; kk<N; kk++)
    {
      for(size_t ll = 0; ll<M; ll++)
      {
        buf[ll] = V[M*N*dd + kk*M + ll];
      }

      assert(TIFFWriteScanline(out, // TIFF
            buf, 
            kk, // row
            0) //sample
          == 1);
    }

    TIFFWriteDirectory(out);
  }

  _TIFFfree(buf);

  TIFFClose(out);
  return 0;
}

float * readtif_asFloat(char * fName, 
    int * M0, int * N0, int * P0)
{
  /* Reads the content of the tif file with fName
   * Puts the images size in M0, N0, P0
   * Fortran or C order?
   * Should be able to read also floating point images etc...
   */

  TIFF * tfile = TIFFOpen(fName, "r");

  if(tfile == NULL) {
    printf("TIFFOpen failed\n");
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
    printf("Warning: TIFFTAG_SAMPLEFORMAT not specified, assuming uint\n");
  }

  if(!(isUint || isFloat))
  {
    printf("Only unsigned integer images are supported\n");
    return NULL;
  }

  if(isUint && (BPS != 16))
  {
    printf("For unsigned images, only 16-bit samples are supported %u\n", BPS);
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

  if(gotB){
    printf(" TIFFTAG_IMAGEDEPTH: %u\n", B);}
  printf(" TIFFTAG_BITSPERSAMPLE: %u\n", BPS);
  printf(" size: %zu x %zu, %zu bits\n", (size_t) M, (size_t) N, (size_t) BPS);
  printf(" # strips: %zu \n", (size_t) nstrips);
  printf(" strip size (ssize): %zu bytes \n", (size_t) ssize);
  printf(" # dirs (slices): %zu\n", (size_t) ndirs);


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
    printf("ReadFloat ...\n");
    readFloat(tfile, V, ssize, ndirs, nstrips, M*N);
  }
  if(isUint)
  {
    printf("ReadUint ...\n");
    readUint(tfile, V, ssize, ndirs, nstrips, M*N);
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
