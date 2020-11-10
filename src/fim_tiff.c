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

#include "fim_tiff.h"
#include "fim.h"

/* see man 3 tifflib
 *
 * From: https://www.cs.rochester.edu/u/nelson/courses/vision/resources/tiff/libtiff.html#Errors
 * Finally, note that the last strip of data in an image may have fewer rows in it than specified by the
 * RowsPerStrip tag. A reader should not assume that each decoded strip contains a full set of rows in it.
 */

//typedef float afloat __attribute__ ((__aligned__(16)));
typedef float afloat;


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
  int64_t M = 1024, N = 2048, P = 2;
  afloat * im = fim_zeros(M*N*P);

  size_t pos1 = 1111;
  size_t pos2 = 2222;
  size_t pos3 = 3333;

  im[pos1] = 1;
  im[pos2] = 2;
  im[pos3] = 3;

  fim_tiff_write(fname, im, M, N, P, stdout);

  int64_t M2 = 0, N2 = 0, P2 = 0;
  afloat * im2 = fim_tiff_read(fname, &M2, &N2, &P2, 0);
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
    printf("Wrote: [%" PRId64 " x %" PRId64 " x %" PRId64 "], Read: [%" PRId64 " x %" PRId64 " x %" PRId64 "]\n",
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


void floatimage_normalize(afloat * restrict I, const size_t N)
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

void floatimage_show_stats(afloat * I, size_t N, size_t M, size_t P)
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


void readUint8_sub(TIFF * tfile, afloat * V,
    const uint32_t ssize,
    const uint32_t ndirs,
    const uint32_t nstrips,
    const uint32_t perDirectory,
    int64_t sM, int64_t sN, int64_t sP, int64_t wM, int64_t wN, int64_t wP)
{
  perror("Not implemented\n");
  exit(1);
}


void readUint8(TIFF * tfile, afloat * V,
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

  for(int64_t dd=0; dd<ndirs; dd++) {
    TIFFSetDirectory(tfile, dd);
    for(int64_t kk=0; kk<nstrips; kk++) {
      int64_t strip = kk;
      tsize_t read = TIFFReadEncodedStrip(tfile, strip, buf, (tsize_t) - 1);
      assert(read>0);

      for(int64_t ii = 0; ii<read/sizeof(uint8_t); ii++) {
        size_t pos = ii+kk*nes + dd*perDirectory;
        V[pos] = (float) buf[ii];
      }
    }
  }
  _TIFFfree(buf);
}

INLINED extern void sub2ind(size_t ind,
    int64_t M, int64_t N, int64_t P,
    int64_t * m, int64_t * n, int64_t * p)
{
  /* If ind is the linear index from a [M,N,P] image, figure out the coordinates of ind
  */

  ldiv_t d = ldiv(ind, M*N);
  p[0] = d.quot;
  int64_t t = d.rem;
  d = ldiv(t, M);
  n[0] = d.quot;
  m[0] = d.rem;
}

void readUint16_sub(TIFF * tfile, afloat * V,
    const uint32_t ssize,
    const uint32_t ndirs,
    const uint32_t nstrips,
    const uint32_t perDirectory,
    int64_t M, int64_t N, int64_t P,
    int64_t sM, int64_t sN, int64_t sP,
    int64_t wM, int64_t wN, int64_t wP)
{

  // V should hold wM*wN*wP uint16

  // Number of elements per strip
  size_t nes = ssize/sizeof(uint16_t);
  uint16_t * buf = _TIFFmalloc(ssize);
  assert(buf != NULL);

  for(int64_t dd=0; dd<ndirs; dd++) {
    TIFFSetDirectory(tfile, dd);
    for(int64_t kk=0; kk<nstrips; kk++) {
      int64_t strip = kk;
      tsize_t read = TIFFReadEncodedStrip(tfile, strip, buf, (tsize_t) - 1);
      assert(read>0);

      for(int64_t ii = 0; ii<read/sizeof(uint16_t); ii++) {
        size_t ipos = ii+kk*nes + dd*perDirectory; // Input position

        int64_t iM = 0; int64_t iN = 0; int64_t iP = 0;

        sub2ind(ipos, M, N, P, &iM, &iN, &iP);

        if(0){
          if(ipos % 100000 == 0){
            printf("%" PRId64 " %" PRId64 " %" PRId64 ", %zu -> (%" PRId64 ", %" PRId64 ", %" PRId64 ")\n", M, N, P, ipos, iM, iN, iP);
            getchar();
          }}

        int64_t oM = (iM-sM);
        int64_t oN = (iN-sN);
        int64_t oP = (iP-sP);

        if(oM >= 0 && oM < wM){
          if(oN >= 0 && oN < wN){
            if(oP >= 0 && oP < wP){
              size_t opos = oM + oN*wM + oP*wM*wN;
              V[opos] = (float) buf[ii];
            }
          }
        }
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
//  uint16_t * buf = _TIFFmalloc(ssize);
  uint16_t * buf = malloc(ssize);

  if(buf == NULL)
  {
      printf("Failed to allocate %u bytes of memory!", ssize);
      fflush(stdout);
      return;
  }

  for(int64_t dd=0; dd<ndirs; dd++) {
    int ok = TIFFSetDirectory(tfile, dd);
    if(ok == 0)
    {
        printf("Failed to choose directory %lu\n", dd);
        printf("Either the file is corrupt or this program has a bug\n");
        fflush(stdout);
    }
    for(int64_t kk=0; kk<nstrips; kk++) {
      int64_t strip = kk;
      tsize_t read = TIFFReadEncodedStrip(tfile,
                                          (tstrip_t) strip, // Strip number
                                          (tdata_t) buf, // target
                                          (tsize_t) - 1); // The entire strip
      if(read == -1)
      {
          printf("Failed to read from tif file\n");
          fflush(stdout);
      }
      if(read > ssize)
      {
          printf("Read to much!\n"); fflush(stdout);
      }
      printf("\r%lu/%lu", dd, kk); fflush(stdout);
      for(size_t ii = 0; ii < read/sizeof(uint16_t); ii++)
      {
          size_t pos = (size_t)( ii+kk*nes + dd*perDirectory);
        V[pos] = buf[ii];
      }
    }
  }
  printf("...\n"); fflush(stdout);
  _TIFFfree(buf);
}

void readFloat_sub(TIFF * tfile, afloat * V,
    uint32_t ssize,
    uint32_t ndirs,
    uint32_t nstrips,
    uint32_t perDirectory,
    int64_t sM, int64_t sN, int64_t sP, int64_t wM, int64_t wN, int64_t wP)
{
  perror("Not implemented!\n");
  exit(1);
}

void readFloat(TIFF * tfile, afloat * V,
    uint32_t ssize,
    uint32_t ndirs,
    uint32_t nstrips,
    uint32_t perDirectory
    )
{
  // Number of elements per strip
  size_t nes = ssize/sizeof(float);
  float * buf = malloc(ssize);
  assert(buf != NULL);

  for(int64_t dd=0; dd<ndirs; dd++) {
    int ok = TIFFSetDirectory(tfile, dd);
    if(ok == 0)
    {
        printf("Failed TIFFSetDirectory\n");
        fflush(stdout);
    }
    for(int64_t kk=0; kk<nstrips; kk++) {
      int64_t strip = kk;
      tsize_t read = TIFFReadEncodedStrip(tfile, (tstrip_t) strip, (tdata_t) buf, (tsize_t)-1);
      if(read == -1)
      {
          printf("Failed to read strip\n"); fflush(stdout);
      }
      if(read>ssize)
      {
          printf("Read to much\n"); fflush(stdout);
      }

      for(int64_t ii = 0; ii < read/sizeof(float); ii++) {
        size_t pos = ii+kk*nes + dd*perDirectory;
        V[pos] = buf[ii];
      }
    }
  }
  _TIFFfree(buf);
}

float raw_file_single_max(const char * rName, const size_t N)
{
  // Get max value from raw data file with N floats
//  printf("Getting max value from %s\n", rName);
  size_t buf_size = 1024*1024;
  float  max = -INFINITY;
  float * buf = malloc(buf_size*sizeof(float));
  FILE * fid = fopen(rName, "r");
  if(fid == NULL)
  {
    printf("ERROR: unable to open %s\n", rName);
  }
  size_t nread = 0;
//  printf("N = %zu\n", N);
  while(nread + buf_size < N)
  {
    size_t nr = fread(buf, buf_size*sizeof(float), 1, fid);
    assert(nr > 0);
    (void) nr;
    for(size_t kk = 0; kk<buf_size; kk++)
    {
      if(buf[kk] > max)
      {
        max = buf[kk];
      }
    }
    nread += buf_size;
  //  printf("%zu\n", nread);
  }
  size_t rem = N-nread;
  nread =  fread(buf, rem*sizeof(float), 1, fid);
  (void) nread;
  for(size_t kk = 0; kk<rem; kk++)
  {
    if(buf[kk] > max)
    {
      max = buf[kk];
    }
  }

  fclose(fid);
  free(buf);
  return max;
}

void uint16toraw(TIFF * tfile, const char * ofile,
    const uint32_t ssize,
    const uint32_t ndirs,
    const uint32_t nstrips,
    const uint32_t perDirectory
    )
{
  uint16_t * buf = _TIFFmalloc(ssize);
  float * wbuf = malloc(ssize/sizeof(uint16_t)*sizeof(float));
  FILE * fout = fopen(ofile, "w");

  assert(buf != NULL);

  for(int64_t dd=0; dd<ndirs; dd++) {
    TIFFSetDirectory(tfile, dd);
    //   printf("\r Directory %d / %u", dd+1, ndirs); fflush(stdout);
    for(int64_t kk=0; kk<nstrips; kk++) {
      tsize_t read = TIFFReadEncodedStrip(tfile, kk, buf, (tsize_t) - 1);
      assert(read>0);

      for(int64_t ii = 0; ii<read/sizeof(uint16_t); ii++) {
        wbuf[ii] = (float) buf[ii];
      }
      fwrite(wbuf, read/sizeof(uint16_t)*sizeof(float), 1, fout);
    }
  }
  //  printf("\n");
  fclose(fout);
  free(wbuf);
  _TIFFfree(buf);
}

void floattoraw(TIFF * tfile, const char * ofile,
    const uint32_t ssize,
    const uint32_t ndirs,
    const uint32_t nstrips,
    const uint32_t perDirectory
    )
{
  float * buf = _TIFFmalloc(ssize);
  float * wbuf = malloc(ssize/sizeof(float)*sizeof(float));
  FILE * fout = fopen(ofile, "w");

  assert(buf != NULL);

  for(int64_t dd=0; dd<ndirs; dd++) {
    TIFFSetDirectory(tfile, dd);
    for(int64_t kk=0; kk<nstrips; kk++) {
      int64_t strip = kk;
      tsize_t read = TIFFReadEncodedStrip(tfile, strip, buf, (tsize_t) - 1);
      assert(read>0);

      for(int64_t ii = 0; ii<read/sizeof(float); ii++) {
        wbuf[ii] = (float) buf[ii];
      }
      fwrite(wbuf, read, 1, fout);

    }
  }
  fclose(fout);
  free(wbuf);
  _TIFFfree(buf);
}


int fim_tiff_to_raw(const char * fName, const char * oName)
{
  // Convert a tif image, fName, to a raw float image, oName

  TIFF * tfile = TIFFOpen(fName, "r");

  if(tfile == NULL) {
    printf("Could not open %s\n", fName);
    return -1;
  }

  uint32_t M = 0, N = 0, BPS = 0, SF = 0, CMP=0;
  TIFFGetField(tfile, TIFFTAG_IMAGELENGTH, &N);
  TIFFGetField(tfile, TIFFTAG_IMAGEWIDTH, &M);
  TIFFGetField(tfile, TIFFTAG_BITSPERSAMPLE, &BPS);
  int gotSF = TIFFGetField(tfile, TIFFTAG_SAMPLEFORMAT, &SF);
  int gotCMP = TIFFGetField(tfile, TIFFTAG_COMPRESSION, &CMP);

  if(gotCMP)
  {
    if(CMP != 1)
    {
      printf("TIFFTAG_COMPRESSION=%u is not supported\n", CMP);
      return -1;
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
    printf("Only unsigned integer and float images are supported\n");
    return -1;
  }

  if(isUint && !(BPS == 16 || BPS == 8))
  {
    printf("Unsigned %d-bit images are not supported, only 8 and 16.\n", BPS);
    return -1;
  }

  if(isFloat && (BPS != 32))
  {
    printf("For floating point images, only 32-bit samples are supported %u\n", BPS);
  }

  // tiffio.h:typedef uint32 ttag_t;
  uint32_t PTAG = 0;
  int gotptag = TIFFGetField(tfile, TIFFTAG_PHOTOMETRIC, &PTAG);

  if(gotptag != 1)
  {
    printf("WARNING: Could not read TIFFTAG_PHOTOMETRIC, assuming minIsBlack\n");
    PTAG=1;
  }
  //printf("PTAG = %u\n", PTAG);
  if(PTAG == 0)
  {
    printf("Tiling mode does not support inverted images\n");
    return -1;
  }

  if(PTAG > 1)
  {
    printf("WARNING: Only WhiteIsZero or BlackIsZero are supported Photometric Interpretations, tag found: %u\n", PTAG);
    PTAG = 1;
    printf("Assuming min-is-black\n");
    fflush(stdout);
  }

  tmsize_t ssize = TIFFStripSize(tfile); // Seems to be in bytes
  uint32_t nstrips = TIFFNumberOfStrips(tfile);
  uint32_t ndirs = TIFFNumberOfDirectories(tfile);


  if(TIFFIsTiled(tfile))
  {
    printf("Tiled files are not supported\n");
    TIFFClose(tfile);
    return -1;
  }

  if(isFloat)
  {
    floattoraw(tfile, oName, ssize, ndirs, nstrips, M*N);
  }
  if(isUint)
  {
    uint16toraw(tfile, oName, ssize, ndirs, nstrips, M*N);
  }

  TIFFClose(tfile);

  return 0;
}

int fim_tiff_from_raw(const char * fName, // Name of tiff file to be written
    int64_t M, int64_t N, int64_t P, // Image dimensions
    const char * rName) // name of raw file
  /* Convert a raw float file to uint16 tif */
{

  size_t bytesPerSample = sizeof(uint16_t);
  char formatString[4] = "w";
  size_t MNP = (size_t) M * (size_t) N * (size_t ) P;
  if(MNP*sizeof(uint16) >= pow(2, 32))
  {
    sprintf(formatString, "w8\n");
    printf("WARNING: File is > 2 GB, using BigTIFF format\n");
  }

  TIFF * out = TIFFOpen(fName, formatString);
 // printf("Opened %s for writing using formatString: %s\n", fName, formatString);

  if(out == NULL)
  {
    printf("ERROR: Failed to open %s for writing using formatString: %s\n", fName, formatString);
    exit(1);
  }

  size_t linbytes = M*bytesPerSample;
  uint16_t * buf = _TIFFmalloc(linbytes);
  float * rbuf = malloc(M*sizeof(float));
  memset(buf, 0, linbytes);

  // Determine max value

  float scaling = 1;
  float rawmax = raw_file_single_max(rName, MNP);
  if(rawmax > 0)
  {
    scaling = 65535/rawmax;
  }
//  printf("Max value of file: %f\n", rawmax);

  FILE * rf = fopen(rName, "r");
  if(rf == NULL)
  {
    printf("Failed to open %s for writing\n", rName);
  }

// Can TIFFCheckpointDirectory be used to speed up this?
// http://maptools-org.996276.n3.nabble.com/help-writing-thumbnails-to-TIFF-file-td3824.html
// No, that isn't the trick to use.
// The "fast" writers puts the metadata last, that should be the way to go...

  for(size_t dd = 0; dd<P; dd++)
  {

    TIFFSetField(out, TIFFTAG_IMAGEWIDTH, M);  // set the width of the image
    TIFFSetField(out, TIFFTAG_IMAGELENGTH, N);    // set the height of the image
    TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1);   // set number of channels per pixel
    TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8*bytesPerSample);    // set the size of the channels
    TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image.
    //   Some other essential fields to set that you do not have to understand for now.
    TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);

    /* We are writing single page of the multipage file */
    TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
    /* Set the page number */
    TIFFSetField(out, TIFFTAG_PAGENUMBER, dd, P);


    for(size_t nn = 0; nn<N; nn++)
    {
      size_t nread = fread(rbuf, M*sizeof(float), 1, rf);
      (void)(nread);

      for(size_t mm = 0; mm<M; mm++)
      {
        buf[mm] = (uint16_t) (rbuf[mm]*scaling);
      }

      int ok = TIFFWriteScanline(out, // TIFF
          buf,
          nn, // row
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
  fclose(rf);
  free(rbuf);
  return 0;
}

int fim_tiff_write_zeros(const char * fName, int64_t N, int64_t M, int64_t P, FILE * fout)
{
  return fim_tiff_write(fName, NULL, N, M, P, fout);
}

int fim_tiff_write_float(const char * fName, const afloat * V,
    int64_t N, int64_t M, int64_t P, FILE * fout)
{

 fprintf(fout, "scaling: 1\n");
  size_t bytesPerSample = sizeof(float);
  char formatString[4] = "w";
  if(M*N*P*sizeof(uint16) >= pow(2, 32))
  {
    sprintf(formatString, "w8\n");
    fprintf(fout, "WARNING: File is > 2 GB, using BigTIFF format\n");
  }

  TIFF* out = TIFFOpen(fName, formatString);
  assert(out != NULL);

  size_t linbytes = (M+N)*bytesPerSample;
  float * buf = _TIFFmalloc(linbytes);
  memset(buf, 0, linbytes);

  for(size_t dd = 0; dd<P; dd++)
  {

    TIFFSetField(out, TIFFTAG_IMAGEWIDTH, N);  // set the width of the image
    TIFFSetField(out, TIFFTAG_IMAGELENGTH, M);    // set the height of the image
    TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1);   // set number of channels per pixel
    TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8*bytesPerSample);    // set the size of the channels
    TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image.
    //   Some other essential fields to set that you do not have to understand for now.
    TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);

    /* We are writing single page of the multipage file */
    TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
    /* Set the page number */
    TIFFSetField(out, TIFFTAG_PAGENUMBER, dd, P);


    for(size_t kk = 0; kk<M; kk++)
    {
      if(V != NULL)
      {
        for(size_t ll = 0; ll<N; ll++)
        {
          float value = V[M*N*dd + kk*N + ll];
          if(!isfinite(value))
          { value = 0; }

          buf[ll] = value;
        }
      }

      int ok = TIFFWriteScanline(out, // TIFF
          buf,
          kk, // row
          0); //sample
      if(ok != 1)
      {
        fprintf(fout, "TIFFWriteScanline failed\n");
      }
    }

    TIFFWriteDirectory(out);
  }

  _TIFFfree(buf);

  TIFFClose(out);
  return 0;
}

int fim_tiff_write(const char * fName, const afloat * V,
    int64_t N, int64_t M, int64_t P, FILE * fout)
{
  // if V == NULL and empty file will be written

  float scaling = 1;

  if(V != NULL)
  {
    float imax = 0;
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
    scaling = 1.0/imax*(pow(2,16)-1.0);
  }

  if(!isfinite(scaling))
  {
    fprintf(fout, "Non-finite scaling value, changing to 1\n");
    scaling = 1;
  }
  fprintf(fout, "scaling: %f\n", scaling);

  size_t bytesPerSample = sizeof(uint16_t);
  char formatString[4] = "w";
  if(M*N*P*sizeof(uint16) >= pow(2, 32))
  {
    sprintf(formatString, "w8\n");
    fprintf(fout, "WARNING: File is > 2 GB, using BigTIFF format\n");
  }

  TIFF* out = TIFFOpen(fName, formatString);
  assert(out != NULL);

  size_t linbytes = (M+N)*bytesPerSample;
  uint16_t * buf = _TIFFmalloc(linbytes);
  memset(buf, 0, linbytes);

  for(size_t dd = 0; dd<P; dd++)
  {

    TIFFSetField(out, TIFFTAG_IMAGEWIDTH, N);  // set the width of the image
    TIFFSetField(out, TIFFTAG_IMAGELENGTH, M);    // set the height of the image
    TIFFSetField(out, TIFFTAG_SAMPLESPERPIXEL, 1);   // set number of channels per pixel
    TIFFSetField(out, TIFFTAG_BITSPERSAMPLE, 8*bytesPerSample);    // set the size of the channels
    TIFFSetField(out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);    // set the origin of the image.
    //   Some other essential fields to set that you do not have to understand for now.
    TIFFSetField(out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(out, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
    // TODO TIFFSSetFieldTIFFTAG_SOFTWARE

    /* We are writing single page of the multipage file */
    TIFFSetField(out, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
    /* Set the page number */
    TIFFSetField(out, TIFFTAG_PAGENUMBER, dd, P);


    for(size_t kk = 0; kk<M; kk++)
    {
      if(V != NULL)
      {
        for(size_t ll = 0; ll<N; ll++)
        {
          float value = V[M*N*dd + kk*N + ll]*scaling;
          if(!isfinite(value))
          { value = 0; }

          buf[ll] = (uint16_t) round(value);
        }
      }

      int ok = TIFFWriteScanline(out, // TIFF
          buf,
          kk, // row
          0); //sample
      if(ok != 1)
      {
        fprintf(fout, "TIFFWriteScanline failed\n");
      }
    }

    TIFFWriteDirectory(out);
  }

  _TIFFfree(buf);

  TIFFClose(out);
  return 0;
}

int fim_tiff_get_size(char * fname, int64_t * M, int64_t * N, int64_t * P)
{
  TIFF * tiff = TIFFOpen(fname, "r");
  uint32_t m, n, p;

  if(tiff == NULL) {
    return -1;
  }

  int ok = 1;
  ok *= TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &m);
  ok *= TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &n);
  if(ok != 1)
  {
    return -1;
    TIFFClose(tiff);
  }

  p = TIFFNumberOfDirectories(tiff);
  TIFFClose(tiff);

  M[0] = m; N[0] = n; P[0] = p;
  return 0;
}


afloat * fim_tiff_read(const char * fName,
    int64_t * N0, int64_t * M0, int64_t * P0, int verbosity)
{
  return fim_tiff_read_sub(fName, N0, M0, P0, verbosity,
      0, // sub disabled
      0,0,0, // start
      0,0,0); // width
}

afloat * fim_tiff_read_sub(const char * fName,
    int64_t * M0, int64_t * N0, int64_t * P0, int verbosity,
    int subregion,
    int64_t sM, int64_t sN, int64_t sP, // start
    int64_t wM, int64_t wN, int64_t wP) // width
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
  uint32_t _M = 0, _N = 0, BPS = 0, SF = 0, CMP=0;
  TIFFGetField(tfile, TIFFTAG_IMAGELENGTH, &_N);
  TIFFGetField(tfile, TIFFTAG_IMAGEWIDTH, &_M);
  TIFFGetField(tfile, TIFFTAG_BITSPERSAMPLE, &BPS);
  int gotSF = TIFFGetField(tfile, TIFFTAG_SAMPLEFORMAT, &SF);
  int gotCMP = TIFFGetField(tfile, TIFFTAG_COMPRESSION, &CMP);
  int64_t M = _M, N = _N;

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
    printf("Only unsigned integer and float images are supported\n");
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

  // tiffio.h:typedef uint32 ttag_t;
  uint32_t PTAG = 0;
  int gotptag = TIFFGetField(tfile, TIFFTAG_PHOTOMETRIC, &PTAG);

  int inverted = 0;
  if(gotptag != 1)
  {
    printf("WARNING: Could not read TIFFTAG_PHOTOMETRIC, assuming minIsBlack\n");
    PTAG=1;
  }
  //printf("PTAG = %u\n", PTAG);
  if(PTAG == 0)
  {
    inverted = 1;
  }

  if(PTAG > 1)
  {
    printf("WARNING: Only WhiteIsZero or BlackIsZero are supported Photometric Interpretations, tag found: %u\n", PTAG);
    PTAG = 1;
    printf("Assuming min-is-black\n");
    fflush(stdout);
  }

  tmsize_t ssize = TIFFStripSize(tfile); // Seems to be in bytes
  uint32_t nstrips = TIFFNumberOfStrips(tfile);
  uint32_t ndirs = TIFFNumberOfDirectories(tfile);
  P0[0] = (size_t) ndirs;
  int64_t P = P0[0];

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
  afloat * V = NULL;
  size_t nel = M*N*ndirs;
  if(subregion)
  {
    nel = wM*wN*wP;
  }

  V = fftwf_malloc(nel*sizeof(float));
  memset(V, 0, nel*sizeof(float));

  if(isFloat)
  {
    if(verbosity > 1)
    {
      printf("ReadFloat ...\n");
    }
    if(subregion)
    {
      readFloat_sub(tfile, V, ssize, ndirs, nstrips, M*N,
          sM, sN, sP, wM, wN, wP);
    } else {
      readFloat(tfile, V, ssize, ndirs, nstrips, M*N);
    }
  }
  if(isUint)
  {
    if(verbosity > 1)
    {
      printf("ReadUint ...\n");
    }
    if(BPS == 16)
    {
      if(subregion)
      {
        printf("%" PRId64 " %" PRId64 " %" PRId64 ", %" PRId64 " %" PRId64 " %" PRId64 ", %" PRId64 " %" PRId64 " %" PRId64 "\n", M, N, P, sM, sN, sP, wM, wN, wP);
        readUint16_sub(tfile, V, ssize, ndirs, nstrips, M*N,
            M, N, (int) P, sM, sN, sP, wM, wN, wP);
      } else {
        readUint16(tfile, V, ssize, ndirs, nstrips, M*N);
      }
    }
    if(BPS == 8)
    {
      if(subregion)
      {
        readUint8_sub(tfile, V, ssize, ndirs, nstrips, M*N,
            sM, sN, sP, wM, wN, wP);
      } else {
        readUint8(tfile, V, ssize, ndirs, nstrips, M*N);
      }
    }
  }

  TIFFClose(tfile);
  if(verbosity>1)
  {printf("Done reading\n");
  }

  if(inverted == 1)
  {
    fim_invert(V, M*N*P);
  }

  return V;
}

#ifdef unittest
int main(int argc, char ** argv)
{

    printf("This test will read a file from command line, normalize it and write back to 'foo.tif'\n");
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

  int64_t M = 0, N = 0, P = 0;

  afloat * I = (afloat *) fim_tiff_read(inname, &M, &N, &P, 1);

  if(I == NULL)
  {
    printf("Failed to read %s\n", inname);
    exit(1);
  }

  floatimage_show_stats(I, M, N, P);
  // floatimage_normalize(I, M*N*P);

  fim_tiff_write(outname, I, M, N, P, stdout);

  free(I);

  return 0;
}
#endif


char * tiff_is_supported(TIFF * tiff)
{
  char * errStr = malloc(1024);
  if(tiff == NULL) {
    sprintf(errStr, "Can't be opened!");
    return errStr;
  }

  // Tags: ImageWidth and ImageLength
  uint32_t _M = 0, _N = 0, BPS = 0, SF = 0, CMP=0;
  TIFFGetField(tiff, TIFFTAG_IMAGELENGTH, &_N);
  TIFFGetField(tiff, TIFFTAG_IMAGEWIDTH, &_M);
  TIFFGetField(tiff, TIFFTAG_BITSPERSAMPLE, &BPS);
  int gotSF = TIFFGetField(tiff, TIFFTAG_SAMPLEFORMAT, &SF);
  int gotCMP = TIFFGetField(tiff, TIFFTAG_COMPRESSION, &CMP);

  if(gotCMP)
  {
    if(CMP != 1)
    {
      sprintf(errStr, "TIFFTAG_COMPRESSION=%u is not supported\n", CMP);
      return errStr;
    }
  }

  int isUint = 0;
  int isFloat = 0;

  if(gotSF)
  {
    int okFormat = 0;
    if( SF == SAMPLEFORMAT_UINT)
    {
      isUint = 1;
      okFormat = 1;
    }
    if( SF == SAMPLEFORMAT_IEEEFP)
    {
      isFloat = 1;
      okFormat = 1;
    }
    if(okFormat == 0)
    {
      sprintf(errStr, "Neither SAMPLEFORMAT_UINT or SAMPLEFORMAT_IEEEFP\n");
      return errStr;
    }
  }
  else {
    printf("Warning: TIFFTAG_SAMPLEFORMAT not specified, assuming uint but that could be wrong!\n");
    isUint = 1;
  }

  if(isFloat && !(BPS == 32))
  {
    sprintf(errStr, "Float %d-bit images are not supported, only 32-bit.\n", BPS);
    return errStr;
  }

  if(isUint && !(BPS == 16))
  {
    sprintf(errStr, "Unsigned %d-bit images are not supported, only 16-bit.\n", BPS);
    return errStr;
  }

  // tiffio.h:typedef uint32 ttag_t;
  uint32_t PTAG = 0;
  int gotptag = TIFFGetField(tiff, TIFFTAG_PHOTOMETRIC, &PTAG);

  if(gotptag != 1)
  {
    printf("WARNING: Could not read TIFFTAG_PHOTOMETRIC, assuming minIsBlack\n");
    PTAG=1;
  }

  if(PTAG != 1)
  {
    sprintf(errStr, "Only BlackIsZero are supported Photometric Interpretations, tag found: %u\n", PTAG);
    return errStr;
  }

  free(errStr);
  return NULL;
}

int fim_tiff_maxproj(char * in, char * out)
{
  int64_t M, N, P;

  if(fim_tiff_get_size(in, &M, &N, &P))
  {
    printf("Can't open %s to get image dimension\n", in);
    return -1;
  }

  TIFF * input = TIFFOpen(in, "r");
  char * errStr = tiff_is_supported(input);
  if(errStr != NULL)
  {
    printf("Can't process %s\n", in);
    printf("Error: %s\n", errStr);
    return -1;
  }

  uint32_t SF;
  int gotSF = TIFFGetField(input, TIFFTAG_SAMPLEFORMAT, &SF);
  if(gotSF != 1)
  {
    printf("Unable to determine the sample format of %s\n", in);
    return 1;
  }

  TIFF * output = TIFFOpen(out, "w");
  TIFFSetField(output, TIFFTAG_IMAGEWIDTH, M);
  TIFFSetField(output, TIFFTAG_IMAGELENGTH, N);
  TIFFSetField(output, TIFFTAG_SAMPLESPERPIXEL, 1);
  if(SF == SAMPLEFORMAT_UINT)
  {
    TIFFSetField(output, TIFFTAG_BITSPERSAMPLE, 16);
    TIFFSetField(output, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_UINT);
  }

  if(SF == SAMPLEFORMAT_IEEEFP)
  {
    TIFFSetField(output, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField(output, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP);
  }

  TIFFSetField(output, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
  TIFFSetField(output, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField(output, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
  TIFFSetField(output, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
  TIFFSetField(output, TIFFTAG_PAGENUMBER, 1, 1);

  tmsize_t ssize = TIFFStripSize(input); // Seems to be in bytes
//  printf("Strip size: %zu b\n", (size_t) ssize);
  uint32_t nstrips = TIFFNumberOfStrips(input);



  // Input is 16 bit unsigned int.
  if(SF == SAMPLEFORMAT_UINT)
  {
  uint16_t * mstrip = _TIFFmalloc(ssize); // For max over all directories
  uint16_t * strip    = _TIFFmalloc(ssize);

  for(int64_t nn = 0; nn<nstrips; nn++) // Each strip
  {
    memset(mstrip, 0, ssize);
    tsize_t read = 0;
    for(int64_t pp = 0; pp<P; pp++) // For each directory
    {
      TIFFSetDirectory(input, pp); // Does it keep the location for each directory?
      read = TIFFReadEncodedStrip(input, nn, strip, (tsize_t)-1);
      for(int64_t kk = 0; kk<read/2; kk++)
      {
        if(strip[kk] > mstrip[kk])
        {
          mstrip[kk] = strip[kk];
        }
      }
    }
    TIFFWriteRawStrip(output, 0, mstrip, read);
  }

  _TIFFfree(strip);
  _TIFFfree(mstrip);
  }

  // Input is 32-bit float
  if(SF == SAMPLEFORMAT_IEEEFP)
  {
  float * mstrip = _TIFFmalloc(ssize); // For max over all directories
  float * strip    = _TIFFmalloc(ssize);

  for(int64_t nn = 0; nn<nstrips; nn++) // Each strip
  {
    memset(mstrip, 0, ssize);
    tsize_t read = 0;
    for(int64_t pp = 0; pp<P; pp++) // For each directory
    {
      TIFFSetDirectory(input, pp); // Does it keep the location for each directory?
      read = TIFFReadEncodedStrip(input, nn, strip, (tsize_t)-1);
      for(int64_t kk = 0; kk<read/4; kk++)
      {
        if(strip[kk] > mstrip[kk])
        {
          mstrip[kk] = strip[kk];
        }
      }
    }
    TIFFWriteRawStrip(output, 0, mstrip, read);
  }
  _TIFFfree(strip);
  _TIFFfree(mstrip);

  }


  TIFFClose(input);
  TIFFClose(output);

  return 0;
}
