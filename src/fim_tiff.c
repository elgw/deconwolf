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
 * From:
 * https://www.cs.rochester.edu/u/nelson/courses/vision/resources/tiff/libtiff.html#Errors
 * Finally, note that the last strip of data in an image may have
 * fewer rows in it than specified by the RowsPerStrip tag. A reader
 * should not assume that each decoded strip contains a full set of
 * rows in it.
 */

/*
 * Forward declarations
 */

/* Read a uint16 volumetric image */
void readUint16(TIFF * tfile, float * V,
                const tsize_t ssize, // strip size in bytes
                const uint32_t ndirs, // number of directories (z-planes)
                const uint32_t nstrips, // number of strips
                const uint32_t perDirectory); // elements per plane


int fim_tiff_write_opt(const char * fName, const float * V,
                       ttags * T,
                       int64_t N, int64_t M, int64_t P, int scaling);

FILE * fim_tiff_log = NULL; /* Indicates stdout */
void fim_tiff_init(void)
{
    fim_tiff_log = stdout;
}

/* Used to redirect errors from libtiff */
void tiffErrHandler(const char* module, const char* fmt, va_list ap)
{
    fprintf(fim_tiff_log, "libtiff: Module: %s\n", module);
    fprintf(fim_tiff_log, "libtiff: ");
    vfprintf(fim_tiff_log, fmt, ap);
    fprintf(fim_tiff_log, "\n");
}


void fim_tiff_set_log(FILE * fp)
{
    if(fp != NULL)
    {
        fim_tiff_log = fp;
    } else {
        fim_tiff_log = stdout;
    }
    TIFFSetWarningHandler(tiffErrHandler);
    TIFFSetErrorHandler(tiffErrHandler);
}

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
    float * im = fim_zeros(M*N*P);

    size_t pos1 = 1111;
    size_t pos2 = 2222;
    size_t pos3 = 3333;

    im[pos1] = 1;
    im[pos2] = 2;
    im[pos3] = 3;

    fim_tiff_write(fname, im, NULL, M, N, P);

    int64_t M2 = 0, N2 = 0, P2 = 0;
    float * im2 = fim_tiff_read(fname, NULL, &M2, &N2, &P2, 0);
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


void floatimage_normalize(float * restrict I, const size_t N)
{
    // Scale image to span the whole 16 bit range.
    float imax = I[0]; float imin = I[0];
    int ok = 1;
//#pragma omp parallel for reduction(min:imin) reduction(max:imax) shared(I)
    for(size_t kk=0; kk<N; kk++)
    {
        if(!isfinite(I[kk]))
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
        fim_mult_scalar(I, N, (pow(2,16)-2)/imax);
    }

}

void floatimage_show_stats(float * I, size_t N, size_t M, size_t P)
{
    float isum = 0;
    float imin = I[0];
    float imax = I[0];
    for(size_t kk = 0; kk<M*N*P; kk++)
    {
        I[kk] > imax ? imax = I[kk] : 0 ;
        I[kk] < imin ? imin = I[kk] : 0 ;
        isum += I[kk];
    }
    printf("min: %f max: %f mean: %f\n", imin, imax, isum/(M*N*P));
}


INLINED extern void sub2ind(size_t ind,
                            int64_t M, int64_t N, __attribute__((unused)) int64_t P,
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

void readUint16_sub(TIFF * tfile, float * V,
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

    for(size_t dd=0; dd < ndirs; dd++) {
        TIFFSetDirectory(tfile, dd);
        for(size_t kk=0; kk<nstrips; kk++) {
            size_t strip = kk;
            tsize_t read = TIFFReadEncodedStrip(tfile, strip, buf, (tsize_t) - 1);
            assert(read>0);

            for(size_t ii = 0; ii < read/sizeof(uint16_t); ii++) {
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

void readUint8_sub(TIFF * tfile, float * V,
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
    size_t nes = ssize/sizeof(uint8_t);
    uint8_t * buf = _TIFFmalloc(ssize);
    assert(buf != NULL);

    for(size_t dd=0; dd < ndirs; dd++) {
        TIFFSetDirectory(tfile, dd);
        for(size_t kk=0; kk<nstrips; kk++) {
            size_t strip = kk;
            tsize_t read = TIFFReadEncodedStrip(tfile, strip, buf, (tsize_t) - 1);
            assert(read>0);

            for(size_t ii = 0; ii < read/sizeof(uint8_t); ii++) {
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
                const tsize_t ssize,
                const uint32_t ndirs,
                const uint32_t nstrips,
                const uint32_t perDirectory
    )
{
//    printf("readUint16\n"); fflush(stdout);
    // Number of elements per strip
    size_t nes = ssize/sizeof(uint16_t);
    //  uint16_t * buf = _TIFFmalloc(ssize);
    //uint16_t * buf = malloc(ssize);

    uint16_t * buf = _TIFFmalloc(TIFFStripSize(tfile));
    if(buf == NULL)
    {
        printf("Failed to allocate %" PRId64 " bytes of memory!", ssize);
        fflush(stdout);
        return;
    }

    for(int64_t dd=0; dd<ndirs; dd++) {
        int ok = TIFFSetDirectory(tfile, dd);
        if(ok == 0)
        {
            printf("Failed to choose directory %" PRIu64 "\n", dd);
            printf("Either the file is corrupt or this program has a bug. "
                   "Please verify with another program that the file is ok before "
                   "filing a bug report.\n");
            exit(EXIT_FAILURE);
            fflush(stdout);
        }


        size_t nread = 0;

        for(int64_t ss=0; ss < TIFFNumberOfStrips(tfile); ss++)
        {
            //printf("Reading strip %ld of size %ld B\n", ss, ssize); fflush(stdout);
            // vs TIFFReadEncodedStrip vs TIFFReadRawStrip ?
            tsize_t read = TIFFReadRawStrip(tfile,
                                            (tstrip_t) ss, // Strip number
                                            (tdata_t) buf, // target
                                            (tsize_t) -1); // The entire strip
            //printf("Got %ld\n", read);
            if(read == -1)
            {
                fprintf(stderr, "Failed to read from tif file and can't continue.\n");
                fprintf(stderr, "ERROR at FILE: %s FUNCTION: %s LINE: %d\n",
                        __FILE__, __FUNCTION__, __LINE__);
                tmsize_t ssize2 = TIFFStripSize(tfile);
                printf("ssize2: %ld\n", ssize2);
                fflush(stdout);
                exit(EXIT_FAILURE);
            }

            //printf("Will write %d elements from strip %d\n", read/sizeof(uint16_t), ss);
            for(size_t ii = 0; ii < read/sizeof(uint16_t); ii++)
            {
                //size_t pos = (size_t)( ii+ss*nes + dd*perDirectory);
                size_t pos = (size_t)( nread/sizeof(uint16_t) + ii + dd*perDirectory);
                //printf("pos=%zu\n", pos);
                V[pos] = buf[ii];
            }
            nread += read;
        } // strip
    } // directory

    _TIFFfree(buf);
    return;
}

void readUint8(TIFF * tfile, float * V,
               const tsize_t ssize,
               const uint32_t ndirs,
               const uint32_t nstrips,
               const uint32_t perDirectory
    )
{
    // Number of elements per strip
    size_t nes = ssize/sizeof(uint8_t);
    //  uint16_t * buf = _TIFFmalloc(ssize);
    uint8_t * buf = malloc(ssize);

    if(buf == NULL)
    {
        printf("Failed to allocate %" PRId64 " bytes of memory!", ssize);
        fflush(stdout);
        return;
    }

    for(int64_t dd=0; dd<ndirs; dd++) {
        int ok = TIFFSetDirectory(tfile, dd);
        if(ok == 0)
        {
            printf("Failed to choose directory %" PRIu64 "\n", dd);
            printf("Either the file is corrupt or this program has a bug. "
                   "Please verify with another program that the file is ok before "
                   "filing a bug report.\n");
            exit(EXIT_FAILURE);
            fflush(stdout);
        }
        for(int64_t kk=0; kk<nstrips; kk++) {
            int64_t strip = kk;
            tsize_t read = TIFFReadEncodedStrip(tfile,
                                                (tstrip_t) strip, // Strip number
                                                (tdata_t) buf, // target
                                                ssize); // The entire strip
            if(read == -1)
            {
                printf("Failed to read from tif file\n");
                fflush(stdout);
            }
            if(read > ssize)
            {
                printf("Read to much!\n"); fflush(stdout);
            }
            //printf("\r%lu/%lu", dd, kk); fflush(stdout);
            for(size_t ii = 0; ii < read/sizeof(uint8_t); ii++) {
                size_t pos = (size_t)( ii+kk*nes + dd*perDirectory);
                //printf("pos=%zu\n", pos);
                V[pos] = buf[ii];
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

            for(size_t ii = 0; ii < read/sizeof(float); ii++) {
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
                 const uint32_t nstrips)
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

            for(size_t ii = 0; ii < read/sizeof(uint16_t); ii++) {
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
                const uint32_t nstrips)
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

            for(size_t ii = 0; ii < read/sizeof(float); ii++) {
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
        fprintf(fim_tiff_log, "Warning: TIFFTAG_SAMPLEFORMAT not specified, assuming uint but that could be wrong!\n");
        isUint = 1;
    }

    if(!(isUint || isFloat))
    {
        fprintf(stderr, "Only unsigned integer and float images are supported\n");
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
        fprintf(fim_tiff_log, "fim_tiff: WARNING: Could not read TIFFTAG_PHOTOMETRIC, assuming minIsBlack\n");
        PTAG=1;
    }
    //printf("PTAG = %u\n", PTAG);
    if(PTAG == 0)
    {
        fprintf(stderr, "Tiling mode does not support inverted images\n");
        return -1;
    }

    if(PTAG > 1)
    {
        fprintf(fim_tiff_log, "fim_tiff: WARNING: Only WhiteIsZero or BlackIsZero are supported Photometric Interpretations, tag found: %u\n", PTAG);
        PTAG = 1;
        fprintf(fim_tiff_log, "fim_tiff: Assuming min-is-black\n");
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
        floattoraw(tfile, oName, ssize, ndirs, nstrips);
    }
    if(isUint)
    {
        uint16toraw(tfile, oName, ssize, ndirs, nstrips);
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
    if(MNP*sizeof(uint16_t) >= pow(2, 32))
    {
        sprintf(formatString, "w8\n");
        fprintf(fim_tiff_log, "fim_tiff WARNING: File is > 2 GB, using BigTIFF format\n");
    }

    TIFF * out = TIFFOpen(fName, formatString);
    // printf("Opened %s for writing using formatString: %s\n", fName, formatString);

    if(out == NULL)
    {
        fprintf(stderr, "fim_tiff ERROR: Failed to open %s for writing using formatString: %s\n", fName, formatString);
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
        fprintf(stderr, "fim_tiff ERROR: Failed to open %s for writing\n", rName);
        exit(EXIT_FAILURE);
    }

    // Can TIFFCheckpointDirectory be used to speed up this?
    // http://maptools-org.996276.n3.nabble.com/help-writing-thumbnails-to-TIFF-file-td3824.html
    // No, that isn't the trick to use.
    // The "fast" writers puts the metadata last, that should be the way to go...

    for(size_t dd = 0; dd < (size_t) P; dd++)
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


        for(size_t nn = 0; nn < (size_t) N; nn++)
        {
            size_t nread = fread(rbuf, M*sizeof(float), 1, rf);
            (void)(nread);

            for(size_t mm = 0; mm < (size_t) M; mm++)
            {
                buf[mm] = (uint16_t) (rbuf[mm]*scaling);
            }

            int ok = TIFFWriteScanline(out, // TIFF
                                       buf,
                                       nn, // row
                                       0); //sample
            if(ok != 1)
            {
                fprintf(stderr, "TIFFWriteScanline failed\n");
                exit(EXIT_FAILURE);
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

int fim_tiff_write_zeros(const char * fName, int64_t N, int64_t M, int64_t P)
{
    return fim_tiff_write(fName, NULL, NULL, N, M, P);
}

int fim_tiff_write_float(const char * fName, const float * V,
                         ttags * T,
                         int64_t N, int64_t M, int64_t P)
{

    if(fim_tiff_log == NULL)
    {
        fim_tiff_log = stdout;
    }
    size_t bytesPerSample = sizeof(float);
    char formatString[4] = "w";
    if(M*N*P*sizeof(uint16_t) >= pow(2, 32))
    {
        sprintf(formatString, "w8\n");
        fprintf(fim_tiff_log, "fim_tiff: File is > 2 GB, using BigTIFF format\n");
    }

    TIFF* out = TIFFOpen(fName, formatString);
    if(T != NULL)
    {
        ttags_set(out, T);
    }
    assert(out != NULL);

    size_t linbytes = (M+N)*bytesPerSample;
    float * buf = _TIFFmalloc(linbytes);
    memset(buf, 0, linbytes);

    for(size_t dd = 0; dd < (size_t) P; dd++)
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


        for(size_t kk = 0; kk < (size_t) M; kk++)
        {
            if(V != NULL)
            {
                for(size_t ll = 0; ll < (size_t) N; ll++)
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
                fprintf(stderr, "TIFFWriteScanline failed\n");
                exit(EXIT_FAILURE);
            }
        }

        TIFFWriteDirectory(out);
    }

    _TIFFfree(buf);

    TIFFClose(out);
    return 0;
}

int fim_tiff_write(const char * fName, const float * V,
                   ttags * T,
                   int64_t N, int64_t M, int64_t P)
{
    /* Default, use scaling */
    return fim_tiff_write_opt(fName, V, T, N, M, P, 1);
}

int fim_tiff_write_noscale(const char * fName, const float * V,
                           ttags * T,
                           int64_t N, int64_t M, int64_t P)
{
    /* Default, use scaling */
    return fim_tiff_write_opt(fName, V, T, N, M, P, 0);
}


int fim_tiff_write_opt(const char * fName, const float * V,
                       ttags * T,
                       int64_t N, int64_t M, int64_t P, int scale)
{
    if(fim_tiff_log == NULL)
    {
        fim_tiff_log = stdout;
    }
    // if V == NULL and empty file will be written

    float scaling = 1;

    if(scale)
    {
        if(V != NULL)
        {
            float imax = fim_max(V, M*N*P);
            scaling = 1.0/imax*(pow(2,16)-2.0);
        }
    }

    if(!isfinite(scaling))
    {
        fprintf(fim_tiff_log, "fim_tiff WARNING: Non-finite scaling value, changing to 1\n");
        scaling = 1;
    }
    fprintf(fim_tiff_log, "scaling: %f\n", scaling);

    size_t bytesPerSample = sizeof(uint16_t);
    char formatString[4] = "w";
    if(M*N*P*sizeof(uint16_t) >= pow(2, 32))
    {
        sprintf(formatString, "w8\n");
        fprintf(fim_tiff_log, "tim_tiff: File is > 2 GB, using BigTIFF format\n");
    }

    TIFF* out = TIFFOpen(fName, formatString);
    assert(out != NULL);

    if(T)
    {
        ttags_set(out, T);
    }

    size_t linbytes = (M+N)*bytesPerSample;
    uint16_t * buf = _TIFFmalloc(linbytes);
    memset(buf, 0, linbytes);

    for(size_t dd = 0; dd < (size_t) P; dd++)
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


        for(size_t kk = 0; kk < (size_t) M; kk++)
        {
            if(V != NULL)
            {
                for(size_t ll = 0; ll < (size_t) N; ll++)
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
                fprintf(stderr, "fim_tiff ERROR: TIFFWriteScanline failed\n");
                exit(EXIT_FAILURE);
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


float * fim_tiff_read(const char * fName,
                      ttags * T,
                      int64_t * N0, int64_t * M0, int64_t * P0, int verbosity)
{
    return fim_tiff_read_sub(fName, T, N0, M0, P0, verbosity,
                             0, // sub disabled
                             0,0,0, // start
                             0,0,0); // width
}


ttags * ttags_new()
{
    ttags * T  = malloc(sizeof(ttags));
    T->xresolution = 1;
    T->yresolution = 1;
    T->zresolution = 1;
    T->imagedescription = NULL;
    T->software = NULL;
    T-> resolutionunit = RESUNIT_NONE;
    T->IJIJinfo = NULL;
    T->nIJIJinfo = 0;
    // Image size MxNxP
    T->M = 0;
    T->N = 0;
    T->P = 0;
    return T;
}

void ttags_set_imagesize(ttags * T, int M, int N, int P)
{
    T->M = M;
    T->N = N;
    T->P = P;
}

void ttags_set_pixelsize(ttags * T, double xres, double yres, double zres)
{
    if(T->M == 0)
    {
        fprintf(stderr, "use ttags_set_imagesize before ttags_set_pixelsize");
        exit(EXIT_FAILURE);
    }
    if(xres != yres)
    {
        fprintf(stderr, "Only supports isotropic pixels in x-y\n");
        exit(EXIT_FAILURE);
    }
    T->xresolution = 1/xres; // Pixels per nm
    T->yresolution = 1/yres;
    T->zresolution = zres; // nm
    if(T->imagedescription)
    {
        free(T->imagedescription);
    }
    T->imagedescription = malloc(1024);
    //sprintf(T->imagedescription, "ImageJ=1.11a.images=%d.slices=%d.hyperstack=true.mode=grayscale.unit=nm.spacing=%.3f.", T->P, T->P, T->zresolution);
    sprintf(T->imagedescription, "ImageJ=1.52r\nimages=%d\nslices=%d\nunit=nm\nspacing=%.1f\nloop=false.", T->P, T->P, T->zresolution);
}

void ttags_free(ttags ** Tp)
{
    ttags * T = Tp[0];
    if(T == NULL)
    {
        fprintf(stderr, "fim_tiff: T = NULL, on line %d in %s\n", __LINE__, __FILE__);
        exit(EXIT_FAILURE);
    }
    if(T->software != NULL)
    {
        free(T->software);
    }

    if(T->imagedescription != NULL)
    {
        free(T->imagedescription);
    }
    free(T);
}

void ttags_set_software(ttags * T, char * sw)
{
    if(T->software != NULL)
    {
        free(T->software);
    }
    T->software = malloc(strlen(sw)+2);
    sprintf(T->software, "%s", sw);
}

void ttags_get(TIFF * tfile, ttags * T)
{
    // https://docs.openmicroscopy.org/ome-model/5.6.3/ome-tiff/specification.html
    // a string of OME-XML metadata embedded in the ImageDescription tag of the first IFD (Image File Directory) of each file. The XML string must be UTF-8 encoded.

    T->xresolution = 1;
    T->yresolution = 1;
    T->imagedescription = NULL;
    T->software = NULL;
    T->resolutionunit = RESUNIT_NONE;


    uint16_t runit;
    if(TIFFGetField(tfile, TIFFTAG_RESOLUTIONUNIT, &runit))
    {
        T->resolutionunit = runit;
    }

    float xres = 0, yres = 0;
    if(TIFFGetField(tfile, TIFFTAG_XRESOLUTION, &xres))
    {
        T->xresolution = xres;
    }

    if(TIFFGetField(tfile, TIFFTAG_XRESOLUTION, &yres))
    {
        T->yresolution = yres;
    }

    char * desc = NULL;
    if(TIFFGetField(tfile, TIFFTAG_IMAGEDESCRIPTION, &desc) == 1)
    {
        T->imagedescription = malloc(strlen(desc)+2);
        strcpy(T->imagedescription, desc);
    }

    char * software = NULL;
    if(TIFFGetField(tfile, TIFFTAG_SOFTWARE, &software) == 1)
    {
        T->software = malloc(strlen(software)+2);
        strcpy(T->software, software);
        //    printf("! Got software tag: %s\n", T->software);
    }

#if 0
    /*
     * From https://github.com/imagej/ImageJA/blob/master/src/main/java/ij/io/TiffDecoder.java
     * 	public static final int META_DATA_BYTE_COUNTS = 50838; // private tag registered with Adobe
     *      public static final int META_DATA = 50839; // private tag registered with Adobe
     */

    uint32_t count;
    void * data;
    if(TIFFGetField(tfile, XTAG_IJIJUNKNOWN, &count, &data))
    {
        uint32_t * dvalue  = (uint32_t*) data;

        for(int kk = 0; kk<count; kk++)
        {
            fprintf(fim_tiff_log, "Tag %d: %d, count %d\n", XTAG_IJIJUNKNOWN, dvalue[kk], count);
        }
    }

    T->nIJIJinfo = 0;
    T->IJIJinfo = NULL;
    if(TIFFGetField(tfile, XTAG_IJIJINFO, &count, &data))
    {
        T->nIJIJinfo = count;
        T->IJIJinfo = malloc(count);
        memcpy(T->IJIJinfo, data, count);
        uint8_t * udata = (uint8_t*) data;
        for(int kk = 0; kk<count; kk++)
        { fprintf(fim_tiff_log, "%02d %c %u\n ", kk, udata[kk], udata[kk]);};
    }
#endif
}

void ttags_set(TIFF * tfile, ttags * T)
{
    TIFFSetDirectory(tfile, 0);
    if(T->software != NULL)
    {
        TIFFSetField(tfile, TIFFTAG_SOFTWARE, T->software);
    }

    if(T->imagedescription != NULL)
    {
        TIFFSetField(tfile, TIFFTAG_IMAGEDESCRIPTION, T->imagedescription);
    }

#if 0
    if(T->nIJIJinfo > 0)
    {
        static TIFFFieldInfo xtiffFieldInfo[] =
            {
                { XTAG_IJIJUNKNOWN,TIFF_VARIABLE, TIFF_VARIABLE, TIFF_LONG, FIELD_CUSTOM,
                  0, 1, "IJIJunknown" },
                { XTAG_IJIJINFO, TIFF_VARIABLE, TIFF_VARIABLE, TIFF_ASCII, FIELD_CUSTOM,
                  1, 0, "IJIJinfo" }
            };

        // Guessing that ImageJ requires IJIJinfo to be of type ASCII
        // Can't write the full string using ASCII since \0 are used
        // to determine the length of the data.
        // Using TIFF_BYTE confuses ImageJ
        // Can't get direct access with TIFFGetField; it returns a copy

        TIFFMergeFieldInfo(tfile, xtiffFieldInfo, 2*sizeof(TIFFFieldInfo));

        uint32_t  avalue[] = {12, 80};
        TIFFSetField(tfile, XTAG_IJIJUNKNOWN, 2, &avalue);
        fprintf(fim_tiff_log, "fim_tiff: -> nIJINFO: %d bytes", T->nIJIJinfo);
        TIFFSetField(tfile, XTAG_IJIJINFO, T->IJIJinfo, T->nIJIJinfo);

        char * readback = NULL;
        TIFFGetField(tfile, XTAG_IJIJINFO, &readback);
        fprintf(fim_tiff_log, "fim_tiff: readback: %s\n", readback);
    }
#endif

    TIFFSetField(tfile, TIFFTAG_XRESOLUTION, T->xresolution);
    TIFFSetField(tfile, TIFFTAG_YRESOLUTION, T->yresolution);
    TIFFSetField(tfile, TIFFTAG_RESOLUTIONUNIT, T->resolutionunit);
    return;
}

void ttags_show(FILE * fout, ttags* T)
{

    fprintf(fout, "Resolution unit: ");
    switch(T->resolutionunit)
    {
    case RESUNIT_NONE:
        fprintf(fout, "RESUNIT_NONE\n");
        break;
    case RESUNIT_INCH:
        fprintf(fout, "RESUNIT_INCH\n");
        break;
    case RESUNIT_CENTIMETER:
        fprintf(fout, "RESUNIT_CENTIMETER\n");
        break;
    default:
        fprintf(fout, "UNKNOWN\n");
        break;
    }

    fprintf(fout, "RESOLUTION: %f %f (pixels per <unit>)\n",
            T->xresolution, T->yresolution);

    if(T->imagedescription != NULL)
    {
        fprintf(fout, "Image description: '%s'\n", T->imagedescription);
    }

    if(T->software != NULL)
    {
        fprintf(fout, "Software: '%s'\n", T->software);
    }

    return;
}

float * fim_tiff_read_sub(const char * fName,
                          ttags * T,
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

    if(T != NULL)
    {
        ttags_get(tfile, T);
    }

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
            fprintf(fim_tiff_log, "TIFFTAG_COMPRESSION=%u is not supported\n", CMP);
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
        fprintf(fim_tiff_log, "Warning: TIFFTAG_SAMPLEFORMAT not specified, "
                "assuming uint but that could be wrong!\n");
        isUint = 1;
    }

    if(!(isUint || isFloat))
    {
        fprintf(fim_tiff_log, "fim_tiff: Only unsigned integer and float images are supported\n");
        return NULL;
    }

    if(isUint && !(BPS == 16 || BPS == 8))
    {
        fprintf(fim_tiff_log, "fim_tiff: Unsigned %d-bit images are not supported, only 8 and 16.\n", BPS);
        return NULL;
    }

    if(isFloat && (BPS != 32))
    {
        fprintf(fim_tiff_log, "ERROR: For floating point images, only 32-bit "
                "samples are supported (this image has %u bits per sample)\n",
                BPS);
        exit(EXIT_FAILURE);
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
        fprintf(fim_tiff_log, "fim_tiff: WARNING: Could not read "
                "TIFFTAG_PHOTOMETRIC, assuming minIsBlack\n");
        PTAG=1;
    }
    //printf("PTAG = %u\n", PTAG);
    if(PTAG == 0)
    {
        inverted = 1;
    }

    if(PTAG > 1)
    {
        fprintf(fim_tiff_log, "fim_tiff WARNING: Only WhiteIsZero or BlackIsZero "
                "are supported Photometric Interpretations, tag found: %u\n",
                PTAG);
        PTAG = 1;
        fprintf(fim_tiff_log, "Assuming min-is-black\n");
        fflush(stdout);
    }

    tmsize_t ssize = TIFFStripSize(tfile); // Seems to be in bytes
    uint32_t nstrips = TIFFNumberOfStrips(tfile);
    uint32_t ndirs = TIFFNumberOfDirectories(tfile);
    P0[0] = (size_t) ndirs;
    int64_t P = P0[0];

    if(verbosity > 1)
    {
        if(gotB)
        {
            fprintf(fim_tiff_log, " TIFFTAG_IMAGEDEPTH: %u\n", B);
        }
        fprintf(fim_tiff_log, " TIFFTAG_BITSPERSAMPLE: %u\n", BPS);
        fprintf(fim_tiff_log, " size: %zu x %zu, %zu bits\n", (size_t) M, (size_t) N, (size_t) BPS);
        fprintf(fim_tiff_log, " # strips: %zu \n", (size_t) nstrips);
        fprintf(fim_tiff_log, " strip size (ssize): %zu bytes \n", (size_t) ssize);
        fprintf(fim_tiff_log, " # dirs (slices): %zu\n", (size_t) ndirs);
    }

    if(TIFFIsTiled(tfile))
    {
        fprintf(stderr, "fim_tiff: Tiled files are not supported\n");
        TIFFClose(tfile);
        return NULL;
    }

    // assert(M*N*BPS/8 == nstrips * ssize);

    //  size_t nel = nstrips * ssize * ndirs;
    float * V = NULL;
    size_t nel = M*N*ndirs;
    if(subregion)
    {
        nel = wM*wN*wP;
    }

    V = fftwf_malloc(nel*sizeof(float));
    assert(V != NULL);
    memset(V, 0, nel*sizeof(float));

    if(isFloat)
    {
        if(verbosity > 1)
        {
            fprintf(fim_tiff_log, "ReadFloat ...\n");
        }
        if(subregion)
        {
            perror("readFloat_sub is not implemented.\n");
            exit(EXIT_FAILURE);
            //readFloat_sub(tfile, V, ssize, ndirs, nstrips, M*N,
            //sM, sN, sP, wM, wN, wP);
        } else {
            readFloat(tfile, V, ssize, ndirs, nstrips, M*N);
        }
    }
    if(isUint)
    {
        if(verbosity > 1)
        {
            fprintf(fim_tiff_log, "ReadUint ...\n");
        }
        if(BPS == 16)
        {
            if(subregion)
            {
                fprintf(fim_tiff_log, "%" PRId64 " %" PRId64 " %" PRId64 ", %"
                        PRId64 " %" PRId64 " %" PRId64 ", %" PRId64 " %"
                        PRId64 " %" PRId64 "\n",
                        M, N, P, sM, sN, sP, wM, wN, wP);
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
                fprintf(fim_tiff_log, "%" PRId64 " %" PRId64 " %" PRId64 ", %"
                        PRId64 " %" PRId64 " %" PRId64 ", %" PRId64 " %"
                        PRId64 " %" PRId64 "\n",
                        M, N, P, sM, sN, sP, wM, wN, wP);
                readUint8_sub(tfile, V, ssize, ndirs, nstrips, M*N,
                              M, N, (int) P, sM, sN, sP, wM, wN, wP);
            } else {
                readUint8(tfile, V, ssize, ndirs, nstrips, M*N);
            }
        }
    }

    TIFFClose(tfile);
    if(verbosity>1)
    {
        fprintf(fim_tiff_log, "Done reading\n");
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

    ttags * T = malloc(sizeof(ttags));
    float * I = (float *) fim_tiff_read(inname, T, &M, &N, &P, 1);

    ttags_show(stdout, T);

    if(I == NULL)
    {
        printf("Failed to read %s\n", inname);
        exit(1);
    }

    floatimage_show_stats(I, M, N, P);
    // floatimage_normalize(I, M*N*P);

    ttags_show(stdout, T);

    fim_tiff_write(outname, I, T, M, N, P);
    ttags_free(&T);

    free(I);
    free(outname);

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
    int verbose = 0;

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

    ttags * T = malloc(sizeof(ttags));;
    ttags_get(input, T);
    ttags_set_software(T, "deconwolf " deconwolf_version);

    uint32_t SF = SAMPLEFORMAT_UINT;
    int gotSF = TIFFGetField(input, TIFFTAG_SAMPLEFORMAT, &SF);
    if(gotSF != 1)
    {
        printf("Warning: Unable to determine the sample format "
               "of %s\n, assuming uint but that could be wrong.\n", in);
    }

    uint32_t BPS;
    int gotBPS = TIFFGetField(input, TIFFTAG_BITSPERSAMPLE, &BPS);
    if(gotBPS != 1)
    {
        fprintf(stderr, "Unable to get the number of bits per sample\n");

    }

    TIFF * output = TIFFOpen(out, "w");
    ttags_set(output, T);
    ttags_free(&T);

    TIFFSetField(output, TIFFTAG_IMAGEWIDTH, M);
    TIFFSetField(output, TIFFTAG_IMAGELENGTH, N);
    TIFFSetField(output, TIFFTAG_SAMPLESPERPIXEL, 1);
    TIFFSetField(output, TIFFTAG_ROWSPERSTRIP, M);
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
    // Irrelevant if samples per pixel is set to 1
    TIFFSetField(output, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField(output, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK);
    TIFFSetField(output, TIFFTAG_SUBFILETYPE, FILETYPE_PAGE);
    TIFFSetField(output, TIFFTAG_PAGENUMBER, 0, 1);

    tmsize_t ssize = TIFFStripSize(input); // Seems to be in bytes
    //  printf("Strip size: %zu b\n", (size_t) ssize);

    uint32_t nstrips = TIFFNumberOfStrips(input);

    if(verbose > 1)
    {
        printf("strip size: %ld\n", ssize);
        printf("number of strips: %u\n", nstrips);
        switch(SF)
        {
        case SAMPLEFORMAT_UINT:
            printf("Sample format: UINT\n");
            break;
        case SAMPLEFORMAT_IEEEFP:
            printf("Sample format: float\n");
            break;
        default:
            printf("Unsupported sample format\n");
            break;
        }
        printf("Bits per sample: %u\n", BPS);
        printf("Number of directories (Z-planes): %ld\n", P);
    }

    TIFFSetDirectory(output, 0);

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
            tsize_t written = TIFFWriteRawStrip(output, nn, mstrip, read);
            assert(written == read);

        }

        _TIFFfree(strip);
        _TIFFfree(mstrip);
    }

    // Input is 32-bit float
    if(SF == SAMPLEFORMAT_IEEEFP)
    {
        //TIFFSetField(output, TIFFTAG_STRIPBYTECOUNTS, M*N*sizeof(float));
        //TIFFSetField(output, TIFFTAG_STRIPBYTECOUNTS, 1);

        float * mstrip = _TIFFmalloc(ssize); // For max over all directories
        float * strip    = _TIFFmalloc(ssize);
#if 0
        float * outbuffer = malloc(M*N*sizeof(float));
        size_t outbufferpos = 0;
#endif
        for(int64_t ss = 0; ss < nstrips; ss++) // Each strip
        {
            memset(mstrip, 0, ssize);
            tsize_t read = 0;
            for(int64_t pp = 0; pp<P; pp++) // For each directory
            {
                TIFFSetDirectory(input, pp); // Does it keep the location for each directory?
                read = TIFFReadEncodedStrip(input, ss, strip, (tsize_t)-1);
                for(int64_t kk = 0; kk<read/4; kk++)
                {
                    if(strip[kk] > mstrip[kk])
                    {
                        mstrip[kk] = strip[kk];
                    }
                }
            }


            tsize_t written = TIFFWriteRawStrip(output, ss, mstrip, read);
            assert(written == read);

        }
        _TIFFfree(strip);
        _TIFFfree(mstrip);

    }

    TIFFClose(input);

    TIFFWriteDirectory(output);
    TIFFClose(output);

    return 0;
}


int fim_tiff_extract_slice(char * in, char * out, int slice)
{

    if(slice < 1)
    {
        printf("Slice %d does not make sense\n", slice);
        return -1;
    }

    int64_t M, N, P;

    if(fim_tiff_get_size(in, &M, &N, &P))
    {
        printf("Can't open %s to get image dimension\n", in);
        return -1;
    }

    if(slice > P)
    {
        printf("Can't extract slice %d from an image with %" PRId64 " slices\n", slice, P);
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

        TIFFSetDirectory(input, slice-1); // Does it keep the location for each directory?

        for(int64_t nn = 0; nn<nstrips; nn++) // Each strip
        {
            memset(mstrip, 0, ssize);
            tsize_t read = 0;


            read = TIFFReadEncodedStrip(input, nn, strip, (tsize_t)-1);
            for(int64_t kk = 0; kk<read/2; kk++)
            {
                mstrip[kk] = strip[kk];

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

        TIFFSetDirectory(input, slice-1); // Does it keep the location for each directory?
        for(int64_t nn = 0; nn<nstrips; nn++) // Each strip
        {
            memset(mstrip, 0, ssize);
            tsize_t read = 0;


            read = TIFFReadEncodedStrip(input, nn, strip, (tsize_t)-1);
            for(int64_t kk = 0; kk<read/4; kk++)
            {
                mstrip[kk] = strip[kk];
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
