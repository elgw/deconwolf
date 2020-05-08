#define deconwolf_version "alpha-0.002"

/* deconwolf
 * Erik Wernersson, 2020
 *
 * Configuration files.
 * If $home/.config/ exists deconwolf will store its settings there
 * specifically the fftw3 wisdom data is stored in
 * ./config/deconwolf/wisdom/
 *
 * The choice of folder is recommended by "XDG Base Directory Specification", see:
 * https://specifications.freedesktop.org/basedir-spec/basedir-spec-latest.html
 */

typedef struct{
  int nThreads;
  int nIter;
  char * imFile;
  char * psfFile;
  char * outFile;
  char * logFile;
  char * prefix;
  FILE * log;
  int tiling_maxSize;
  int tiling_padding;
  int overwrite; // overwrite output tif file?

  int verbosity;
  fftwf_plan fft_plan;
  fftwf_plan ifft_plan;
} opts;

int deconwolf(opts *);

void usage(const int argc, char ** argv, const opts * );
void unittests();
void shift_vector_ut();
void fArray_flipall_ut();
float * autocrop_psf(float *, int *, int *, int *, int, int , int, opts * );
