#define deconwolf_version "alpha-0.001"

typedef struct{
  int nThreads;
  int nIter;
  char * imFile;
  char * psfFile;
  char * outFile;
  char * logFile;
  FILE * log;

  int verbosity;
  fftwf_plan fft_plan;
  fftwf_plan ifft_plan;
} opts;

int deconwolf(opts *);

void usage(int argc, char ** argv);
void unittests();
void shift_vector_ut();
void fArray_flipall_ut();
