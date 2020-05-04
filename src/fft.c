#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>

const char swf[] = "fftw_wisdom_float_threads.dat";

void myfftw_start(void)
{
  int nThreads = 8;
  fftwf_init_threads();
  fftwf_plan_with_nthreads(nThreads);
  fftwf_import_wisdom_from_filename(swf);
}

void myfftw_stop(void)
{
  fftwf_export_wisdom_to_filename(swf);
  fftwf_cleanup_threads();
  fftwf_cleanup();
}

void dim3_real_float_inverse(fftwf_complex * in, float * out,
    const int n1, const int n2, const int n3)
{
  myfftw_start();

  fftwf_plan p;

  p = fftwf_plan_dft_c2r_3d(n3, n2, n1, 
      in, out, 
      FFTW_PATIENT);
  fftwf_execute(p);
  fftwf_destroy_plan(p);

  myfftw_stop();
}

void dim3_real_float(float * in, fftwf_complex* out,
    const int n1, const int n2, const int n3)
{

  if(0){
  FILE * fh = fopen("/tmp/log.txt", "w");
  fprintf(fh, "fftwf_version=%s\n", fftwf_version);
  fclose(fh);
  }

  myfftw_start();
  fftwf_plan p = fftwf_plan_dft_r2c_3d(n3, n2, n1, 
      in, // Float
      out, // fftwf_complex 
      FFTW_PATIENT);

  fftwf_execute(p); 
  fftwf_destroy_plan(p);
  myfftw_stop();
}

