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
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>
#include <string.h>
#include <dirent.h>
#include <errno.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "fim.h"

// FFTW_MEASURE, FFTW_PATIENT or FFTW_EXHAUSTIVE
const unsigned int MYPLAN = FFTW_MEASURE;
//const unsigned int MYPLAN = FFTW_PATIENT;

static int isdir(char * dir)
{
  /* Check if directory exist, do not create if missing 
   * returns 1 if it exist
   * */

  DIR* odir = opendir(dir);
  if (odir) {
    /* Directory exists. */
    closedir(odir);
    return 1;
  } else if (ENOENT == errno) {
    /* Directory does not exist. */
    return 0;
  } else {
    /* opendir() failed for some other reason. */
    return 0;
  }
}

static int ensuredir(char * dir)
  /* Create dir if it does not exist. 
   * Returns 0 if the dir already existed or could be created
   * returns non-zeros if the dir can't be created
   */
{
  if(isdir(dir) == 1)
  {
    return 0;
  }

  if(mkdir(dir, 0700) == 0)
  {
    return 0;
  }

  return 1;
}

static char * get_swf_file_name(int nThreads)
{
  char * dir_home = getenv("HOME");
  char * dir_config = malloc(1024*sizeof(char));
  char * swf = malloc(1024*sizeof(char));
  sprintf(swf, "fftw_wisdom_float_threads_%d.dat", nThreads);    

  sprintf(dir_config, "%s/.config/", dir_home);
  // printf("dir_config = %s\n", dir_config);
  if( !isdir(dir_config) )
  {     
    free(dir_config);
    return swf;
  }

  sprintf(dir_config, "%s/.config/deconwolf/", dir_home);
  //  printf("dir_config = %s\n", dir_config);
  if( ensuredir(dir_config) == 0 )
  {
    char * prefered = malloc(1024*sizeof(char));
    sprintf(prefered, "%s%s", dir_config, swf);
    free(dir_config);
    free(swf);
    return prefered;
  } else {
    free(dir_config);
    return swf;
  }
}

void myfftw_start(const int nThreads)
{  
  //  printf("\t using %s with %d threads\n", fftwf_version, nThreads);
  fftwf_init_threads();
  fftwf_plan_with_nthreads(nThreads);
  char * swf = get_swf_file_name(nThreads);
  if(swf == NULL)
  { 
    assert(0);
  } 
  else
  {
    fftwf_import_wisdom_from_filename(swf);
    free(swf);
  }
}

void myfftw_stop(void)
{
  fftwf_cleanup_threads();
  fftwf_cleanup();
  // Note: wisdom is only exported by fft_train
}


fftwf_complex * fft(float * restrict in, const int n1, const int n2, const int n3)
{
  int n3red = (n3+3)/2;
  fftwf_complex * out = fftwf_malloc(n1*n2*n3red*sizeof(fftwf_complex));
  memset(out, 0, n1*n2*n3red*sizeof(fftwf_complex));

  fftwf_plan p = fftwf_plan_dft_r2c_3d(n3, n2, n1, 
      in, // Float
      out, // fftwf_complex 
      MYPLAN);
  fftwf_execute(p); 
  fftwf_destroy_plan(p);
  return out;
}

void fft_mul(fftwf_complex * restrict C, 
    fftwf_complex * restrict A, 
    fftwf_complex * restrict B, 
    const size_t n1, const size_t n2, const size_t n3)
{
  const int n3red = (n3+3)/2;
  const size_t N = n1*n2*n3red;
  // C = A*B
  for(size_t kk = 0; kk<N; kk++)
  {
    float a = A[kk][0]; float ac = A[kk][1];
    float b = B[kk][0]; float bc = B[kk][1];
    C[kk][0] = a*b - ac*bc;
    C[kk][1] = a*bc + b*ac;
  }
  return;
}


void fft_mul_conj(fftwf_complex * restrict C, 
    fftwf_complex * restrict A, 
    fftwf_complex * restrict B, 
    const size_t n1, const size_t n2, const size_t n3)
  /* Multiply and conjugate the elements in the array A 
   * i.e. C = conj(A)*B
   * All inputs should have the same size [n1 x n2 x n3]
   * */
{
  const int n3red = (n3+3)/2;
  const size_t N = n1*n2*n3red;
  // C = A*B
  for(size_t kk = 0; kk<N; kk++)
  {
    float a = A[kk][0]; float ac = -A[kk][1];
    float b = B[kk][0]; float bc = B[kk][1];
    C[kk][0] = a*b - ac*bc;
    C[kk][1] = a*bc + b*ac;
  }
  return;
}

float * fft_convolve_cc(fftwf_complex * A, fftwf_complex * B, 
    const int M, const int N, const int P)
{
  size_t n3red = (P+3)/2;
  fftwf_complex * C = fftwf_malloc(M*N*n3red*sizeof(fftwf_complex));
  fft_mul(C, A, B, M, N, P); 

  float * out = fftwf_malloc(M*N*P*sizeof(float));

  fftwf_plan p = fftwf_plan_dft_c2r_3d(P, N, M, 
      C, out, 
      MYPLAN);
  fftwf_execute(p);
  fftwf_destroy_plan(p);
  fftwf_free(C);

  for(size_t kk = 0; kk<M*N*P; kk++)
  {
    out[kk]/=(M*N*P);
  }
  return out;
}

float * fft_convolve_cc_conj(fftwf_complex * A, fftwf_complex * B, 
    const int M, const int N, const int P)
{
  size_t n3red = (P+3)/2;
  fftwf_complex * C = fftwf_malloc(M*N*n3red*sizeof(fftwf_complex));
  fft_mul_conj(C, A, B, M, N, P); 

  float * out = fftwf_malloc(M*N*P*sizeof(float));

  fftwf_plan p = fftwf_plan_dft_c2r_3d(P, N, M, 
      C, out, 
      MYPLAN);
  fftwf_execute(p);
  fftwf_destroy_plan(p);
  fftwf_free(C);

  for(size_t kk = 0; kk<M*N*P; kk++)
  {
    out[kk]/=(M*N*P);
  }
  return out;
}

void fft_train(const size_t M, const size_t N, const size_t P, const int verbosity, const int nThreads)
{
  if(verbosity > 1){
    printf("fftw3 training ... \n"); fflush(stdout);
  }
  size_t MNP = M*N*P;
  fftwf_complex * C = fftwf_malloc(MNP*sizeof(fftwf_complex));
  float * R = fftwf_malloc(MNP*sizeof(float));

  fftwf_plan p0 = fftwf_plan_dft_c2r_3d(P, N, M, 
      C, R, MYPLAN | FFTW_WISDOM_ONLY);

  if(p0 == NULL)
  {    
    if(verbosity > 0)
    {
      printf("> generating c2r plan\n");
    }
    fftwf_plan p1 = fftwf_plan_dft_c2r_3d(P, N, M, 
        C, R, MYPLAN);
    fftwf_execute(p1);
    fftwf_destroy_plan(p1);
  } else {
    if(verbosity > 1)
    {
      printf("\tc2r -- ok\n");
    }
  }
  fftwf_destroy_plan(p0);

  p0 = fftwf_plan_dft_r2c_3d(P, N, M, 
      R, C, MYPLAN | FFTW_WISDOM_ONLY);

  if(p0 == NULL)
  {
    if(verbosity > 0){
      printf("> generating r2c plan \n");
    }
    fftwf_plan p2 = fftwf_plan_dft_r2c_3d(P, N, M, 
        R, C, MYPLAN);
    fftwf_execute(p2);
    fftwf_destroy_plan(p2);
  } else {
    if(verbosity > 1)
    {
      printf("\tr2c -- ok\n");
    }
  }
  fftwf_destroy_plan(p0);

  fftwf_free(R);
  fftwf_free(C);

  char * swf = get_swf_file_name(nThreads);
  if(swf == NULL)
  { assert(0); } 
  else {
    fftwf_export_wisdom_to_filename(swf);
  }

  return;
}

void fft_ut_wisdom_name(void){
  /* Wisdom file names
   * Try this when $HOME/.config/deconwolf/ does not exist
   * and when it does ... 
   * could also test it when that dir isn't writeable.
   * */

  int nThreads = 4;
  char * swf = get_swf_file_name(nThreads);
  printf("swf = '%s'\n", swf);
  free(swf);
}

void fft_ut_flipall_conj()
{
  myfftw_start(2);
  /* Test the identity 
   * flip(X) = ifft(conj(fft(X)))
   */
  int M = 12, N = 13, P = 15;
  float * A = fftwf_malloc(M*N*P*sizeof(float));
  for(int kk = 0; kk<M*N*P; kk++)
  { A[kk] = (float) rand() / (float) RAND_MAX; }
  fim_stats(A, M*N*P);
  float * B = fftwf_malloc(M*N*P*sizeof(float));
  memcpy(B, A, M*N*P*sizeof(float));
  float * B_flipall = fftwf_malloc(M*N*P*sizeof(float));
  fim_flipall(B_flipall, B, M, N, P);


  fftwf_complex * FA = fft(A, M, N, P);
  fftwf_complex * FB = fft(B, M, N, P);
  fftwf_complex * FB_flipall = fft(B_flipall, M, N, P);

  float * Y1 = fft_convolve_cc(FA, FB_flipall, M, N, P);
  float * Y2 = fft_convolve_cc_conj(FA, FB, M, N, P);

  float mse = fim_mse(Y1, Y2, M*N*P);
  printf("mse=%f ", mse);
  if(mse < 1e-5)
  { printf("ok!\n"); } else 
  { printf("BAD :(\n"); }

  fftwf_free(A); fftwf_free(FA);
  fftwf_free(B); fftwf_free(FB);
  fftwf_free(B_flipall); fftwf_free(FB_flipall);

  fftwf_free(Y1);
  fftwf_free(Y2);

  myfftw_stop();
}

void fft_ut(void)
{
  fft_ut_wisdom_name();
  fft_ut_flipall_conj();
  return;
}
