#include <fftw3.h>

void myfftw_start(void);
void myfftw_stop(void);

void dim3_real_float_inverse(fftwf_complex * in, float * out,
    const int n1, const int n2, const int n3);

void dim3_real_float(float * in, fftwf_complex* out,
    const int n1, const int n2, const int n3);

fftwf_complex * fft(float * in, int n1, int n2, int n3);
void fft_mul(fftwf_complex * C, fftwf_complex * A, fftwf_complex * B, size_t N);

float * fft_convolve_cc(fftwf_complex * A, fftwf_complex * B, int M, int N, int P);


void fft_train(size_t, size_t, size_t);
