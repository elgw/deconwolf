#include <fftw3.h>

void myfftw_start(void);
void myfftw_stop(void);

void dim3_real_float_inverse(fftwf_complex * in, float * out,
    const int n1, const int n2, const int n3);

void dim3_real_float(float * in, fftwf_complex* out,
    const int n1, const int n2, const int n3);


