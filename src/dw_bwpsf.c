#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>
#include "fim_tiff.h"

// Inspiration:
//https://github.com/Biomedical-Imaging-Group/PSFGenerator/tree/master/src/psf/bornwolf
// gcc -Wall dw_bwpsf.c fim.c fim_tiff.c -lm -ltiff -lfftw3f

typedef struct {
  float lambda; // Emission maxima
  float NA; // Numerical aperture of lens
  float ni;
  float TOL;
  size_t K; // Number of iterations
  float resLateral;
  float resAxial;
  // shape of output image
  int M;
  int N;
  int P;
} bw_conf;

void bw_conf_printf(bw_conf * conf, FILE * out)
{
  fprintf(out, "lambda = %.2f nm\n", conf->lambda*1e9);
  fprintf(out, "NA = %f\n", conf->NA);
  fprintf(out, "ni = %f\n", conf->ni);
  fprintf(out, "TOL = %f\n", conf->TOL);
  fprintf(out, "resLateral = %.2f nm\n", conf->resLateral*1e9);
  fprintf(out, "resAxial = %.2f nm\n", conf->resAxial*1e9);
  fprintf(out, "out size: [%d x %d x %d] pixels\n", conf->M, conf->N, conf->P);
}

float complex integrand(float rho, float r, float defocus, bw_conf * conf) {

  // 'rho' is the integration parameter.
  // 'r' is the radial distance of the detector relative to the optical
  // axis.
  // NA is assumed to be less than 1.0, i.e. it assumed to be already
  // normalized by the refractive index of the immersion layer, ni.
  // The return value is a complex number.

  assert(rho<=1); assert(rho>=0);

  float k0 = 2.0 * M_PI / conf->lambda;
  float BesselValue = j0(k0 * conf->NA * r * rho);

// Optical path difference
  float OPD = pow(conf->NA,2) * defocus * pow(rho,2) / (2.0 * conf->ni);
// Phase aberrations.
  float W = k0 * OPD;

  float complex y = BesselValue*cos(W)*rho - I*BesselValue*sin(W)*rho;
  return y;
}



// Simpson approximation for the Kirchhoff diffraction integral
// 'r' is the radial distance of the detector relative to the optical axis.
float calculate(float r, float defocus, bw_conf * conf) {
// is the stopping criterion really good?
// why del^2 -- isn't it a 1D integral?
// Better set a tolerance relative to the wanted precision (never more than 32-bit)
// doubling the number of points 9 times!!

  float del = 0.5; // integration interval
  float curDifference = conf->TOL; // Stopping criterion

  complex float value = 0;

  float curI = 0.0, prevI = 0.0;

  // Initialization of the Simpson sum (first iteration)
  int N = 2; // number of sub-intervals
  int k = 0; // number of consecutive successful approximations
  float rho = 0.5;

  complex float sumOddIndex = integrand(rho, r, defocus, conf);
  complex float sumEvenIndex = 0;

  complex float valueX0 = integrand(0.0, r, defocus, conf);
  complex float valueXn = integrand(1.0, r, defocus, conf);

  float complex sum = valueX0 + 2*sumEvenIndex + 4*sumOddIndex + valueXn;
  curI = (pow(creal(sum),2) + pow(cimag(sum), 2)) * pow(del,2);

  prevI = curI;

  // Finer sampling grid until we meet the TOL value with the specified
  // number of repetitions, K
  size_t iteration = 1;
  while (k < conf->K && iteration < 10000) {
    iteration++;
    N *= 2;
    del /= 2;

    sumEvenIndex = sumEvenIndex + sumOddIndex;
    sumOddIndex = 0 + I*0;

    for (int n = 1; n < N; n = n + 2) {
      rho = n * del;
      value = integrand(rho, r, defocus, conf);
      sumOddIndex += value;
    }

    complex float sum = valueX0 + 2*sumEvenIndex + 4*sumOddIndex + valueXn;
    curI = (pow(creal(sum),2) + pow(cimag(sum), 2)) * pow(del, 2);

    // Relative error between consecutive approximations
    if (prevI == 0.0)
    {
      curDifference = fabs((prevI - curI) / 1E-5);
    } else {
      curDifference = fabs((prevI - curI) / curI);
    }
    if (curDifference <= conf->TOL)
    {
      k++;
    } else {
      k = 0;
    }
//    if(fabs(curI-prevI) < 1e-12)
//    {
//      break;
//    }
    prevI = curI;
//    printf("Iteration: %zu, curDifference: %e, curI: %e\n", iteration, curDifference, curI);
  }
//  printf("\n");
//  printf("@%d: curI = %f\n", k, curI);
  return curI;
}


void BornWolf(float * V, float z, bw_conf * conf)
{

  // The center of the image in units of [pixels]
  float x0 = (conf->M - 1) / 2.0;
  float y0 = (conf->N - 1) / 2.0;

  // Radial locations.
  // float xpAbs = Math.abs(xp), ypAbs = Math.abs(yp);
  // float maxRadialDistanceInPixels =
  // Math.round(Math.sqrt((xpAbs+nx-x0)*(xpAbs+nx-x0)+(ypAbs+ny-y0)*(ypAbs+ny-y0)))+1;
  int maxRadius = (int) round(sqrt(pow(conf->M - x0, 2) + pow(conf->N - y0, 2))) + 1;
  int OVER_SAMPLING = 1;
  size_t nr = maxRadius*OVER_SAMPLING;
  float * r = malloc(nr * sizeof(float));
  float * h = malloc(nr * sizeof(float));

  for (int n = 0; n < nr; n++) {
    r[n] = ((float) n) / ((float) OVER_SAMPLING);
    h[n] = calculate(r[n] * conf->resLateral, z, conf);
//    printf("r[%d=%e] = %e, h[%d] = %e\n", n, r[n]*conf->resLateral, r[n], n, h[n]);
  }

  assert(r[0] == 0);
  //exit(1);
  for (int x = 0; 2*x <= conf->M; x++) {
    for (int y = 0; 2*y <= conf->N; y++) {
      // radius of the current pixel in units of [pixels]
      float rPixel = sqrt(pow(x-x0, 2) + pow(y-y0, 2));
      // Index of nearest coordinate from below (replace by pixel integration)
      int index = (int) floor(rPixel * OVER_SAMPLING);
      assert(index < nr);
      // Interpolated value.
      float v = h[index] + (h[index + 1] - h[index]) * (rPixel - r[index]) * OVER_SAMPLING;
      
      int xf = conf->M-x-1;
      int yf = conf->N-y-1;

      V[x + conf->M * y] = v;
      V[y + conf->M * x] = v;
      V[xf + conf->M * y] = v;
//      V[yf + conf->M * x] = v;
      V[x + conf->M * yf] = v;
//      V[y + conf->M * xf] = v;
      V[xf + conf->M * yf] = v;
      V[yf + conf->M * xf] = v;


    }
  }
  free(r);
  free(h);
  return;
}

int main(int argc, char ** argv)
{
  bw_conf conf;
  conf.lambda = 600*1e-9;
  conf.NA = 1.4;
  conf.ni = 1.5;
  conf.TOL = 1e-1;
  conf.K = 9; // corresponding to "best" in PSFGenerator
  conf.M = 181;
  conf.N = 181;
  conf.P = 181;
  conf.resAxial = 300*1e-9;
  conf.resLateral = 130*1e-9;

  bw_conf_printf(&conf, stdout);

//  conf.NA/=conf.ni;
//  conf.ni = 1;

  // Just testing ...
  if(0){
    float complex cy = integrand(200, 200, 0, &conf);
    printf("Integrand = %.2f + i%.2f\n", creal(cy), cimag(cy));
    float y = calculate(0, 0, &conf);
    printf("calculate(0) = %f\n", y);
  }

  float * V = malloc(conf.M*conf.N*conf.P*sizeof(float));
  for (int z = 0; z <= (conf.P-1)/2; z++) {
    float defocus = conf.resAxial * (z - (conf.P - 1.0) / 2.0);
    BornWolf(V + z*conf.M*conf.N, defocus, &conf);
  }

  // symmetry in Z
  size_t MN = conf.M*conf.N;
  for(int z = 0; z<(conf.P-1)/2; z++)
  {
    memcpy(V+(conf.P-z-1)*MN, V+z*MN, MN*sizeof(float));
  }

  printf("Writing to psf.tif\n");
  fim_tiff_write("psf.tif", V, conf.M, conf.N, conf.P);
  free(V);

  return 0;
}
