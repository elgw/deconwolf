#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <complex.h>
#include <stdlib.h>


// Inspiration:
//https://github.com/Biomedical-Imaging-Group/PSFGenerator/tree/master/src/psf/bornwolf

typedef struct {
  double lambda; // Emission maxima
  double k0;
  double NA; // Numerical aperture of lens
  double ni;
  double TOL;
  size_t K; // Number of iterations
  double resLateral;
  double resAxial;
  // shape of output image
  int M;
  int N;
  int P;
} bw_conf;


double complex integrand(double rho, double r, double defocus, bw_conf * conf) {

  // 'rho' is the integration parameter.
  // 'r' is the radial distance of the detector relative to the optical
  // axis.
  // NA is assumed to be less than 1.0, i.e. it assumed to be already
  // normalized by the refractive index of the immersion layer, ni.
  // The return value is a complex number.

  double k0 = 2 * M_PI / conf->lambda;
  double BesselValue = j0(conf->k0 * conf->NA * r * rho);

  double OPD; // Optical path difference
  double W; // Phase aberrations.

  OPD = pow(conf->NA,2) * defocus * pow(rho,2) / (2.0 * conf->ni);
  W = k0 * OPD;

  double complex y = BesselValue * cos(W) * rho 
    -I*BesselValue * sin(W) * rho;
  return y;
}



// calculate()
// Simpson approximation for the Kirchhoff diffraction integral
// 'r' is the radial distance of the detector relative to the optical axis.
double calculate(double r, double defocus, bw_conf * conf) {

  // IJ.log("(p.ti0, p.ti) = (" + p.ti0 + ", " + p.ti + ")");
  double del = 0.5; // integration interval
  double curDifference = conf->TOL; // Stopping criterion

  complex double valueX0 = 0, valueXn = 0;
  complex double value = 0;

  double curI = 0.0, prevI = 0.0;

  // Initialization of the Simpson sum (first iteration)
  int N = 2; // number of sub-intervals
  int k = 0; // number of consecutive successful approximations
  double rho = 0.5;
  printf("r = %f, rho = %f\n", r, rho);
  complex double sumOddIndex = integrand(rho, r, defocus, conf);
  complex double sumEvenIndex = 0;

  valueX0 = integrand(0.0, r, defocus, conf);
  valueXn = integrand(1.0, r, defocus, conf);

  double complex sum = valueX0 + 2*sumEvenIndex + 4*sumOddIndex + valueXn;
  curI = pow(creal(sum),2) + pow(cimag(sum), 2) * del * del;

  prevI = curI;
  
  // Finer sampling grid until we meet the TOL value with the specified
  // number of repetitions, K
  size_t iteration = 0;
  while (k < conf->K && iteration < 10000) {
    iteration++;
    N *= 2;
    del = del / 2;

    sumEvenIndex = sumEvenIndex + sumOddIndex;

    sumOddIndex = 0;
    for (int n = 1; n < N; n = n + 2) {
      rho = n * del;
      value = integrand(rho, r, defocus, conf);
      sumOddIndex += value;
    }

    complex double sum = valueX0 + 2*sumEvenIndex + 4* sumOddIndex + valueXn;
    curI = (pow(creal(sum),2) + pow(sum, 2)) * del * del;

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
    prevI = curI;
  }

  return curI;
}


void BornWolf(double * V, float z, bw_conf * conf)
{

  // The center of the image in units of [pixels]
  double x0 = (conf->M - 1) / 2.0;
  double y0 = (conf->N - 1) / 2.0;

  // Radial locations.
  // double xpAbs = Math.abs(xp), ypAbs = Math.abs(yp);
  // double maxRadialDistanceInPixels =
  // Math.round(Math.sqrt((xpAbs+nx-x0)*(xpAbs+nx-x0)+(ypAbs+ny-y0)*(ypAbs+ny-y0)))+1;
  int maxRadius = (int) round(sqrt(pow(conf->M - x0,2) + pow(conf->N - y0,2))) + 1;
  int OVER_SAMPLING = 1;
  size_t nr = maxRadius*OVER_SAMPLING;
  double * r = malloc(nr * sizeof(double));
  double * h = malloc(nr * sizeof(double));

  for (int n = 0; n < nr; n++) {
    r[n] = ((double) n) / ((double) OVER_SAMPLING);
    h[n] = calculate(r[n] * conf->resLateral * 1E-9, z, conf);
  }

  // Linear interpolation of the pixels values
  int midx = (conf->M+1)/2 - 1;
  int midy = (conf->N+1)/2 -1;

  for (int x = 0; x < conf->M; x++) {
    for (int y = 0; y < conf->N; y++) {
      // radius of the current pixel in units of [pixels]
      double rPixel = sqrt(pow(x-midx, 2) + pow(y-midy, 2));
      // Index of nearest coordinate from bellow
      int index = (int) floor(rPixel * OVER_SAMPLING);
      assert(index < nr);
      // Interpolated value.
      V[x + conf->M * y] = h[index] + (h[index + 1] - h[index]) * (rPixel - r[index]) * OVER_SAMPLING;
    }
  }
  return;
}

  int main(int argc, char ** argv)
  {
    bw_conf conf;
    conf.lambda = 600;
    conf.k0 = 1;
    conf.NA = 1.4;
    conf.ni = 1.5;
    conf.TOL = 1e-1;
    conf.K = 9; // corresponding to "best" in PSFGenerator
    conf.M = 11;
    conf.N = 11;
    conf.P = 11;
    conf.resAxial = 300;
    conf.resLateral = 130;

    // Just testing ...
    double complex cy = integrand(200, 200, 0, &conf);
    printf("Integrand = %.2f + i%.2f\n", creal(cy), cimag(cy));
    double y = calculate(0, 0, &conf);
    printf("calculate(0) = %f\n", y);

    double * V = malloc(conf.M*conf.N*conf.P*sizeof(double));
    for (int z = 0; z < conf.P; z++) {
      double defocus = conf.resAxial * 1E-9 * (z - (conf.P - 1.0) / 2.0);
      printf("defocus: %f\n", defocus);
      BornWolf(V + z*conf.M*conf.N, defocus, &conf);
    }

  return 0;
}
