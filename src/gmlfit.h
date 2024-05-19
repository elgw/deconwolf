#ifdef NDEBUG
#define GSL_RANGE_CHECK_OFF
#endif

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <omp.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_statistics.h>


/** Fit a 3D gaussian to dots in images.
 * - Fitting function - Log Likelihood for Poissonian noise.
 * - Optimization Algorithm: Broyden-Fletcher-Goldfarb-Shanno algoritm (bfgs2 in GSL)
 * - Parallelization with OpenMP.
 *
 * @param V: The volume where the dots resides
 * @param M, N, P: The size of V
 * @param sigma_xy lateral size of the dots, this will be used as the starting
 *        for optimization and will also determine how large regions that should
 *        be investigated for each dot.
 * @param sigma_z the axial size of the dots.
 * @param X: 3xnX list of dot coordinates to start from
 * @param onbased: set to 1 if the coordinates in X are 1-based (i.e. for
 *        matlab interop).
 * @return Returns a table of size [nX, 9] with columns:
 *         0: background level
 *         1: number of photons in signal
 *         2: Peak signal intensity (according to the model)
 *         4: x
 *         5: y
 *         6: z
 *         7: sigma_xy (same sigma for x and y / isotropic in the lateral plane)
 *         8: sigma_z
 *         9: Status from GSL (27 = Found a local minima)
 *         9: final error
 *
 * You might want to discard dots that moved more than sqrt(3) pixels.
 *
 * Updated: 2024-05-06.
 **/

typedef float ifloat;

typedef struct{
    const float * image;
    size_t M; /* First dimension, stride 1 */
    size_t N; /* 2nd dimension, stride M */
    size_t P; /* 3rd dimension, stride M*N */
    const double * X; /* 3 x nX coordinates of start points */
/* Optional per dot scaling. sigma is multiplied with this if available.*/
    const double * DS;

    size_t nX;
    int x_base; /* 0-based or 1-based coordinates in X */
    double sigma_xy; /* Start size in lateral plane */
    double sigma_z; /* Start size in axial plane */
    FILE * log; /* File to write to */
    int verbose; /* 0=quiet, 1=some info, 2=debug */
    size_t max_iter; /* Number of iterations to run */
    int num_threads; /* to set the number of omp threads */
} gmlfit;

/* Returns a struct with default settings. After you are done, just free it */
gmlfit * gmlfit_new();

/* Run the problem. Returns NULL if the settings are invalid. */
double * gmlfit_run(gmlfit *);

/* Minimal set of test routines, mainly to be able to run the code
 * isolated and test with valgrind. Also provides example usage.
 */
int
gmlfit_ut(int argc, char ** argv);
