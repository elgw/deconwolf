#include "lanczos3.h"

static double sinc(double x0)
{
    double x = M_PI*x0;
    if(fabs(x) > 1e-6)
    {
        return sin(x)/(x);
    } else {
        return 1 - pow(x,2)/6 + pow(x,4)/120;
    }

}

static double lanczos3_weight(double x)
{
    return sinc(x)*sinc(x/3.0);
}

double lanczos3(const double * v, size_t nV, double x)
{
    /*
     * Lanczos-3 interpolation of v at floating point indices
     * assumes symmetry around 0
     * assumes that x >= 0 and that x+3 < numel(v)
     */

    assert(x>= 0);

    /* Integer index before x */
    int n = (int) floor(x);
    /* Can handled by some other boundary condition or
     * interpolation method */
    assert( (size_t) n+2 < nV);

    if((size_t) n+3 >= nV)
    {
        printf("Can't interpolate x=%f in a %zu vector", x, nV);
        exit(EXIT_FAILURE);
    }

    /* Distance to integer index before x */
    double d = x - (double) n;

    /* If very close to integer point, don't interpolate */
    if(d < 1e-8)
    {
        return v[n];
    }

    if( d > 1-1e-8)
    {
        return v[n+1];
    }

    double y = 0;
    double sumw = 0;
    /* three points before */
    for(int kk = 0; kk<3; kk++)
    {
        int idx = abs(n-kk);
        double w = lanczos3_weight(d + (double) kk);
        y += v[idx]*w;
        sumw+=w;
    }

    /* three points after */
    for(int kk = 1; kk<4; kk++)
    {
        double w = lanczos3_weight((double) kk - d);
        y += v[n+kk]*w;
        sumw+=w;
    }

    /* Ad-hoc restoration of the partition of unity property
     * Numerical experiments confirms that this is useful in this application
     */
    y /= sumw;

    /* Lanczos-3 can produce negative output when all input points
     * are positive. We don't want that.*/

    if(y > 0)
    {
        return y;
    } else {
        return 0;
    }
}
