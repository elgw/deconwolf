#include "entropy.h"
#include <math.h>
#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

/* Return log(n!) using Stirling's formula with the peculiar behavior
 * that log_p(0) returns 0. For small values the function is exact.
 *
 * For Stirling's formula, see
 * https://doi.org/10.21468/SciPostPhysLectNotes.76
 *
 * Note:
 *
 * A direct lookup is used for small numbers since
 *
 * - the approximation is not that good in that range
 *
 * - benchmarks says that it is faster which can be understood since
 *    the function is called more often for small numbers which
 *    weights up the extra costs of branching.
 *
 * It is not tested how large the lookup table should be for optimal
 *   performance.
 */
static double log_p(const size_t n)
{

    switch(n)
    {
    case 0:
        return 0;
    case 1:
        return 0;
    case 2:
        return log(2);
    case 3:
        return log(2) + log(3);
    case 4:
        return log(2)+log(3)+log(4);
    case 5:
        return log(2)+log(3)+log(4)+log(5);
    case 6:
        return log(2)+log(3)+log(4)+log(5)+log(6);
    case 7:
        return log(2)+log(3)+log(4)+log(5)+log(6)+log(7);
    case 8:
        return log(2)+log(3)+log(4)+log(5)+log(6)+log(7)+log(8);
    default:
        // This takes around 40 cpu cycles on AMD 3700X
        return n*log(n) - n + log(sqrt(2.0*M_PI*n));
    }
}


/* Entropy calculated as ln( n! / \prod n_i! ) where n_i is the number
 * of elements of class i.
 *
 * In other words, the log of the number of possible configurations
 * taking symmetry into account.
 */


static double
entro_from_histogram(const size_t * restrict H,
                     const u32 npoint,
                     const u32 max_label)
{
    double entro = log_p(npoint);

    for(size_t kk = 0; kk < max_label; kk++)
    {
        entro -= log_p(H[kk]);
    }
    return entro;
}


double
entropy_split(const u32 * restrict class,
              const f64 * restrict feature,
              const u32 npoint,
              const u32 max_label,
              u32  * restrict _nleft,
              u32  * restrict _nright,
              f64* restrict _eleft,
              f64* restrict _eright)
{
    if(npoint < 2)
    {
        return 0;
    }

    /* Histogram of the number of elements of each class */
    size_t HL[max_label + 1];
    size_t HR[max_label + 1];

    memset(HL, 0, (max_label+1)*sizeof(size_t));
    memset(HR, 0, (max_label+1)*sizeof(size_t));

    assert(HL[0] == 0);
    /* Initially the "right" histogram contain all points */
    for(size_t kk = 0; kk < npoint; kk++)
    {
        assert(class[kk] <= max_label);
        HR[class[kk]]++;
    }

    /* Base line, all in one partition */
    double entro_min =
        entro_from_histogram(HR, npoint, max_label);
    u32 best_nleft = 0;
    double best_lentro = 0;
    double best_rentro = entro_min;

    //printf("All in one S: %f\n", entro_min);
    for(size_t kk = 1; kk < npoint; kk++)
    {
        // move one point from right to left
        HL[class[kk-1]]++;
        HR[class[kk-1]]--;


        if(feature[kk-1] == feature[kk])
        {
            /* There is no threshold that splits the data with kk
               points to the left. */
            continue;
        }
        double left_entro = entro_from_histogram(HL, kk, max_label);
        double right_entro = entro_from_histogram(HR, npoint-kk, max_label);
#if 0
        printf("nleft: %4u left: %4f, right: %4f, sum: %4f\n",
               kk,
               left_entro, right_entro,
               left_entro+right_entro);
#endif
        if(left_entro + right_entro < entro_min)
        {
            entro_min = left_entro + right_entro;
            best_nleft = kk;
            best_lentro = left_entro;
            best_rentro = right_entro;
        }
    }

    *_nleft = best_nleft;
    *_nright = npoint - best_nleft;
    *_eleft = best_lentro;
    *_eright = best_rentro;

    return entro_min;
}

double
entropy_evaluate(const u32 * class, const u32 npoint, const u32 max_label)
{
    if(npoint < 2)
    {
        return 0;
    }

    /* Histogram of the number of elements of each class */
    size_t H[max_label + 1];
    memset(H, 0, (max_label+1)*sizeof(size_t));

    /* Initially the "right" histogram contain all points */
    for(size_t kk = 0; kk < npoint; kk++)
    {
        assert(class[kk] <= max_label);
        H[class[kk]]++;
    }

    return entro_from_histogram(H, npoint, max_label);
}
