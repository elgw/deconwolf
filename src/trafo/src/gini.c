#include <math.h>
#include "gini.h"

/* Calculates
 *
 * $$ \sum \left( \frac{V_i}{T} \right)^2 $$
 *
 * Here it is assumed that $$\sum V_i == T$$
 */
static double
squaresum(const size_t * V, size_t T, size_t N)
{
    double s = 0;
    for(size_t kk = 0; kk<N; kk++)
    {
        assert( V[kk] <= T );
        s+= pow((double) V[kk]/ (double) T, 2);
    }
    assert(s <= 1);
    return s;
}


static double
gini_from_2H(const size_t * HL,
             const size_t * HR,
             size_t nclass,
             const size_t nL,
             const size_t nR,
             double * lgini, double * rgini)
{
    double gL = 1.0 - squaresum(HL, nL, nclass);
    double gR = 1.0 - squaresum(HR, nR, nclass);
    double N = nL + nR;
    *lgini = gL;
    *rgini = gR;

    return  (double) nL/N * gL + (double) nR/N * gR;
}

double
gini_split(const u32 * restrict class,
           const f64 * restrict feature,
           const u32 npoint,
           const u32 max_label,
           u32  * restrict _nleft,
           u32  * restrict _nright,
           f64* restrict _gleft,
           f64* restrict _gright)
{
    if(npoint  < 2)
    {
        return 0;
    }
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

    double best_lgini, best_rgini;
    double mingini = gini_from_2H(HR, HR, max_label,
                                  npoint, npoint,
                                  &best_lgini, &best_rgini);

    u32 best_nleft = 0;

    for(size_t kk = 1; kk < npoint; kk++)
    {
        HL[class[kk-1]]++;
        HR[class[kk-1]]--;

        if(feature[kk-1] == feature[kk])
        {
            /* There is no threshold that splits the data with kk
               points to the left. */
            continue;
        }
        double lgini, rgini;
        double g = gini_from_2H(HL, HR, max_label,
                                kk, npoint-kk,
                                &lgini, &rgini);
        if(g < 0)
        {
            printf("g = %f (%f, %f)\n", g, lgini, rgini);
            printf("left: %zu, right: %zu\n", kk, npoint-kk);
            for(size_t ll = 0; ll <= max_label; ll++)
            {
                printf("%zu, %zu\n", HL[ll], HR[ll]);
            }
        }
        assert(g>=0);

        if(g < mingini)
        {
            best_nleft = kk;
            mingini = g;
            best_lgini = lgini;
            best_rgini = rgini;
        }
    }

    /* Note: The case when all points are to the left
       is identical to the case where all points are to the right
       and not considered */
    *_gleft = best_lgini;
    *_gright = best_rgini;
    *_nleft = best_nleft;
    *_nright = npoint - best_nleft;
    return mingini;
}

double gini_evaluate(const u32 * restrict class,
                     const u32 npoint,
                     const u32 max_label)
{
    if(npoint  < 2)
    {
        return 0;
    }

    size_t H[max_label + 1];

    memset(H, 0, (max_label+1)*sizeof(size_t));

    /* Initially the "right" histogram contain all points */
    for(size_t kk = 0; kk < npoint; kk++)
    {
        assert(class[kk] <= max_label);
        H[class[kk]]++;
    }

    return 1.0 - squaresum(H, npoint, max_label);
}
