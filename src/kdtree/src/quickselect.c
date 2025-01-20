#include "quickselect.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#define max(x,y) (x>y ? x : y)
#define min(x,y) (x<y ? x : y)

static double med5(double a, double b, double c, double d, double e)
{
    double f=max(min(a,b),min(c,d)); // discards lowest from first 4
    double g=min(max(a,b),max(c,d)); // discards biggest from first 4
    return max(min(e,f),min(g,max(e,f))); /* median of 3 elements */
}

static double med3(double a, double b, double c)
{
    return max(min(a,b),min(c,max(a,b)));
}

static void
partition(double * restrict X, const size_t n,
          const double pivot,
          size_t * nLow, size_t * nHigh)
{
    int64_t low = -1;
    int64_t high = n;
    int64_t n2 = n;
    while(1)
    {
        do { low++; } while ( low < n2 && X[low] <= pivot );

        do { high--; } while ( high > 0 && X[high] > pivot);

        if(low >= high)
        { *nLow = low;  *nHigh = n-*nLow;
#ifndef NDEBUG
            assert(*nLow + *nHigh == n );
            for(int64_t kk = 0; kk < low; kk++)
            {
                assert(X[kk] <= pivot);
            }
            for(int64_t kk = low; kk < n2; kk++)
            {
                assert(X[kk] > pivot);
            }
#endif
            return;
        }
        double t = X[low];
        X[low] = X[high];
        X[high] = t;
    }

    return;
}


static double _quickselect(double * restrict X,
                           const size_t N, const size_t s,
                           size_t recursion_depth, size_t last_N)
{

    assert(N!=0);

    /* Only one element left. Has to be the one we are looking for */
    if(N == 1)
    {
        assert(s == 0);
        return(X[0]);
    }

    if(N == 2)
    {
        if(s == 0)
        {
            return min(X[0], X[1]);
        } else {
            return max(X[0], X[1]);
        }
    }

    if(N == 3)
    {
        if(s == 0)
        {
            return min(min(X[0], X[1]), X[2]);
        }
        if(s == 1)
        {
            return med3(X[0], X[1], X[2]);
        }
        return max(max(X[0], X[1]), X[2]);
    }
    //printf("_quickselect(X, %zu, %zu)\n", N, s);


    double  pivot = med5(X[0], X[(N-1)/4],
                     X[(N-1)/2], X[3*(N-1)/4], X[N-1]);

    //printf("pivot = X[%zu] = %f\n", pivot_idx, pivot);

    /* Solve the dutch flag problem,
     *  split into
     *  A: Lower than the pivot,
     *  B: equal to the pivot
     *  C Higher than the pivot
     */

    size_t nA=0, nC=0;
    //dutch_flag(X, N, pivot, &nA,  &nB, &nC);
    partition(X, N, pivot, &nA, &nC);
    assert(nA + nC == N);
    double *A = X;
    double *C = X+nA;

#ifndef NDEBUG
    for(size_t kk = 0; kk<nA; kk++)
    {
        assert(A[kk] <= pivot);
    }
    for(size_t kk = 0; kk<nC; kk++)
    { assert(C[kk] > pivot); }
#endif

    /* If we could not zoom in any in two iterations
    * it is likely that the remaing data contains only duplicates.
    * to avoid infinite recursion we end here. */
    if(nA == 0 || nA == N)
    {
        if(last_N == N)
        {
            return X[0];
        }
    }

    /* Which of the partitions contain the median? */
    if(s < nA)
    {
        assert(nA != 0);
        return _quickselect(A, nA, s, recursion_depth+1, N);
    }

    if(nC == 0)
    {
        return X[N-1];
    }
    return _quickselect(C, nC, s-nA, recursion_depth, N);
}

double quickselect(double * X, size_t N, size_t s)
{
    return _quickselect(X, N, s, 0, 2*N);
}
