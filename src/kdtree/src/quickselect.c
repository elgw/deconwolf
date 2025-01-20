/*
  Copyright 2025 Erik Wernersson

  Permission is hereby granted, free of charge, to any person
  obtaining a copy of this software and associated documentation files
  (the “Software”), to deal in the Software without restriction,
  including without limitation the rights to use, copy, modify, merge,
  publish, distribute, sublicense, and/or sell copies of the Software,
  and to permit persons to whom the Software is furnished to do so,
  subject to the following conditions:

  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
  BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
  ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#include "quickselect.h"

#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifndef QUICKSELECT_F32
#ifndef QUICKSELECT_F64
#define QUICKSELECT_F32
#endif
#endif

#ifdef QUICKSELECT_F32
typedef float etype;
#define FNAME(x) x ## _f32
#endif

#ifdef QUICKSELECT_F64

#ifdef QUICKSELECT_F32
#error "both QUICKSELECT_F32 and QUICKSELECT_F64 are defined"
#endif

typedef double etype;
#define FNAME(x) x ## _f64

#endif

#ifndef FNAME
#error "quickselect precision not defined, use -DQUICKSELECT_FP32 or -DQUICKSELECT_FP64"
#endif

#define qs_max(x,y) (x>y ? x : y)
#define qs_min(x,y) (x<y ? x : y)

static etype
med3(etype a, etype b, etype c)
{
    return qs_max(qs_min(a,b),qs_min(c,qs_max(a,b)));
}

static etype
med5(etype a, etype b, etype c, etype d, etype e)
{
    etype f=qs_max(qs_min(a,b),qs_min(c,d)); // discards lowest from first 4
    etype g=qs_min(qs_max(a,b),qs_max(c,d)); // discards biggest from first 4
    return qs_max(qs_min(e,f),qs_min(g,qs_max(e,f))); /* median of 3 elements */
}

static void
array_minmax(const etype * restrict A, const size_t N,
             etype * restrict min, etype * restrict max)
{
    *min = A[0];
    *max = A[0];
    for(size_t kk = 1; kk < N; kk++)
    {
        *min > A[0] ? *min = A[0] : 0;
        *max < A[0] ? *max = A[0] : 0;
    }
    return;
}

static etype
_quickselect(etype * restrict X, const size_t N, const size_t s,
             const int use_pb, etype ** restrict PB);

static etype get_pivot_random(const etype * X, size_t N)
{
    size_t idx = rand() % N;
    return X[idx];
}

static etype
get_pivot(const etype * restrict X, const size_t N, etype ** restrict PB)
{
    const int use_pb = 1;
    assert(N > 1);
    if(N == 2)
    {
        return X[0];
    }

    if(N < 5)
    {
        return med3(X[0], X[(N-1)/2], X[N-1]);
    }

    size_t nM = N / 5;


    etype * M = PB[0];
    PB[0] += nM; /* Next time, use a higher address */
    for(size_t mm = 0; mm < nM; mm++)
    {
        const etype * P = X+5*mm;
        M[mm] = med5(P[0], P[1], P[2], P[3], P[4]);
    }

    etype pivot = _quickselect(M, nM, nM/2, use_pb, PB);

    PB[0] -= nM; /* Done, return the data */

    return pivot;
}


static void
partition(etype * restrict X, const size_t n,
          const etype pivot,
          size_t * nLow, size_t * nHigh)
{
    assert(n > 0);
    int64_t low = -1;
    int64_t high = n;
    int64_t n2 = n;

    while(1)
    {
        do
        {
            low++;
        } while ( low < n2 && X[low] <= pivot );

        do
        {
            high--;
        } while ( high > 0 && X[high] > pivot );

        if(low >= high)
        {
            *nLow = low;
            *nHigh = n-*nLow;
#ifndef NDEBUG
            assert(low >= 0);
            assert(high < (int64_t) n);
            assert(*nLow + *nHigh == n );
            for(int64_t kk = 0; kk < low; kk++)
            {
                assert(X[kk] <= pivot);
            }
            for(int64_t kk = low; kk < n2; kk++)
            {
                assert(X[kk] > pivot);
            }
            assert(high > -1);
#endif
            return;
        } else {
            // Make a swap
#ifndef NDEBUG
            assert(low >= 0);
            assert(high < (int64_t) n);
            assert(low < high);
#endif
            etype t = X[low];
            X[low] = X[high];
            X[high] = t;
        }
    }

    return;
}

static etype
_quickselect(etype * restrict X, const size_t N, const size_t s,
             const int use_pb, etype ** restrict PB)
{

    assert( s < N );
    assert( N != 0 );
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
            return qs_min(X[0], X[1]);
        } else {
            return qs_max(X[0], X[1]);
        }
    }

    if(N == 3)
    {
        if(s == 0)
        {
            return qs_min(qs_min(X[0], X[1]), X[2]);
        }
        if(s == 1)
        {
            return med3(X[0], X[1], X[2]);
        }
        return qs_max(qs_max(X[0], X[1]), X[2]);
    }

    etype pivot;
    if(use_pb) {
        pivot = get_pivot(X, N, PB);
    } else {
        pivot = med5(X[0], X[(N-1)/4],
                     X[(N-1)/2], X[3*(N-1)/4], X[N-1]);
    }

    /* Solve the dutch flag problem,
     *  split into
     *  A: Lower than the pivot,
     *  B: equal to the pivot
     *  C Higher than the pivot
     * However we skip the B partition ...
     */

    size_t nA=0, nC=0;
    partition(X, N, pivot, &nA, &nC);

    /* Avoid infinite recursion */
    if(nA == 0 || nC == 0)
    {
        pivot = get_pivot_random(X, N);
        partition(X, N, pivot, &nA, &nC);
    }

    /* Possibly constant array, fallback */
    if(nA == 0 || nC == 0)
    {
        etype min_value, max_value;
        array_minmax(X, N, &min_value, &max_value);

        if(min_value == max_value)
        {
            return min_value;
        } else {
            pivot = (min_value + max_value) / 2;
            partition(X, N, pivot, &nA, &nC);
        }
    }

    assert(nA != 0);
    assert(nC != 0);
    assert(nA + nC == N);

    etype *A = X;
    etype *C = X+nA;

#ifndef NDEBUG
    for(size_t kk = 0; kk<nA; kk++)
    { assert(A[kk] <= pivot); }
    for(size_t kk = 0; kk<nC; kk++)
    { assert(C[kk] > pivot); }
#endif

    /* Recurse on the partition that can contain the median */

    if(s < nA)
    {
        return _quickselect(A, nA, s, use_pb, PB);
    }

    return _quickselect(C, nC, s-nA, use_pb, PB);
}


static etype
quickselect_start(const etype * X, size_t N, size_t s, int use_pb)
{
    assert(s < N);

    if(N == 1)
    {
        assert(s == 0);
        return(X[0]);
    }

    etype * X2 = calloc(N, sizeof(etype));
    assert(X2 != NULL);
    memcpy(X2, X, N*sizeof(etype));


    etype * PB = NULL;
    if(use_pb) {
        PB = calloc(N/2, sizeof(etype)); /* Pivot buffer */
        assert(PB != NULL);
    }

    etype element = _quickselect(X2, N, s, use_pb, &PB);

    if(use_pb){
        free(PB);
    }

    free(X2);
    return element;
}

etype FNAME(qselect)(const etype * X, size_t N, size_t s)
{
    if( X == NULL )
    {
        return 0;
    }
    return quickselect_start(X, N, s, 0);
}
