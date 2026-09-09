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
#include <stdio.h>

typedef int64_t i64;

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

#define qs_max(x,y) ((x>y) ? x : y)
#define qs_min(x,y) ((x<y) ? x : y)

static etype
med3(etype a, etype b, etype c)
{
    return qs_max(qs_min(a,b),qs_min(c,qs_max(a,b)));
}

static inline etype
med5(const etype * X)
{
    etype f=qs_max(qs_min(X[0],X[1]),qs_min(X[2],X[3])); // discards lowest from first 4
    etype g=qs_min(qs_max(X[0],X[0]),qs_max(X[2],X[3])); // discards biggest from first 4
    return qs_max(qs_min(X[4],f),qs_min(g,qs_max(X[4],f))); /* median of 3 elements */
}

static inline etype
min5(const etype * X)
{
    etype min = X[0];
    for(int kk = 1; kk < 5; kk++){
        X[kk] < min ? min = X[kk] : 0;
    }
    return min;
}

static inline etype
max5(const etype * X)
{
    etype max = X[0];
    for(int kk = 1; kk < 5; kk++){
        X[kk] > max ? max = X[kk] : 0;
    }
    return max;
}


static void
array_minmax(const etype * restrict A, const size_t N,
             etype * restrict min, etype * restrict max)
{
    *min = A[0];
    *max = A[0];
    for(size_t kk = 1; kk < N; kk++) {
        *min > A[kk] ? *min = A[kk] : 0;
        *max < A[kk] ? *max = A[kk] : 0;
    }
    return;
}

#if 0
#define SWAPE(a, b) etype t = a; a = b; b = t;
static void sort4(etype X[4])
{
    if(X[0] > X[2]){SWAPE(X[0], X[2])};
    if(X[1] > X[3]){SWAPE(X[1], X[3])};
    if(X[0] > X[1]){SWAPE(X[0], X[1])};
    if(X[2] > X[3]){SWAPE(X[2], X[3])};
    if(X[1] > X[2]){SWAPE(X[1], X[2])};
}
#endif

static etype
_quickselect(etype * restrict X, const size_t N, const size_t s);

static void
partition(etype * restrict X, const size_t n,
          const etype pivot,
          size_t * nLow, size_t * nHigh)
{
    assert(n > 0);
    i64 low = -1;
    i64 high = (i64) n;
    i64 n2 = (i64) n;

    while(1){
        do {
            low++;
        } while ( low < n2 && X[low] <= pivot );

        do {
            high--;
        } while ( high > 0 && X[high] > pivot );

        if(low >= high) {
            goto done;
        } else { // swap
            assert(low < high);
            etype t = X[low];
            X[low] = X[high];
            X[high] = t;
        }
    }

done:
    *nLow = (size_t) low;
    *nHigh = n - (size_t) low;

#ifdef QUICKSELECT_DEBUG
    assert(low >= 0);
    assert(high < (i64) n);
    assert(*nLow + *nHigh == n );
    for(i64 kk = 0; kk < low; kk++) {
        assert(X[kk] <= pivot);
    }
    for(i64 kk = low; kk < n2; kk++) {
        assert(X[kk] > pivot);
    }
    assert(high > -1);
#endif
    return;
}


// Returns 1 if a pivot could be found
// or 0 if all elements are equal, which means that we found our number
static int
get_pivot(const etype * X, size_t N, etype * pivot)
{
    // Strategy 1
    // Something more dynamic would be compelling.
    // please benchmark on much data!
    etype vals[5];
    if(N < 1024){ // median of 5 maximally spread numbers
        for(size_t kk = 0; kk < 5; kk++){
            vals[kk] = X[kk*(N-1)/4];
        }
    } else { // med5(med5, med5, ...)
        for(size_t kk = 0; kk < 5; kk++){
            etype lvals[5];
            for(size_t ll = 0; ll < 5; ll++){
                size_t ii = kk*5 + ll;
                lvals[ll] = X[ii*(N-1)/25];
            }
            vals[kk] = med5(lvals);
        }
    }

    *pivot = med5(vals);

    if(! (*pivot == min5(vals) || *pivot == max5(vals))){
        return 0;
    }

#if 0
    // Strategy 2: median of 5 values from random locations
    // Probably better to do kk*p % N, where gcd(p, N) == 1
    // and p large
    if(N > 100)
    {
        for(int kk = 0; kk < 5; kk++){
            vals[0] = X[(size_t) rand() % N];
        }

        *pivot = med5(vals);
        if(! (*pivot == min5(vals) || *pivot == max5(vals))){
            return 0;
        }
    }
#endif
    // Last resort, check if all numbers are equal
    etype min_value, max_value;
    array_minmax(X, N, &min_value, &max_value);
    if(min_value == max_value) {
        *pivot = max_value;
        return 1; // all equal and we are done
    }
    // If not all equal, select something in-between
    *pivot = (min_value + max_value) / (etype) 2.0;
    // there might not by any value between them
    // nextafter(min_value, max_value) == max_value
    // so we simply rely on this
    if(*pivot >= max_value || *pivot <= min_value){
        *pivot = min_value;
    }
    return 0;
}

static etype
_quickselect(etype * restrict X, const size_t N, const size_t s)
{
    assert( s < N );
    assert( N != 0 );
    /* Only one element left. Has to be the one we are looking for */
    if(N == 1){
        return(X[0]);
    }

    if(N == 2){
        if(s == 0) {
            return qs_min(X[0], X[1]);
        } else {
            return qs_max(X[0], X[1]);
        }
    }

    if(N == 3) {
        if(s == 0) {
            return qs_min(qs_min(X[0], X[1]), X[2]);
        }
        if(s == 1) {
            return med3(X[0], X[1], X[2]);
        }
        return qs_max(qs_max(X[0], X[1]), X[2]);
    }

#if 0 // would make it slower
    if(N == 4){
        sort4(X);
        return X[s];
    }
#endif

    etype pivot = X[0];
    if(get_pivot(X, N, &pivot)){
        // All elements equal, no pivot could be found
        return pivot;
    }

    //  split into
    //  A: Lower than the pivot,
    //  C: Higher or equal to the pivot
    size_t nA=0, nC=0;
    partition(X, N, pivot, &nA, &nC);
    if(nA + nC != N){
        fprintf(stderr, "Internal error, L%d\n", __LINE__);
        return 0;
    }

    if(nA == 0){
        fprintf(stderr, "Internal error, L%d\n", __LINE__);
        return 0;
    }
    if(nC == 0){
        fprintf(stderr, "Internal error, L%d\n", __LINE__);
        return 0;
    }

    etype *A = X;
    etype *C = X+nA;

#ifdef QUICKSELECT_DEBUG
    for(size_t kk = 0; kk<nA; kk++)
    { assert(A[kk] <= pivot); }
    for(size_t kk = 0; kk<nC; kk++)
    { assert(C[kk] > pivot); }
#endif

    // Recurse on the partition that can contain element with index s
    if(s < nA) {
        return _quickselect(A, nA, s);
    } else {
        return _quickselect(C, nC, s-nA);
    }
}

etype
FNAME(qselect)(const etype * X, size_t N, size_t s)
{
    if(X == NULL){ return 0; }
    if(s >= N){ return 0; }

    etype * X2 = malloc(N*sizeof(etype));
    assert(X2 != NULL);
    memcpy(X2, X, N*sizeof(etype));
    etype element = _quickselect(X2, N, s);
    free(X2);
    return element;
}
