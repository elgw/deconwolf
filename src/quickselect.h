#pragma once
// Copyright 2025 Erik Wernersson
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation files
// (the “Software”), to deal in the Software without restriction,
// including without limitation the rights to use, copy, modify, merge,
// publish, distribute, sublicense, and/or sell copies of the Software,
// and to permit persons to whom the Software is furnished to do so,
// subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//  SOFTWARE.

#define ELGW_QS_VERSION_MAJOR 1
#define ELGW_QS_VERSION_MINOR 0
#define ELGW_QS_VERSION_PATCH 3

// An implementation of the quicksort algorithm.
//
// Only one of the functions declared below will be built at a time,
// the default is for double precision. For single precision compile
// with -DQUICKSELECT_F32
//
// Version history
//
// 1.0.3
// - Fix: pivot finding would cause infinite loops when
//   the min and max values were very similar.
//
// 1.0.2
// - Fix: bug in array_minmax
// - better test suite
//
// 1.0.1
// - renames min and max macros to qs_min and qs_max since min and
//   max might already be defined elsewhere.
//
// Web page:
// https://github.com/elgw/quickselect/


#include <stddef.h>

// Use to turn on expensive internal debugging
// #define QUICKSELECT_DEBUG


// A quickselect routine
//
// if XS is a sorted copy of the array X (smallest first)
// then
// qselect_f64(X, nX, s) == XS[s]
//
// Arguments:
//
//  X the data points
//  nX number of elements of X
//  s what element to select from sorted X
//
// If nX == 0, or X==NULL the function will return 0
double
qselect_f64(const double * X, size_t nX, size_t s);

float
qselect_f32(const float * X, size_t nX, size_t s);
