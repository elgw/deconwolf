#pragma once

#include <assert.h>
#include <math.h> // tgmath.h would be an alternative
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef u64
typedef uint64_t u64;
#endif

float * sparse_preprocess(const float * image, u64 M, u64 N, u64 P,
                          double lambda, double lambda_s,
                          int periodic, int directions, u64 iter,
                          int verbose, FILE * log);
