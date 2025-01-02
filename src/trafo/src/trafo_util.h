#pragma once

#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

typedef uint8_t u8;
typedef uint32_t u32;
typedef float f32;
typedef double f64;
typedef int32_t i32;
typedef int32_t i64;


double
timespec_diff(struct timespec* end, struct timespec * start);


int
get_peakMemoryKB(size_t * VmPeak, size_t * VmHWM);

void print_peak_memory(void);

void
print_section(const char * msg);
