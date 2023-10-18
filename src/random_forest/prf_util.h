#ifndef __prf_util_h_
#define __prf_util_h_

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include "../dw_util.h"

// For Column-Major

// Show row n from the table T with NxM items
void show_row_n(const float * T , size_t N, size_t M, size_t n);

// Index of element with largest value
size_t vector_size_t_argmax(const size_t * V, size_t N);


// Returns a binary list where n elements are set to 1
// at random. Useful for random subselection without replacement
uint8_t * get_random_selection(size_t N, size_t n);

// Shuffle a row-major table
void shuffle_rm_table(float * X, size_t N, size_t M);

// Get training and validation sets for k-fold cross-validation
void get_subset_k(const float * T, const size_t N, const size_t M,
                 const int K, const int k,
                 float ** Tr, size_t * nTr,
                 float ** Va, size_t * nVa);

/* For row-major tables,
 * Copy M elements from table B to table A
 * the stride (number of rows) is specified by sA and sB
 */
void copy_row(float * A, const size_t sA,
              const float * B, const size_t sB,
              const size_t M);

/* Replace one feature by random values */
void scramble_feature(float * X, const size_t N, const size_t M, const int f);

void vector_show(float * v, int n);

/* Print a section header to the terminal */
void print_section(char * msg);
void print_ok(void);
void print_todo(char * msg);
void print_warning(char * msg);
void print_error(char * msg);
#endif
