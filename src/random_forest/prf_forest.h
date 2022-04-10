#ifndef __prf_forest_h__
#define __prf_forest_h__

#include "prf_tree.h"
#include "prf_util.h"

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

/* A random forest implementation
*/
typedef struct PrfForest{
    size_t ntrees; // Number of trees
    // Number of samples to use when constructing a tree
    int samples_per_tree; // settable
    // Number of features to use when constructing a gree
    int min_node_size; // settable
    int features_per_tree; // settable
    int nthreads; // settable
    PrfTree ** trees;
    int verbose; // 0 = quiest, 1 = some, 2 = ...
    int64_t maxClass;
    size_t ** H; // Histogram buffers, one per thread
    size_t nH_allocated; // Number of allocated histogram buffers
} PrfForest;

/* Create a new random forest
 * Initially empty until trained.
 * Set the appropriate parameters before training.
 */

PrfForest * prf_forest_new(size_t ntrees);


// Copy settings from one tree to another
void prf_forest_copy_settings(PrfForest *, const PrfForest *);

/* Train a new forest based on the data in X
 * X: The input data in X should be in column-major format,
 *    i.e., first all data for feature0, then all data for feature1, ...
 *    ending with the classes as the last N values.
 *
 * N: the number of samples
 * M: the number of columns, i.e., number of features + 1 (for class)
 *
 * Returns EXIT_FAILURE or EXIT_SUCCESS
 */
int prf_forest_train(PrfForest * F,
                      const float * X,
                      const size_t N,
                      const size_t M);

int prf_forest_classify(PrfForest *, const float *);
int * prf_forest_classify_table(PrfForest * F, float * X, size_t N, size_t M);

// Print a summary of the settings
void prf_forest_print(FILE * f, PrfForest *);

/* Free a forest, i.e. everything that was allocated with prf_forest_new
 * and prf_forest_train
 */

void prf_forest_free(PrfForest *);


/* Perform k-fold cross validation on the table X.
 * If F is supplied, the settings in F are used.
 * Returns the average percentage of correctly classified
 * samples.  */
float prf_forest_cross_validate_k(const PrfForest * F,
                                   const float * X, const size_t N, const size_t M,
                                   const int k);

/* For each feature, determine how important it is by replacing it by
 * random values. A[kk] tells the accuracy when feature kk is replaced by
 * random values. A[M-1] tells the accuracy when all features are used. */
void prf_forest_feature_importance(const PrfForest * F,
                              const float * X, const size_t N, const size_t M,
                              float ** FI);

// FORWARD DECLARATIONS
void prf_forest_free_H(PrfForest *);
void prf_forest_init_H(PrfForest *);

// TODO /* Alternative representation, as a simpler table */


typedef struct{
    size_t Ntrees;
    PrfTreeTable * tree_tables;
} PrfForestTable;

void prf_forest_table_to_file(PrfForestTable * FT, char * filename);
PrfForestTable * prf_forest_table_from_file(char * filename);
int prf_forest_table_classify(PrfForestTable * FT, float * X);
void prf_forest_table_to_code(PrfForestTable * FT, char * filename);


#endif
