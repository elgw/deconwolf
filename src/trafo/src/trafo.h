#pragma once

#include <stdint.h>

#define TRAFO_VERSION_MAJOR 0
#define TRAFO_VERSION_MINOR 1
#define TRAFO_VERSION_PATCH 5


/* These are the configuration options for the Random Forest
 * fitting/creation
 *
 * The input data, i.e. the features should be supplied setting either
 * F_col_major or F_row_major depending on the memory layout.
 *
 * n_tree: The number of trees to construct.
 * n_feature: The number of features
 * n_sample: The number of samples
 *
 * tree_n_sample: Says how many samples that should be used when
 * constructing each tree. If not set (= 0) 0.623*n_sample will
 * be used.
 *
 * tree_n_feature: How many feature to sample for each tree. If not
 * set (=0) sqrt(n_feature) will be used.
 *
 * verbose:
 * The verbosity level. 0=minimal, 1=some, 2 or above is only
 * for debugging
 *
 **/

/* To do: enable support for building for float */
#ifdef TRAFO_FLOAT32
#define C(x) x##f
#define fpnumber float
#else
#define C(x) x
#define fpnumber double
#endif

typedef struct {
    const fpnumber * F_col_major;
    const fpnumber * F_row_major;
    const uint32_t * label;
    uint32_t n_tree;
    uint32_t n_sample;
    uint32_t n_feature;

    /* Optional settings that have defaults */

    /* Minimum number of samples per node which converts it into a
     * leaf. I.e. nodes are not split if their number of samples are
     * less or equal to min_samples_leaf
     *
     * Default: 1 will be used if this is set to 0.
     */
    uint32_t min_samples_leaf;

    /* Fraction of samples to use for each tree
     *
     * Valid settings: ]0, 1]. The actual number of samples per tree
     * will be rounded to be in the integer interval [1, n_sample]
     *
     * Default: will be set to 0.632 if set to 0
     */
    float tree_f_sample;

    /* Number of feature per tree.
     *
     * Default: sqrt(n_feature) will be used if this parameter is set to 0
     */
    uint32_t tree_n_feature;
    uint32_t verbose;
    int entropy; // Split criterion: 0 = gini, 1 = entropy

} trafo_settings;

typedef struct _trf trf;

/* Train and return a forest using the supplied settings. Both the
 * settings and all the data that it points to can be freed after this
 * call as the function copies what it needs.
 *
 * Return NULL if the settings are invalid or on failure
 */
trf * trafo_fit(trafo_settings * settings);

/* Free up all memory associated with T */
void trafo_free(trf * T);

/* Predict the class for each point/sample in X. It is the
 * responsibility of the caller to free the returned array.
 *
 * Could return NULL on failure.
 **/
uint32_t * trafo_predict(trf * T,
                       const fpnumber * X_cm,
                       const fpnumber * X_rm,
                       uint64_t n_point);

/* Print out a summary of the settings */
void
trafo_print(FILE * fid, const trf * s);

/* Return a vector with the importance for each feature. Returns NULL
 * if this information is not available.
 *
 * When splitting by Gini
 * at each successful split the importance is increased by
 * importance(feature) += n*G_0 - (n_left+G_left + n_right*G_right)
 *
 * When splitting by Entroty
 * at each successful split the importance is increased by
 * importance(feature) += E_0 - (E_left + E_right)
 *
 * Please note that the feature importance will be diluted for forests
 * since each tree normally contain only a subset of the features
 *
 * The caller is responsible for freeing the returned memory
 */
fpnumber *
trafo_importance(trf * T);


/* Save to disk */
int trafo_save(trf * s,
             const char * filename);

/* Load from disk */
trf *
trafo_load(const char * filename);

/* Run some unit tests */
int
trafo_ut(void);
