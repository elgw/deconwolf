#pragma once
#include <stdint.h>
#include <stddef.h>
#include "qsort.h"

/* maintain sorted feature vectors by auxiliary rank vectors.
* It would be possible to do the same using order vectors as well.
* with similar tradeofs.
*
* To update ranks could be delayed until needed.
* Some temp buffers could be shared ...
* gprof it !
*
* Splitting could also be delayed until needed?
*/

typedef struct _sortbox sortbox;

/* F and classes can be freed afterwards since the sortbox will copy
   them
 */
sortbox *
sortbox_init(const double * F,
             const uint32_t * classes,
             size_t nsamples, size_t nfeatures);

void
sortbox_free(sortbox*);


/* Returns pointers to the feature of interest -- which will be sorted
 * as well as the corresponding class vector. The called is responsible
 * for freeing the arrays. */
void
sortbox_get_feature(const sortbox *,
                    uint32_t feature,
                    uint32_t start,
                    uint32_t len,
                    double ** ret_feature,
                    uint32_t ** ret_class);


/* Splits all vectors in the sortbox where
 * where F[:, feature] < threshold
 *
 * Sets the number of elements to the left and
 * to the right of the threshold value.
 *
 * This operates on the start:start+len segment of each array.
 */

void sortbox_split(sortbox *,
                   const int feature_id,
                   const double th,
                   // Index of first element
                   uint32_t start,
                   // Number of elements
                   uint32_t len,
                   // Return values:
                   uint32_t * ret_nleft,
                   uint32_t * ret_nright);

/* Same as above but splits at a certain number of elements to the
 * left for the given feature. Avoids ambiguities in thresholds. */
void sortbox_split_n(sortbox * P,
                     const int feature_id,
                     uint32_t nleft,
                     // Index of first element
                     uint32_t start,
                     // Number of elements
                     uint32_t len,
                     // Return values:
                     uint32_t * _nleft, uint32_t * _nright);

sortbox *
sortbox_clone(sortbox * B);

/* Subselect samples and features features
 * i.e. subselect for a "random" Tree construction.
*/
sortbox *
sortbox_subsample(const sortbox *,
                 size_t n_samples,
                 size_t n_features);

void sortbox_free(sortbox * B);

size_t
sortbox_get_nsample(const sortbox * B);

size_t
sortbox_get_nfeature(const sortbox * B);

/* Return the number of classes in the sortbox */
uint32_t
sortbox_get_nclass(const sortbox * B);

/* Return the original feature id from a subsampled feature id */
uint32_t
sortbox_map_feature(const sortbox * B, uint32_t id);

const uint32_t *
sortbox_get_class_array_unsorted(const sortbox * B);
