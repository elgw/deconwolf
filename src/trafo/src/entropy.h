#pragma once

#include "trafo_util.h"


/* Entropy is defined as the log of the number of possible configurations
*
* S = log W
*
* hence, if there are two partions the entropy is defined as the sum
* of their entropies,
*
* S = log W_1 + log W_2
*
*/

/* Evaluate the entropy of a set of labels */
double
entropy_evaluate(const u32 * class, const u32 npoint, const u32 max_label);

/* Suggests how to split the vector _class_ containing class labels
 * into two partitions of size _nleft and _nright.
 *
 * Since a split only is possible between points with different
 * values, the feature vector is used for exactly this purpose. For example if
 * feature = [0.1, 0.2, 0.2, 0.3, 0.3] the only possible outcomes are
 * nleft \in {0, 1, 3, 5}.
 *
 * _eleft and _eright is the entropy of the left and right partition.
*/

double
entropy_split(const u32 * restrict class,
              const f64 * restrict feature,
              const u32 npoint,
              const u32 max_label,
              u32  * restrict _nleft,
              u32  * restrict _nright,
              f64* restrict _eleft,
              f64* restrict _eright);
