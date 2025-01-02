#pragma once

#include "trafo_util.h"

/* Suggests how to split the vector _class_ containing class labels
 * into two partitions of size _nleft and _nright.
 *
 * Since a split only is possible between points with different
 * values, the feature vector is used for exactly this purpose. For example if
 * feature = [0.1, 0.2, 0.2, 0.3, 0.3] the only possible outcomes are
 * nleft \in {0, 1, 3, 5}.
 *
 * The Gini impurity, $G$, for a partition is defined as
 *
 * $$G = 1 - \sum p_i^2 $$
 *
 * where $p_i$ is the fraction of elements that belongs to class $i$
 *
 * 0, or "no impurity" is achieved if all elements belong to the same
 * class, if the partition has a mix of classes this value will be larger
 *
 * Returns:
 * the return value is the gini of the split.
 * G_left and G_right is gini purity of the left and right partition
 * n_left and n_right is the number of elements in the left and right partition
 *
 */

double
gini_split(const u32 * restrict class,
           const f64 * restrict feature,
           const u32 npoint,
           const u32 max_label,
           u32  * restrict n_left,
           u32  * restrict n_right,
           f64* restrict G_left,
           f64* restrict G_right);

/* Evaluate the gini impurity of a set of class labels */
double gini_evaluate(const u32 * restrict class,
                     const u32 npoint,
                     const u32 max_label);
