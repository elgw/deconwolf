#pragma once
#include <stddef.h>

/** @brief quickselect
 * @param X the data points which will be scrambled
 * @param nX number of elements of X
 * @param s what element to select from sorted X
 * @return the value of element s of sorted X.
 */
double quickselect(double * X, size_t nX, size_t s);
