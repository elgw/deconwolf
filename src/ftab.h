/*    Copyright (C) 2020 Erik L. G. Wernersson
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


/* Floating point-only table */

#ifndef __ftab_h__
#define __ftab_h__

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>


/* row-major table */
typedef struct {
    float * T;
    size_t nrow;
    size_t ncol;
    size_t nrow_alloc; /* To know if we need to extend the size */
    char ** colnames; /* Name of columns */
} ftab_t;

/* Create a new table with a fixed number of columns
 * Set column names with ftab_set_colname */
ftab_t * ftab_new(int ncol);

/* Set the name of a column
* name can be freed */
void ftab_set_col_name(ftab_t *, int col, const char * name);

/* Load a TSV file */
ftab_t * ftab_from_tsv(char * fname);

/* Free a ftab and all associated data */
void ftab_free(ftab_t * T);

/* Append a row */
void ftab_insert(ftab_t * T, float * row);

/* Get the index of a certain column name
 * Returns -1 on failure */
int ftab_get_col(ftab_t * T, char * name);

#endif
