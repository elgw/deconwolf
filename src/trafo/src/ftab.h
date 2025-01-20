#pragma once

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

/* Floating point-only table stored in row-major format.
 * Original repository: github.com/elgw/ftab/
 *
 * The uggly:
 *
 * - Missing values will be parsed as 0
 * - Always interprets the first line as a header
 *
 *
 * TODO
 *
 * Use strtof instead of atof
 * Possibly support writing using %a for exactness
 * Decide what to do about missing values
 * Handle NAN and Inf
 * Parsing header or not option.
 * Option to ignoring comment lines starting with #
 *
 *
 * CHANGELOG
 *
 * 0.1.1 : switched from strtok to strsep to handle also empty values
 * 0.1.2 : added ftab_compare to compare two tables.
 */

#define FTAB_VERSION_MAJOR 0
#define FTAB_VERSION_MINOR 1
#define FTAB_VERSION_PATCH 2

#include <stdint.h>
#include <stdio.h>


/* row-major table */
typedef struct {
    /* Pointer to the table data. Please note that the address can change between calls to the API */
    float * T;
    /* Number of rows */
    size_t nrow;
    /* Number of columns */
    size_t ncol;
    size_t nrow_alloc; /* To know if we need to extend the size */
    char ** colnames; /* Name of columns can be NULL. Also the pointer can be NULL */
} ftab_t;

/* Create a new table with a fixed number of columns
 * Set column names with ftab_set_colname */
ftab_t * ftab_new(int ncol);

/* Create a new table from raw data. The data has to be in row major format */
ftab_t * ftab_new_from_data(int nrow, int ncol, const float * data);


/* Load a TSV file. The first line is interpreted as
 * containing the column names. Everything else is interpreted
 * as float values.
 */
ftab_t * ftab_from_tsv(const char * fname);

ftab_t * ftab_from_csv(const char * fname);

/* Write tsv file do disk */
int ftab_write_tsv(const ftab_t * T, const char * fname);

/* Write tsv file do disk */
int ftab_write_csv(const ftab_t * T, const char * fname);

/* Write tsv file do disk */
int ftab_write_csv(const ftab_t * T, const char * fname);

/** Print table to file
 * @param[in] fid An open FILE to write to
 * @param[in] T the table to write
 * @param[in] sep The separator or deliminator to use, e.g., "," or "\t"
 * @return EXIT_SUCCESS if file could be written
 */
int ftab_print(FILE * fid, const ftab_t * T, const char * sep);

/* Set the name of a column
* the name is copied, i.e., can be freed */
void ftab_set_colname(ftab_t *, int col, const char * name);

/* Free a ftab and all associated data */
void ftab_free(ftab_t * T);

/* Append a single row. */
void ftab_insert(ftab_t * T, float * row);

/* Get the index of a certain column name
 * Returns -1 on failure. Undefined behavior if
 * multiple columns have the same name. */
int ftab_get_col(const ftab_t * T, const char * name);


/** @brief Set the data for one column.
 * @param T: table to receive data
 * @param col: column to write to
 * @param data: pointer to data to insert
 */

int ftab_set_coldata(ftab_t * T, int col, const float * data);

/** @brief Horizontal concatenation
 *
 * @param L : data on the left side
 * @param R : data on the right side
*/
ftab_t * ftab_concatenate_columns(const ftab_t * L, const ftab_t * R);

/* Concatenate two tables vertically with T on the top and B on the bottom */
ftab_t * ftab_concatenate_rows(const ftab_t * T, const ftab_t * B);

/* Subselect rows where row_selector > 0 */
void ftab_subselect_rows(ftab_t * T, const uint8_t * row_selector);

/* Keep n head rows */
void ftab_head(ftab_t * T, int64_t n);

/* Subselect rows where row_selector > 0
 * The row_selector array needs to have as many elements as there are rows.
 * The table is modified.
 */
void ftab_subselect_rows(ftab_t * T, const uint8_t * row_selector);

/** Create a deep copy */
ftab_t * ftab_copy(const ftab_t * T);

/* Compare two tables. Returns 0 if the are equal returns a positive
* number if either or both pointer are null or if the tables are
* different
* Column names are also compared
*/
int ftab_compare(const ftab_t *, const ftab_t * );

/* Run some unit tests */
int ftab_ut(int argc, char ** argv);
