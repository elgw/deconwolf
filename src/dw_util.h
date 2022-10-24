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

#ifndef __dw_util_h__
#define __dw_util_h__

#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#define tictoc struct timespec tictoc_start, tictoc_end;
#define tic clock_gettime(CLOCK_REALTIME, &tictoc_start);
#define toc(X) clock_gettime(CLOCK_REALTIME, &tictoc_end); printf(#X); printf(" %f s\n", timespec_diff(&tictoc_end, &tictoc_start)); fflush(stdout);


/* Create dir if it does not exist.
 * Returns 0 if the dir already existed or could be created
 * returns non-zeros if the dir can't be created
 */
int ensuredir(char * dir);

/* Check if directory exist, do not create if missing
 * returns 1 if it exist
 * */

int isdir(char * dir);

/* Return a suggestion for how many threads to use
 * ideally this should be the same as the number of cores
 */
int dw_get_threads(void);

/* Return 1 if the file exist, else 0 */
int dw_file_exist(char * fname);

/* Free if not null */
void dw_nullfree(void * p);

/* Difference between two timepoints */
float timespec_diff(struct timespec* end, struct timespec * start);

/* Get peak memory usage, works on Linux and MacOS
* returns 0 under windows */
size_t get_peakMemoryKB(void);

/* Read the scaling of file from the .log.txt file if exists */
float dw_read_scaling(char * file);

#endif
