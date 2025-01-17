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

#include <assert.h>
#ifndef WINDOWS
#include <dirent.h>
#endif
#include <errno.h>
#include <fcntl.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#ifdef WINDOWS
#include <io.h>
#include <tchar.h>
#include <direct.h>
#else
#include <libgen.h>
#include <unistd.h>
#endif

#include <string.h>


#define tictoc struct timespec tictoc_start, tictoc_end;
#define tic dw_gettime(&tictoc_start);
#define toc(X) dw_gettime(&tictoc_end); printf(#X); printf(" %f s\n", timespec_diff(&tictoc_end, &tictoc_start)); fflush(stdout);

#ifdef WINDOWS
#define FILESEP '\\'
#else
#define FILESEP '/'
#endif

/* Get the current time  */
void dw_gettime(struct timespec *);

/* Create dir if it does not exist.
 * Returns 0 if the dir already existed or could be created
 * returns non-zeros if the dir can't be created
 */
int ensuredir(const char * dir);

/* Check if directory exist
 * returns 1 if it exist
 * */

int isdir(const char * dir);

/* Return a suggestion for how many threads to use
 * ideally this should be the same as the number of cores
 */
int dw_get_threads(void);

/* Return 1 if the file exist, else 0
 * also returns 0 if fname is a folder (some consider that to be a file as well)
 */
int dw_isfile(const char * fname);

/* Difference between two timepoints */
float timespec_diff(struct timespec* end, struct timespec * start);

/* On Linux:
 * Reads /proc/self/status and assumes that the unit will be kb
 *
 * physical:
 * VmHWM = Peak Resident Set Size (peak RSS) i.e. how much physical memory
 * virtual:
 * VmPeak = peak Virtual Memory Size (peak VMZ) i.e. how much that is
 * allocated, although possibly never used.
 *
 * On MACOS:
 * physical = r_usage.ru_maxrss
 * virtual: set to 0
 *
 * On Windows: TODO
 * physical: 0
 * virtual: 0
 * Returns non 0
 *
 * Returns 0 on success.
 */

int get_peakMemoryKB(size_t * physical, size_t * virtual);

/* Print a line about peak memory usage to a file */
void fprint_peak_memory(FILE * fid);

/* Read the scaling of file from the .log.txt file if exists */
float dw_read_scaling(const char * file);

/* Linux : returns dirname() of path.
 * path is not altered and the returned strings has to be freed.
 * WINDOWS: returns drive and path
 * Example: "C:\home\dir\file.ext"" -> "C:\home\dir\"
 */
char * dw_dirname(const char * path);

/*
 * Linux: basename, compare to dw_dirname
* Windows: file name without extension.
* Should only be called for files, not for folders */
char * dw_basename(const char * path);

/* POSIX getcwd or _getcwd on windows */
char * dw_getcwd(char * buf, size_t size);

/** Add a prefix to a file
 * Examples:
 * ("file.tif", "dw") -> "dw_file.tif"
 * ("/dir/file.tif", "dw" -> "/dir/dw_file.tif"
 **/
char * dw_prefix_file(const char * file, const char * prefix);

#ifdef WINDOWS
int getline(char **lineptr, size_t *n, FILE *stream);
#endif


float abbe_res_xy(float lambda, float NA);
float abbe_res_z(float lambda, float NA);

int64_t
float_arg_max(const float * v, size_t N);
