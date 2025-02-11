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

/* Read and write tiff files to/from single precision floats
 * fim_tiff is not thread safe
 * You need to initialize with a call to
 * fim_tiff_init()
 * and should probably also redirect the output by
 * fim_tiff_set_log(FILE *)
 *
 */


#include <tiffio.h>


#include "ftab.h"
#include "dw_version.h"


#define XTAG_IJIJUNKNOWN 50838
#define XTAG_IJIJINFO 50839


/* Tiff tags -- for simple transfer from one image to another */
typedef struct{
    float xresolution;
    float yresolution;
    float zresolution;
    char * imagedescription;
    char * software;
    uint16_t resolutionunit;
    char * IJIJinfo; // Tag 50839 contains a string, used by Imagej.
    uint32_t nIJIJinfo;
    // Image size
    int M;
    int N;
    int P;
} ttags;

#include "fim.h"

/** ttags new with everything set to defaults */
ttags * ttags_new();

/** @brief Parse metadata from a open tif file */
void ttags_get(TIFF *, ttags *);

/** @brief print tags to a file */
void ttags_show(FILE *, ttags *);

/** @brief set tags to open tiff file*/
void ttags_set(TIFF *, const ttags *);

/** @brief Set software tag to S */
void ttags_set_software(ttags * ,
                        const char * S);

/** @brief set image size (in pixels) to tags */
void ttags_set_imagesize(ttags *, int M, int N, int P);

/** @brief set pixel size to tags
 * Note that the resolution unit is not set
*/
void ttags_set_pixelsize(ttags *, double xres, double yres, double zres);

/** @brief Free all data in a ttag* and set it to NULL */
void ttags_free(ttags **);

/* Initialization, sets the output file to stdout */
void fim_tiff_init(void);

/* Redirect all output here, both messages from tif_tiff as well as
 * warnings and errors from libtiff */

void fim_tiff_set_log(FILE * fp);

/* Write to disk, if scaling <= 0 : automatic scaling will be used. Else the provided value. */
int fim_tiff_write_opt(const char * fName, const float * V,
                       const ttags * T,
                       int64_t N, int64_t M, int64_t P, float scaling);


/* Scale between 0 and 2^16-1 and write data */
int fim_tiff_write(const char * fName, const float * V,
                   ttags * T,
                   int64_t M, int64_t N, int64_t P);

/* Don't scale data */
int fim_tiff_write_noscale(const char * fName, const float * V,
                           ttags * T,
                           int64_t N, int64_t M, int64_t P);


int fim_tiff_write_float(const char * fName, const float * V,
                         const ttags * T,
                         int64_t M, int64_t N, int64_t P);

int fim_tiff_write_zeros(const char * fName, int64_t M, int64_t N, int64_t P);

/** @brief Write raw floating point data to a 16-bit tif file.
 *
 * This performs a streaming write, never keeping all the data in memory.
 * @param output_file_name name of tiff file to write to
 * @param M, N, P the size of the image
 * @param raw_data_file_name file containing raw float data
 * @param meta_tiff_file Specify a tif file to copy metadata from. Can be NULL
 * @returns EXIT_SUCCESS or EXIT_FAILURE
 */

int
fim_tiff_from_raw(const char * output_tif_file_name,
                  int64_t M, int64_t N, int64_t P,
                  const char * raw_data_file_name,
                  const char * meta_tiff_file);

/** @brief Convert a tiff image to raw float image
 * Note that the image size is not encoded
*/
int
fim_tiff_to_raw(const char *tif_file_name,
                const char * output_file_name);

/* @brief Read a 3D tif stack as a float array
 * @param fName file name
 * @param verbosity how verbose the function should be
 * @param[out] M0, N0, P0, the image size in pixels
 * @param[out] T tiff tags will be written to T
 * @return The returned image is allocate with fim_malloc
 */
float * fim_tiff_read(const char * fName,
                      ttags * T,
                      int64_t * M0, int64_t * N0, int64_t * P0,
                      int verbosity);



/** @breif Read a sub region of a 3D stack as float array
 * set sub to 1
 * reads sM:sM+wM-1, sN:sN+wN-1, sP:sP+wP-1
 * @return The returned image is allocate with fim_malloc
 */
float * fim_tiff_read_sub(const char * fName,
                          ttags *,
                          int64_t * M0, int64_t * N0, int64_t * P0,
                          int verbosity,
                          int sub,
                          int64_t sM, int64_t sN, int64_t sP, // start
                          int64_t wM, int64_t wN, int64_t wP); // width

/** @brief Run self-tests
 *
*/
void fim_tiff_ut();

/** @brief Get image dimensions from tif file
 * @return 0 upon success.
 * @param fname Name of tiff file
 * @param[out] M, N, P the image size
 */

int fim_tiff_get_size(const char * fname,
                      int64_t * M, int64_t * N, int64_t * P);

/** @brief Max projection from input to output file
 * Not loading the full images into memory.
 * The output file will have the same sample format as the input image.
 */
int fim_tiff_maxproj(const char * in, const char * out);

/* Max projection from input to output file
 *  This version produces XY XZ and YZ */
int fim_tiff_maxproj_XYZ(const char * in, const char * out);


/** @brief Extract a single slice from input to output file */
int fim_tiff_extract_slice(const char *in, const char *out, int slice);
