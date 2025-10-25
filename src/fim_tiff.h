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

/* PROVIDES
 * - Read tiff files as floats
 * - Writes tiff files as uint16 or float
 *
 * TODO
  * - Remove all dependencies
 * - Rename to vtif_* (volumetric tif)
 * - Correct imagej metadata when writing 2D images (and the metadata says 3D)
 *   I suggest parsing the IJ metadata and recrating it in a meaningful way on writing.
 * - Read individual image plane(s).
 * - Read and write back OME-XML metadata ?
 *
 * BUGS:
 * - When writing max projections, the metadata is
 *   transferred from the 3D images. If the source image had an
 *   imagedescription tag it will be wrong for the 2D output.
 *
 * The image description can look like this:
 * ImageJ=1.52r
 * images=71
 * slices=71
 * unit=nm
 * spacing=200.0
 * loop=false
 *
 * In that case the following lines are wrong/irrelevant:
 *
 * images=71
 * slices=71
 * loop=false
 *
 * and should be stripped. At the same time we would like to keep the 'unit=...'
 * and possibly all the other information as well.
 *
 * NOTES
 * _TIFFmalloc can safely be replaced by any other allocator
 * (it simply wraps the system malloc, possibly it did something else
 * a long time ago)
*/

#include <tiffio.h> // https://gitlab.com/libtiff/libtiff

typedef struct {
    FILE * log;
    int verbose;
    // would make sense to include, error, stdout, allocator, ...
} ftif_t;

/* Initialization */
ftif_t * fim_tiff_new(FILE * log, int verbose);

/* Destruction stuff */
void fim_tiff_destroy(ftif_t * ftif);


/** For storing tiff meta data */
typedef struct _ttags ttags;

/** ttags new with everything set to defaults */
ttags * ttags_new();

/** @brief Parse metadata from a open tif file */
void ttags_get(ftif_t *, TIFF *, ttags *);

/** @brief print tags to a file */
void ttags_show(ftif_t *, FILE *, ttags *);

/** @brief set tags to open tiff file*/
void ttags_set(ftif_t *, TIFF *, const ttags *);

/** @brief Set software tag to S */
void ttags_set_software(ftif_t *,
                        ttags * ,
                        const char * S);

/** @brief set image size (in pixels) to tags */
void ttags_set_imagesize(ftif_t *, ttags *, int M, int N, int P);

/** @brief set pixel size to tags
 * Note that the resolution unit is not set
 */
void ttags_set_pixelsize(ftif_t *, ttags *, double xres, double yres, double zres);

/** @brief Free all data in a ttag* and set it to NULL */
void ttags_free(ttags **);

/* TIFF input and output  */

/* Write to disk, if scaling <= 0 : automatic scaling will be used. Else the provided value. */
int fim_tiff_write_opt(ftif_t * ftif,
                       const char * fName, const float * V,
                       const ttags * T,
                       int64_t N, int64_t M, int64_t P, float scaling);


/* Scale between 0 and 2^16-1 and write data */
int fim_tiff_write(ftif_t * ftif,
                   const char * fName, const float * V,
                   ttags * T,
                   int64_t M, int64_t N, int64_t P);

/* Don't scale data */
int fim_tiff_write_noscale(ftif_t * ftif,
                           const char * fName, const float * V,
                           ttags * T,
                           int64_t N, int64_t M, int64_t P);


int fim_tiff_write_float(ftif_t * ftif,
                         const char * fName, const float * V,
                         const ttags * T,
                         int64_t M, int64_t N, int64_t P);

/** @brief Write a 16-bit tif file from f32 raw data
 *
 * This performs a streaming write, never keeping all the data in memory.
 * @param output_file_name name of tiff file to write to
 * @param M, N, P the size of the image
 * @param raw_data_file_name file containing raw float data
 * @param meta_tiff_file Specify a tif file to copy metadata from. Can be NULL
 * @param scaling Specify a scaling values for all pixels. If <=0 this will be set to use the full dynamic range of the image.
 * @returns EXIT_SUCCESS or EXIT_FAILURE
 */

int
fim_tiff_imwrite_u16_from_raw(ftif_t * ftif,
                              const char * output_tif_file_name,
                  int64_t M, int64_t N, int64_t P,
                  const char * raw_data_file_name,
                              const char * meta_tiff_file,
    float scaling);


/** @brief Write a f32 tif file from f32 raw data
 *
 * This performs a streaming write, never keeping all the data in memory.
 * @param output_file_name name of tiff file to write to
 * @param M, N, P the size of the image
 * @param raw_data_file_name file containing raw float data
 * @param meta_tiff_file Specify a tif file to copy metadata from. Can be NULL
 * @returns EXIT_SUCCESS or EXIT_FAILURE
 */

int
fim_tiff_imwrite_f32_from_raw(ftif_t * ftif,
    const char * fName, // Name of tiff file to be written
    int64_t M, int64_t N, int64_t P, // Image dimensions
    const char * rName,  // name of raw file
    const char * meta_tiff_file);

/** @brief Convert a tiff image to raw float image
 * Note that the image size is not encoded
 */
int
fim_tiff_to_raw_f32(ftif_t * ftif,
                    const char *tif_file_name,
                const char * output_file_name);

/* @brief Read a 3D tif stack as a float array
 * @param fName file name
 * @param[out] M0, N0, P0, the image size in pixels
 * @param[out] T tiff tags will be written to T. Ignored if NULL.
 * @return The returned image is allocate with fim_malloc
 */
float * fim_tiff_read(ftif_t * ftif,
                      const char * fName,
                      ttags * T,
                      int64_t * M0, int64_t * N0, int64_t * P0);



/** @breif Read a sub region of a 3D stack as float array
 * set sub to 1
 * reads sM:sM+wM-1, sN:sN+wN-1, sP:sP+wP-1
 * @return The returned image is allocate with fim_malloc
 */
float * fim_tiff_read_sub(ftif_t * ftif,
                          const char * fName,
                          ttags *,
                          int64_t * M0, int64_t * N0, int64_t * P0,
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
 *
 * Please use fim_tiff_get_info instead.
 */

int fim_tiff_get_size(ftif_t * ftif, const char * fname,
                      int64_t * M, int64_t * N, int64_t * P);

typedef struct {
    int64_t M; // Number of elements along non-strided dimension
    int64_t N;
    int64_t P;
    uint32_t BPS; // Bits per sample
} fim_tiff_info;

/** @brief Get some basic information about a tiff file without reading the image data.
 * @return 0 upon success.
 * @param fname Name of tiff file
 * @param[out] M, N, P the image size
 */

int fim_tiff_get_info(ftif_t * ftif,
                      const char * fname, fim_tiff_info * info);

/** @brief Max projection from input to output file
 * Not loading the full images into memory.
 * The output file will have the same sample format as the input image.
 *
 * return value: 0 on success
 */
int fim_tiff_maxproj(ftif_t * ftif,
                     const char * in, const char * out);

/* Max projection from input to output file
 *  This version produces XY XZ and YZ */
int fim_tiff_maxproj_XYZ(ftif_t * ftif,
                         const char * in, const char * out);


/** @brief Extract a single slice from input to output file */
int fim_tiff_extract_slice(ftif_t * ftif,
                           const char *in, const char *out, int slice);

/** @brief Read the max value of a file consisting of f32 (raw) */
float raw_file_single_max(ftif_t * ftif,
                          const char * rName, const size_t N);
