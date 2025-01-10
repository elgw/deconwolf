#pragma once

#include <stdint.h>

/** Write an array of RGBA values to a png file. */

int
rgb_to_png(const uint8_t * Pixels,
           uint32_t width, uint32_t height,
           const char * filename);

/** Read a png image
 * Always into RGB format. Returns NULL on failure.
 * Can only read RGB and RGBA images, for other formats,
 * there is just one line to change in the code.
 */
uint8_t *
rgb_from_png(const char * filename,
             uint32_t * width, uint32_t * height);

int dw_png_ut(int argc, char ** argv);
