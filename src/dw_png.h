#pragma once

#include <stdint.h>

/** Write an array of RGBA values to a png file. */

int
rgba_to_png(const uint8_t * Pixels,
            uint32_t height, uint32_t width,
            const char * filename);

uint8_t *
rgba_from_png(const char * filename,
                        uint32_t * height, uint32_t * width);
