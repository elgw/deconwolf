#include <assert.h>
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "dw_png.h"

typedef uint32_t u32;
typedef uint64_t u64;
typedef int64_t i64;
typedef uint8_t u8;


int
rgb_to_png(const u8 * pixels,
           u32 width, u32 height,
           const char * filename)
{
    assert(filename != NULL);
    assert(pixels != NULL);
    assert(width*height > 0);

    FILE * fid = fopen(filename, "wb");
    if(fid == NULL)
    {
        fprintf(stderr, "Unable to open %s for writing\n", filename);
        return EXIT_FAILURE;
    }

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                              NULL, NULL, NULL);
    if(png == NULL)
    {
        fprintf(stderr, "png_create_write_struct failed\n");
        fclose(fid);
        return EXIT_FAILURE;
    }

    png_infop info = png_create_info_struct(png);
    if(info == NULL)
    {
        fprintf(stderr, "png_create_info_struct failed\n");
        fclose(fid);
        return EXIT_FAILURE;
    }

    png_init_io(png, fid);

    png_set_IHDR(
        png,
        info,
        width, height,
        8,
        PNG_COLOR_TYPE_RGB,
        PNG_INTERLACE_NONE,
        PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT
        );

    //printf("Writing with width=%u, height=%u\n", width, height);
    const u8 ** row_pointers = calloc(height, sizeof(u8*));
    assert(row_pointers != NULL);
    for(u64 kk = 0; kk < height; kk++)
    {
        row_pointers[kk] = &pixels[3*kk*width];
    }

    png_write_info(png, info);
    png_write_image(png, (u8 **) row_pointers);
    png_write_end(png, NULL);
    png_destroy_write_struct(&png, &info);
    fclose(fid);
    free(row_pointers);

    return EXIT_SUCCESS;
}

/* Throw away the A component from an RGBA image */
static void
rgba_to_rgb(uint8_t * buffer, u32 M, u32 N)
{
    for(size_t kk = 0; kk < M*N; kk++)
    {
        for(int ii = 0; ii < 3; ii++)
        {
            buffer[3*kk+ii] = buffer[4*kk+ii];
        }
    }
    return;
}

u8 *
rgb_from_png(const char * filename,
             u32 * width, u32 * height)
{
    png_image image;
    memset(&image, 0, sizeof(image));
    image.version = PNG_IMAGE_VERSION;
    u8 * buffer = NULL;
    if(png_image_begin_read_from_file(&image, filename))
    {
        // The png library can perform conversions automatically
        // for us. However that triggers a lot of "Use of uninitialised value"
        // warnings in valgrind. Since those warnings clog up the terminal
        // I choose to only support RGB and RGBA with manual conversion.
        // To be able to read any format, just uncomment the next line:

        // image.format = PNG_FORMAT_RGB;

        if(image.warning_or_error)
        {
            printf("png error: %s\n", image.message);
        }

        buffer = calloc(PNG_IMAGE_SIZE(image), 1);

        if(buffer == NULL)
        {
            return NULL;
        }

        if(png_image_finish_read(&image, NULL, buffer, 0, NULL))
        {
            *height= image.height;
            *width = image.width;
        } else {
            free(buffer);
            return NULL;
        }
    }

    if(image.warning_or_error)
    {
        printf("png error: %s\n", image.message);
        free(buffer);
        return NULL;
    }
    assert(buffer != NULL);

    int known_format = 1;
    switch(image.format)
    {
    case PNG_FORMAT_RGB:
        known_format = 1;
        break;
    case PNG_FORMAT_RGBA:
        rgba_to_rgb(buffer, image.width, image.height);
        known_format = 1;
        break;
    default:
        printf("Warning: Unknown image format\n");
    }

    png_image_free(&image);
    if(known_format)
    {
        return buffer;
    } else {
        free(buffer);
        return NULL;
    }
}

int dw_png_ut(int argc, char ** argv)
{
    if(argc < 3)
    {
        printf("Usage:\n");
        printf("%s input.png output.png\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    u32 width, height;
    u8 * pixels = rgb_from_png(argv[1], &height, &width);
    if(pixels == NULL)
    {
        printf("Unable to read %s\n", argv[1]);
        exit(EXIT_FAILURE);
    }
    printf("Input dimensions: %u (height) x %u (width)\n",
           height, width);
    rgb_to_png(pixels, height, width, argv[2]);
    free(pixels);
    return EXIT_SUCCESS;
}
