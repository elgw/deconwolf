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
rgb_to_png(const uint8_t * pixels,
            u32 height, u32 width,
            const char * filename)
{
    assert(filename != NULL);
    assert(pixels != NULL);
    assert(width*height > 0);

    FILE * fid = fopen(filename, "w");
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

    printf("Writing with width=%u, height=%u\n", width, height);
    const uint8_t ** row_pointers = calloc(height, sizeof(uint8_t*));
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


uint8_t *
rgb_from_png(const char * filename,
              uint32_t * height, uint32_t * width)
{
    png_image image;
    memset(&image, 0, sizeof(image));
    image.version = PNG_IMAGE_VERSION;
    u8 * buffer = NULL;
    if(png_image_begin_read_from_file(&image, filename))
    {
        image.format = PNG_FORMAT_RGB;
        buffer = malloc(PNG_IMAGE_SIZE(image));

        if(buffer == NULL)
        {
            fprintf(stderr, "Unable to allocate memory for the image buffer\n");
            exit(EXIT_FAILURE);
        }

        if(png_image_finish_read(&image, NULL, buffer, 0, NULL))
        {
            *height= image.height;
            *width = image.width;
        } else {
            free(buffer);
            buffer = NULL;
        }
    }

    png_image_free(&image);
    return buffer;
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
