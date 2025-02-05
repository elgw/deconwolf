#pragma once

/* A library to read/write numpy .npy-files
 * see npio_cli.c for example usage.
 *
 * web: https://www.github.com/elgw/npio
 */

#include <stdint.h>
#include <stdio.h>

#define NPIO_VERSION_MAJOR "0"
#define NPIO_VERSION_MINOR "0"
#define NPIO_VERSION_PATCH "6"
#define NPIO_version NPIO_VERSION_MAJOR "."     \
    NPIO_VERSION_MINOR "."                      \
    NPIO_VERSION_PATCH

/* This is what is returned from npio_load */
typedef struct{
    /* Read from the file */
    char * descr; /* data type description */
    char np_byte_order;
    char np_type;
    int np_bytes;
    int fortran_order;
    int ndim; // number of dimensions
    int * shape;
    char * shape_str;
    size_t nel; // number of elements
    void * data; // raw pointer to the data
    size_t data_size; // total size of data, in bytes

    /* Not from the npy file */
    char * filename;
} npio_t;

/** Read a .npy file
 *
 * Returns NULL on failure and might print a message to stdout
 * the returned struct should be freed with npio_free
 *
 * The pointer to the data can be stolen:
 *    double * my_data = (double*) np->data;
 *    np->data = NULL;
 *    npio_free(np); // will ignore np->data
 */
npio_t * npio_load(const char * filename);


/** Save a double array to an npy file */
int npio_save_f64(const char * filename,
                  const int ndim, const int * shape,
                  const double * data);

int npio_save_f32(const char * filename,
                  const int ndim, const int * shape,
                  const float * data);

int npio_save_i8(const char * filename,
                 const int ndim, const int * shape,
                 const int8_t * data);

int npio_save_i16(const char * filename,
                  const int ndim, const int * shape,
                  const int16_t * data);

int npio_save_i32(const char * filename,
                  const int ndim, const int * shape,
                  const int32_t * data);

int npio_save_i64(const char * filename,
                  const int ndim, const int * shape,
                  const int64_t * data);

int npio_save_u8(const char * filename,
                 const int ndim, const int * shape,
                 const uint8_t * data);

int npio_save_u16(const char * filename,
                  const int ndim, const int * shape,
                  const uint16_t * data);

int npio_save_u32(const char * filename,
                  const int ndim, const int * shape,
                  const uint32_t * data);

int npio_save_u64(const char * filename,
                  const int ndim, const int * shape,
                  const uint64_t * data);



/** @brief Print some info about the npio_t object
 *
 */
void npio_print(FILE *, const npio_t * np);

/**  Free an npio_t object
 */
void npio_free(npio_t * np);
