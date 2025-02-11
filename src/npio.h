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
#define NPIO_VERSION_PATCH "8"
#define NPIO_version NPIO_VERSION_MAJOR "."     \
    NPIO_VERSION_MINOR "."                      \
    NPIO_VERSION_PATCH

typedef enum
{ NPIO_F64,
  NPIO_F32,
  NPIO_U8, NPIO_U16, NPIO_U32, NPIO_U64,
  NPIO_I8, NPIO_I16, NPIO_I32, NPIO_I64,
  NPIO_NOSUPPORT
} npio_dtype;


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
    npio_dtype dtype;
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

/* Read the metadata but do not load the data */
npio_t * npio_load_metadata(const char * filename);

/* Write data to a file decriptor, such as retrieved from fopen or fmemopen
 * return the number of bytes written or -1 on failure
 */
int64_t
npio_write_FILE(FILE * fid,
                const int ndim,
                const int * shape,
                void * data,
                npio_dtype in, npio_dtype out);

int64_t
npio_write(const char * fname,
           const int ndim,
           const int * shape,
           void * data,
           npio_dtype in, npio_dtype out);


/** @brief Print some info about the npio_t object
 *
 */
void npio_print(FILE *, const npio_t * np);

/**  Free an npio_t object
 */
void npio_free(npio_t * np);
