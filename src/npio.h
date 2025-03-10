#pragma once

/* npio: a library to read/write numpy .npy-files
 * see npio_cli.c for example usage.
 *
 * Only intended for/tested on numeric arrays.
 *
 * web: https://www.github.com/elgw/npio
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdio.h>

    typedef enum {
        NPIO_F64,
        NPIO_F32,
        NPIO_U8, NPIO_U16, NPIO_U32, NPIO_U64,
        NPIO_I8, NPIO_I16, NPIO_I32, NPIO_I64,
        NPIO_NOSUPPORT
    } npio_dtype;


    /* This is what is returned from npio_load */
    typedef struct {
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
        size_t data_offset; // Where the data starts in the file
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
     * Note: If data == NULL the function will write the metadata only.
     */
    int64_t
    npio_write_FILE(FILE * fid,
                    const int ndim,
                    const int * shape,
                    const void * data,
                    npio_dtype in, npio_dtype out);

    /* Write to a file given by its name. Overwrites existing files by
     * default. See npio_write_FILE */

    int64_t
    npio_write(const char * fname,
               const int ndim,
               const int * shape,
               const void * data,
               npio_dtype in, npio_dtype out);

    /* Write an npy file to a memory buffer
     *
     * On success:
     * returns a memory buffer of mem_size bytes
     *
     * On failure:
     * returns NULL
     */
    void *
    npio_write_mem(const int ndim,
                   const int * shape,
                   const void * data,
                   npio_dtype in,
                   npio_dtype out,
                   int64_t * mem_size);

    /** Print some info about the npio_t object
     */
    void npio_print(FILE *, const npio_t * np);

    /**  Free an npio_t object
     */
    void npio_free(npio_t * np);

    /* npio_version returns the version number of the library in the
       format "$MAJOR.$MINOR.$PATCH"
    */

    const char * npio_version(void);

    int npio_version_major(void);
    int npio_version_minor(void);
    int npio_version_patch(void);

#ifdef __cplusplus
}
#endif
