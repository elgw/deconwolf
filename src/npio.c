#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#ifdef WIN32
#define _CRT_SECURE_CPP_OVERLOAD_STANDARD_NAMES
#define _CRT_SECURE_NO_WARNINGS
#include <sys/types.h>
#endif

#include "npio.h"
#include "npio_config.h"

typedef double f64;
typedef float f32;
typedef int8_t i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;
typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;


#ifdef _WIN32
static char * strndup(const char * S, size_t n)
{
    char * Y = calloc(n+1, 1);
    if(Y == NULL)
    {
        return NULL;
    }
    Y[n] = '\0';
    for(size_t kk = 0; kk < n; kk++)
    {
        Y[kk] = S[kk];
        if(Y[kk] == '\0')
        {
            break;
        }
    }
    return Y;
}
#endif

static int fseek2(FILE *fid, int64_t offset, int origin)
{
    int ret = 0;
#ifdef _WIN32
    ret =  _fseeki64(fid, offset, origin);
#else
    ret = fseek(fid, offset, origin);
#endif
    if(ret)
    {
        perror("dw_fseek error:");
    }
    return ret;
}


/* FORWARD DECLARATIONS FOR DICTIONARY PARSER */

/* A dictionary-string parser,
   parses keywords and values
   and point to the corresponding positions in the dictionary string.
*/

typedef struct {
    size_t pos;
    int toknext;
} dp_t;

/* Token representation, start and end coordinates of the dictionary string */
typedef struct {
    int start; /* start position in dict string */
    int end; /* end position in dict string */
} dptok_t;

/** @brief Parse a simple dictionary in the string dict,
    since dict might not be \0 terminated it also needs to know the
    length of the string.
    The tokens should be pre-allocated to a fixed size supplied by tok_len.
*/
static int
dp_parse(dp_t dp, const char * dict,
         const size_t dict_len, dptok_t * tok, const int tok_len);

/** @brief Comparison of token string to the string s
 */
static int dp_eq(const char *dict, const dptok_t *tok, const char *s);

/* END OF FORWARD DECLARATIONS */
void npio_free(npio_t * np)
{
    free(np->filename);
    free(np->descr);
    free(np->shape_str);
    free(np->shape);
    free(np->data);
    free(np);
    return;
}

static int npd_parse_descr(npio_t * npd)
{
    if(strlen(npd->descr) != 5)
    {
        fprintf(stderr, "npio: unsupported descriptor: >>%s<<\n", npd->descr);
        fprintf(stderr, "      expected a non-nested descriptor\n");
        npd->np_type = NPIO_NOSUPPORT;
        npd->np_byte_order = '?';
        return EXIT_FAILURE;
    }
    char * descr = npd->descr;
    npd->np_byte_order = descr[1];
    npd->np_type = descr[2];
    // ok since the dictionary is at most 2^16-1 B
    npd->np_bytes = (int) strtol(descr+3, NULL, 10);
    if(npd->np_bytes < 1 || npd->np_bytes > 16)
    {
        fprintf(stderr,
                "npd_parse_descr: "
                "Could not figure out the size of the data type\n");
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

static npio_dtype
descr_to_dtype(const char * _descr)
{
    const char * descr = _descr;
    /* Sometimes the descriptor string is enclosed with single quotes
     */
    if(_descr[0] == '\'')
    {
        descr = _descr + 1;
    }

    switch(descr[1])
    {
    case 'f':
        switch(descr[2])
        {
        case '4':
            return NPIO_F32;
            break;
        case '8':
            return NPIO_F64;
            break;
        }
        break;
    case 'i':
        switch(descr[2])
        {
        case '1':
            return NPIO_I8;
            break;
        case '2':
            return NPIO_I16;
            break;
        case '4':
            return NPIO_I32;
            break;
        case '8':
            return NPIO_I64;
            break;
        }
        break;
    case 'u':
        switch(descr[2])
        {
        case '1':
            return NPIO_U8;
            break;
        case '2':
            return NPIO_U16;
            break;
        case '4':
            return NPIO_U32;
            break;
        case '8':
            return NPIO_U64;
            break;
        }
        break;
    }
    return NPIO_NOSUPPORT;
}

static const char * npio_type_to_descr(npio_dtype type)
{
    switch(type)
    {
    case NPIO_F32:
        return "<f4";
    case NPIO_F64:
        return "<f8";
    case NPIO_I8:
        return "<i1";
    case NPIO_I16:
        return "<i2";
    case NPIO_I32:
        return "<i4";
    case NPIO_I64:
        return "<i8";
    case NPIO_U8:
        return "<u1";
    case NPIO_U16:
        return "<u2";
    case NPIO_U32:
        return "<u4";
    case NPIO_U64:
        return "<u8";
    case NPIO_NOSUPPORT:
        return NULL;
    }
    return NULL;
}

static const char * npio_dtype_string(npio_dtype dtype)
{
    switch(dtype)
    {
    case NPIO_F32:
        return "NPIO_F32";
    case NPIO_F64:
        return "NPIO_F64";
    case NPIO_I8:
        return "NPIO_I8";
    case NPIO_I16:
        return "NPIO_I16";
    case NPIO_I32:
        return "NPIO_I32";
    case NPIO_I64:
        return "NPIO_I64";
    case NPIO_U8:
        return "NPIO_U8";
    case NPIO_U16:
        return "NPIO_U16";
    case NPIO_U32:
        return "NPIO_U32";
    case NPIO_U64:
        return "NPIO_U64";
    case NPIO_NOSUPPORT:
        return "NPIO_NOSUPPORT";
    }
    assert(0);
    return "ERROR";
}

static int npio_element_size(npio_dtype dtype)
{
    switch(dtype)
    {
    case NPIO_F32:
        return 4;
    case NPIO_F64:
        return 8;
    case NPIO_I8:
        return 1;
    case NPIO_I16:
        return 2;
    case NPIO_I32:
        return 4;
    case NPIO_I64:
        return 8;
    case NPIO_U8:
        return 1;
    case NPIO_U16:
        return 2;
    case NPIO_U32:
        return 4;
    case NPIO_U64:
        return 8;
    case NPIO_NOSUPPORT:
        return 0;
    }
    return 0;
}

static void print_dtype(FILE * fid, const npio_t * npd)
{
    char byte_order = npd->np_byte_order;
    char type = npd->np_type;;
    int nbytes = npd->np_bytes;
    switch(byte_order)
    {
    case '=':
        fprintf(fid, "native");
        break;
    case '<':
        fprintf(fid, "little endian");
        break;
    case '>':
        fprintf(fid, "big endian");
        break;
    case '|':
        fprintf(fid, "bit/or little endian");
        break;
    default:
        fprintf(fid, "?");
    }

    fprintf(fid, ", ");

    switch(type)
    {
    case 'i':
        fprintf(fid, "integer");
        break;
    case 'b':
        fprintf(fid, "boolean");
        break;
    case 'u':
        fprintf(fid, "unsigned integer");
        break;
    case 'f':
        fprintf(fid, "float");
        break;
    case 'c':
        fprintf(fid, "complex float");
        break;
    case 'm':
        fprintf(fid, "timedelta");
        break;
    case 'M':
        fprintf(fid, "datetime");
        break;
    case 'O':
        fprintf(fid, "object");
        break;
    case 'S':
        fprintf(fid, "string");
        break;
    case 'U':
        fprintf(fid, "unicode string");
        break;
    case 'V':
        fprintf(fid, "void");
        break;
    default:
        fprintf(fid, "?");
        break;
    }
    if(nbytes < 2)
    {
        fprintf(fid, ", %d byte", nbytes);
    } else {
        fprintf(fid, ", %d bytes", nbytes);
    }
}

// TODO: This prints too much, useful for debug though...
void npio_print(FILE * fid, const npio_t * np)
{
    if(fid == NULL)
    {
        fid = stdout;
    }
    fprintf(fid, "filename: %s\n", np->filename);
    if(np->dtype == NPIO_NOSUPPORT)
    {
        fprintf(fid, "descr: %s\n", np->descr);
        fprintf(fid, "The descriptor indicates a structured "
                "array, that is not supported\n");
    } else {
        fprintf(fid, "descr: %s (", np->descr);
        fprintf(fid, "npio_dtype: %s\n", npio_dtype_string(np->dtype));
        print_dtype(fid, np);
        fprintf(fid, ")\n");
        fprintf(fid, "np_byte_order: '%c'\n", np->np_byte_order);
        fprintf(fid, "np_type: '%c'\n", np->np_type);
        fprintf(fid, "fortran_order: %d\n", np->fortran_order);
        fprintf(fid, "ndim: %d\n", np->ndim);
    }
    fprintf(fid, "shape_str: '%s'\n", np->shape_str);
    fprintf(fid, "shape: [");
    for(int kk = 0 ; kk < np->ndim ; kk++)
    {
        if(kk>0)
        {
            fprintf(fid, ", ");
        }
        fprintf(fid, "%d", np->shape[kk]);
    }
    fprintf(fid, "]\n");
    fprintf(fid, "nel: %zu\n", np->nel);

    if(np->data == NULL)
    {
        fprintf(fid, "data is not loaded\n");
    } else {
        fprintf(fid, "size of data: %zu x %d = %zu B\n",
                np->nel, np->np_bytes, np->data_size);
    }
    fprintf(fid, "Data offset: %zu\n", np->data_offset);
}

static char *
gen_dictionary(int ndim, const int * shape, npio_dtype type_in)
{
    const i64 dict_alloc = ndim*12 + 128;
    i64 offset = 0;
    char * dict = calloc(dict_alloc, 1);
    assert(dict != NULL);
    offset += snprintf(dict+offset, dict_alloc-offset,
                       "{'descr': '%s', 'fortran_order': False, 'shape': ",
                       npio_type_to_descr(type_in));

    offset += snprintf(dict+offset, dict_alloc-offset,
                       "(");
    for(int kk = 0; kk+1<ndim; kk++)
    {
        offset += snprintf(dict+offset, dict_alloc-offset,
                           "%d, ", shape[kk]);
    }
    offset += snprintf(dict+offset, dict_alloc-offset,
                       "%d,)", shape[ndim-1]);
    offset += snprintf(dict+offset, dict_alloc-offset,
                       ", }");
    while( (10 + offset) % 64 != 63)
    {
        offset += snprintf(dict+offset, dict_alloc-offset,
                           "\x20");
    }
    offset += snprintf(dict+offset, dict_alloc-offset,
                       "\n");
    assert(dict_alloc > offset);
    dict[offset] = '\0';
    return dict;
}

static int parse_shape_string(npio_t * npd,
                              const char * sstring, const int len)
{
    //printf("To parse shape string: %.*s\n", len, sstring);
    char * str = strndup(sstring, len);

    if(str == NULL)
    {
        return EXIT_FAILURE;
    }
    //printf("<%s>\n", str);
    str[len-1] = ' ';
    str[0] = ' ';
    //printf("<%s>\n", str);

    char * saveptr = NULL;
    int ndim = 0;
    size_t nel = 1;
    int * shape = malloc(len*sizeof(int)); // more than big enough
    if(shape == NULL)
    {
        free(str);
        return EXIT_FAILURE;
    }

    char * str1 = NULL;
    for(str1 = str; ; str1 = NULL)
    {
#ifdef _WIN32
        char * tok = strtok_s(str1, ",", &saveptr);
#else
        char * tok = strtok_r(str1, ",", &saveptr);
#endif
        if(tok == NULL)
        {
            break;
        }
        char * endptr;
        int s = (int) strtol(tok, &endptr, 10);//atoi(tok);
        if(endptr == tok)
        {
            // no digits found
        } else {
            //printf("'%s' -> %d\n", tok, s);
            shape[ndim] = s;
            nel *= s;
            ndim++;
        }
    }
    shape[ndim] = 0; // 0-termination
    npd->ndim = ndim;
    npd->shape = shape;
    npd->nel = nel;
    free(str);
    return EXIT_SUCCESS;
}

static i64 get_file_size(FILE * fid)
{
#ifdef _WIN32
    struct _stat64 info;
    if(_fstat64(fileno(fid), &info))
    {
        return -1;
    }
    return info.st_size;
#else
    struct stat info;
    if(fstat(fileno(fid), &info))
    {
        return -1;
    }
    return info.st_size;
#endif
}

npio_t * npio_load_opts(const char * filename, int load_data)
{
    FILE * fid = fopen(filename, "rb");
    if(fid == NULL)
    {
        fprintf(stderr,
                "npio: failed to open %s\n", filename);
        return NULL;
    }

    i64 filesize = get_file_size(fid);

    if(filesize <= 0)
    {
        fprintf(stderr,
                "npio: could not fstat %s\n", filename);
        fclose(fid);
        return NULL;
    }

    // Check magic number and version
    char magic[] = "123456";

    size_t nread = fread(magic, 1, 6, fid);
    if(nread != 6)
    {
        fprintf(stderr, "npio: Could not read the magic number\n");
        fclose(fid);
        return NULL;
    }

    if(strncmp(magic, "\x93NUMPY", 6) != 0)
    {
        fprintf(stderr, "npio: Invalid magic number\n");
        fclose(fid);
        return NULL;
    }

    uint8_t version[2];
    nread = fread(&version, 1, 2, fid);
    if(nread != 2)
    {
        fclose(fid);
        fprintf(stderr, "npio: Couldn't read version\n");
        return NULL;
    }
    if(version[0] != 1 || version[1] != 0)
    {
        fprintf(stderr, "npio: Numpy file is v %d.%d, only tested for v 1.0\n",
                (int) version[0], (int) version[1]);
        fclose(fid);
        return NULL;
    }

    // Read the size of the dictionary
    uint16_t dsize = 0;
    nread = fread(&dsize, 1, 2, fid);
    if(nread != 2)
    {
        fprintf(stderr, "npio: could not read dictionary size\n");
        fclose(fid);
        return NULL;
    }
    //printf("Dictionary size: %u\n", dsize);

    // Read the dictionary
    char * dict = calloc(dsize+1, 1);
    assert(dict != NULL);

    nread = fread(dict, 1, dsize, fid);
    if(nread != dsize)
    {
        fprintf(stderr, "npio: could not read %d B dictionary\n", dsize);
        free(dict);
        fclose(fid);
        return NULL;
    }
    dict[dsize] = '\0'; // make a proper string
    //printf("Dictionary: %s\n", dict);

    npio_t * npd = calloc(1, sizeof(npio_t));
    assert(npd != NULL);
    npd->filename = strdup(filename);
    assert(npd->filename != NULL);
    npd->fortran_order = 1;

    // Parse the dictionary
    dp_t dp = {0}; // dictionary parser
    dptok_t t[10]; // expect at most 9 tokens
    // printf("DEBUG: >>%s<<\n", dict);
    int r = dp_parse(dp, dict, strlen(dict), t, 10);
    assert(r < 10);

    int found_fortran_order = 0;

    for(int kk = 0; kk+1<r; kk = kk + 2)
    {
        if(dp_eq(dict, t+kk, "'descr'") == 0)
        {

            npd->descr = strndup(dict+t[kk+1].start,
                                 t[kk+1].end-t[kk+1].start);

            npd->dtype = descr_to_dtype(npd->descr);

            if(npd_parse_descr(npd))
            {
                // goto fail;
                printf("there are warnings\n");
            }

        }
        else if(dp_eq(dict, t+kk, "'fortran_order'") == 0)
        {
            found_fortran_order = 1;
            int parsed = 0;
            if(dp_eq(dict, t+kk+1, "False") == 0)
            {
                npd->fortran_order = 0;
                parsed = 1;
            } else if(dp_eq(dict, t+kk+1, "True") == 0)
            {
                npd->fortran_order = 1;
                parsed = 1;
            }
            if(parsed == 0)
            {
                printf("Warning: Could not parse fortran_order\n");
                npd->fortran_order = 0; /* This is the default for numpy */
            }
        }
        else if(dp_eq(dict, t+kk, "'shape'") == 0)
        {
            int ret = parse_shape_string(npd,
                                         dict+t[kk+1].start,
                                         t[kk+1].end-t[kk+1].start);
            npd->shape_str = strndup(dict+t[kk+1].start,
                                     t[kk+1].end-t[kk+1].start);
            if(ret == EXIT_FAILURE)
            {
                fprintf(stderr, "npio: Could not parse the shape string\n");
                free(dict);
                dict = NULL;
                goto fail;
            }
        }
        /*
          printf(" %.*s = ",
          t[kk].end-t[kk].start, dict + t[kk].start);
          printf("%.*s\n",
          t[kk+1].end-t[kk+1].start, dict + t[kk+1].start);
        */
    }

    if(found_fortran_order == 0)
    {
        printf("Warning: fortran_order not specified (or not parsed)\n");
    }
    free(dict);

    if( npd->descr == NULL )
    {
        fprintf(stderr, "Failed to read the data type description\n");
        goto fail;
    }


    if(npd->dtype == NPIO_NOSUPPORT)
    {
        if(load_data == 1)
        {
            fprintf(stderr, "npio: Will not be able to load data\n");
        }
        load_data = 0;
    }




    // Forward to the data
    long pos = ftell(fid);
    while(((size_t) pos % 64) != 0)
    {
        pos++;
    }
    r = fseek2(fid, pos, SEEK_SET);
    if(r != EXIT_SUCCESS)
    {
        fprintf(stderr, "npio: fseek failed on line %d\n", __LINE__);
        goto fail;
    }

    npd->data_offset = ftell(fid);

    if(load_data == 0)
    {
        npd->data_size = 0;
        npd->data = NULL;
        goto post_data;
    }

    /// Read the data
    size_t nBytes = (u64) npd->nel* (u64) npd->np_bytes;
    // Don't be tricked to allocate more memory than the actual
    // file size
    if(nBytes > (size_t) filesize)
    {
        fprintf(stderr, "npio: nBytes=%zu, file size=%zu\n", nBytes, filesize);
        fprintf(stderr, "npio: corrupt npy file?\n");
        goto fail;
    }
    //printf("To read %zu elements in %zu B\n", npd->nel, nBytes);

    uint8_t * data = calloc(nBytes, 1);
    if(data == NULL)
    {
        fprintf(stderr, "npio failed: Could not allocate %zu bytes\n", nBytes);
        goto fail1;
    }
    nread = fread(data, npd->np_bytes, npd->nel, fid);
    if(nread != npd->nel)
    {
        fprintf(stderr, "npio failed: Could not read %zu elements of size %d\n",
                npd->nel, npd->np_bytes);
        goto fail1;
    }
    npd->data_size = npd->np_bytes*npd->nel;
    npd->data = data;

post_data:
    fclose(fid);
    return npd;

fail1:
    free(data);
fail:
    npio_free(npd);
    npd = NULL;
    fclose(fid);
    return NULL;
}

npio_t * npio_load(const char * filename)
{
    return npio_load_opts(filename, 1);
}

npio_t * npio_load_metadata(const char * filename)
{
    return npio_load_opts(filename, 0);
}


int64_t
npio_write_FILE(FILE * fid,
                const int ndim,
                const int * shape,
                const void * data,
                npio_dtype type_in, npio_dtype type_out)
{
    if(type_in != type_out)
    {
        fprintf(stdout, "Input and output format combination not supported\n");
        return -1;
    }
    if(type_in == NPIO_NOSUPPORT)
    {
        fprintf(stderr, "Unknown/unsupported input data type\n");
        fprintf(stderr, "%s %d\n", __FILE__, __LINE__);
        return -1;
    }
    if(type_out == NPIO_NOSUPPORT)
    {
        fprintf(stderr, "Unknown/unsupported output data type\n");
        fprintf(stderr, "%s %d\n", __FILE__, __LINE__);
        return -1;
    }
    size_t nwritten = 0;

    size_t nelements = 1;
    for(int kk = 0; kk<ndim; kk++)
    {
        nelements *= shape[kk];
    }

    char d[] = "\x93NUMPY";
    size_t nw = fwrite(d, 1, 6, fid);
    if(nw != 6)
    {
        printf("Failed to write %s %d\n", __FILE__, __LINE__);
        goto fail;
    }
    nwritten += 6*1;

    /// Write major and minor version
    d[0] = 1;
    d[1] = 0;
    nw = fwrite(d, 1, 2, fid);
    if(nw != 2)
    {
        fprintf(stdout, "Failed to write version\n");
        goto fail;
    }
    nwritten += 1*2;

    // Generate the dictionary
    // write the size of the dictionary
    // write the dictionary
    char * dictionary = gen_dictionary(ndim, shape, type_in);
    u16 dlen = strlen(dictionary);
    f64 hl = fwrite(&dlen, 2, 1, fid);
    assert(hl == 1);
    i64 inw = fwrite(dictionary, 1, dlen, fid);
    assert(inw == dlen);
    free(dictionary);

    nwritten += inw + hl;

    /// write the data
    //printf("Will write %zu elements\n", nelements);
    //printf("First element: %f\n", data[0]);
    int element_size = npio_element_size(type_in);
    assert(element_size > 0);
    if(data != NULL)
    {
        nw = fwrite(data, element_size, nelements, fid);
        if(nw != (size_t) nelements)
        {
            fprintf(stderr, "Failed to write data %s %d\n", __FILE__, __LINE__);
            printf("Expeced nw = %ld == nelements %lu\n", nw, nelements);
            printf("element_size: %u\n", element_size);
            goto fail;
        }
        nwritten += element_size*nelements;
    }
    return nwritten;

fail:
    return -1;
}


i64 npio_write(const char * fname,
               const int ndim,
               const int * shape,
               const void * data,
               npio_dtype dtin, npio_dtype dtout)
{
    FILE * fid = fopen(fname, "wb");
    if(fid == NULL)
    {
        printf("Unable to open %s for writing\n", fname);
        return -1;
    }
    i64 nwritten = npio_write_FILE(fid, ndim, shape, data, dtin, dtout);
    fclose(fid);
    return nwritten;
}

void *
npio_write_mem(const int ndim,
               const int * shape,
               const void * data,
               npio_dtype type_in,
               npio_dtype type_out,
               int64_t * mem_size)
{
    if(type_in != type_out)
    {
        fprintf(stdout, "Input and output format combination not supported\n");
        return NULL;
    }
    if(type_in == NPIO_NOSUPPORT)
    {
        fprintf(stderr, "Unknown/unsupported input data type\n");
        fprintf(stderr, "%s %d\n", __FILE__, __LINE__);
        return NULL;
    }
    if(type_out == NPIO_NOSUPPORT)
    {
        fprintf(stderr, "Unknown/unsupported output data type\n");
        fprintf(stderr, "%s %d\n", __FILE__, __LINE__);
        return NULL;
    }

    size_t nelements = 1;
    for(int kk = 0; kk<ndim; kk++)
    {
        nelements *= shape[kk];
    }
    int element_size = npio_element_size(type_in);

    char * dictionary = gen_dictionary(ndim, shape, type_in);
    uint8_t * buff = calloc(10 + strlen(dictionary) + nelements*element_size, 1);
    size_t dictionary_size = strlen(dictionary);
    char header[6] = "\x93NUMPY";
    memcpy(buff, header, 6);
    buff[6] = 1;
    buff[7] = 0;
    u16 * _dsize = (u16*) (buff+8);
    _dsize[0] = dictionary_size;
    assert( (10 + dictionary_size) % 64 == 0);
    /* Write the dictionary */
    memcpy(buff + 10, dictionary, dictionary_size);
    free(dictionary);

    /* Copy the data */
    memcpy(buff+(10+dictionary_size), data, element_size*nelements);

    /* Set the final size */
    *mem_size = 10+dictionary_size + element_size*nelements;
    return buff;
}

/* DICTIONARY PARSER

   Supported descriptors look like this:
   {'descr': '<f8', 'fortran_order': False, 'shape': (2, 3, 4), }
   Valid but not supported descriptor:
   {'descr': [('name', '<U10'), ('age', '<i4'), ('weight', '<f4')], 'fortran_order': False, 'shape': (2,), }

*/

static int dp_eq(const char *dict, const dptok_t *tok, const char *s)
{
    assert(dict != NULL);
    assert(tok != NULL);
    assert(s != NULL);

    int len = (int) strnlen(s, tok->end - tok->start);

    if(len == 0)
    {
        return -1;
    }

    if (len != tok->end - tok->start)
    {
        return -1;
    }

    if(strncmp(dict + tok->start, s, tok->end - tok->start) == 0)
    {
        return 0;
    }
    return -1;
}

static int
dp_parse(dp_t dp,
         const char * dict,
         const size_t dict_len,
         dptok_t * tok,
         const int tok_len)
{
    dp.pos = 0; /* position in dict string */
    dp.toknext = 0; /* next token to write to */

    /* Find first '{' while ignoring white spaces */
    while(dict[dp.pos] != '{' && dp.pos < dict_len)
    {
        if( dp.pos == dict_len)
        {
            return 0;
        }
        dp.pos++;
    }

    /* States:
     * 0: looking for key to start '
     * 1: looking for key to end ''
     * 2: looking for ':'
     * 3: looking for value to start nonwhite
     * 4: looking for value to end ,
     */
    int state = 0;
    int nest = 0; /* Keep track of parentheses */
    while( dp.pos < dict_len)
    {
        char c = dict[dp.pos];
        if(state == 0) /* looking for key to start */
        {
            if(c == '\'')
            {
                state = 1;
                tok[dp.toknext].start = dp.pos;
                goto nextchar;
            }
        }
        if(state == 1) /* looking for key to end */
        {
            if(c == '\'')
            {
                state = 2;
                tok[dp.toknext].end = dp.pos+1;
                dp.toknext++;
                goto nextchar;
            }
        }
        if(state == 2)
        {
            if(c == ':') /* looking for key:value separator */
            {
                state = 3;
                goto nextchar;
            }
        }
        if(state == 3) /* looking for start of value */
        {
            nest = 0;
            if(c != ' ')
            {
                switch(c)
                {
                case '(':
                    nest++;
                    break;
                case '[':
                    nest++;
                    break;
                case '{':
                    nest++;
                    break;
                }

                state = 4;
                tok[dp.toknext].start = dp.pos;
                goto nextchar;
            }
        }
        if(state == 4) /* looking for end of value */
        {
            switch(c)
            {
            case '(':
                nest++;
                break;
            case '[':
                nest++;
                break;
            case '{':
                nest++;
                break;
            case ')':
                nest--;
                break;
            case ']':
                nest--;
                break;
            case '}':
                nest--;
                break;
            }

            if( nest == 0 && c == ',')
            {
                state = 0;
                tok[dp.toknext].end = dp.pos;
                dp.toknext++;
                goto nextchar;
            }
        }
    nextchar: ;
        dp.pos++;
        if(dp.toknext == tok_len)
        {
            // Could warn here...
            return dp.toknext;
        }
    }
    //printf("\n");
    return dp.toknext;
}

const char * npio_version(void)
{
    return NPIO_VERSION;
}

int npio_version_major(void)
{
    return atoi(NPIO_VERSION_MAJOR);
}

int npio_version_minor(void)
{
    return atoi(NPIO_VERSION_MINOR);
}

int npio_version_patch(void)
{
    return atoi(NPIO_VERSION_PATCH);
}
