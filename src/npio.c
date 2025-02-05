#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "npio.h"

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
static int dp_parse(dp_t dp, const char * dict,
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
        fprintf(stderr, "Invalid descriptor: %s ", npd->descr);
        fprintf(stderr, "Should be of length 3, not %zu\n",
                (ssize_t) strlen(npd->descr) - 2 );
        exit(EXIT_FAILURE);
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


void npio_print(FILE * fid, const npio_t * np)
{
    if(fid == NULL)
    {
        fid = stdout;
    }
    fprintf(fid, "filename: %s\n", np->filename);
    fprintf(fid, "descr: %s (", np->descr);
    print_dtype(fid, np);
    fprintf(fid, ")\n");
    fprintf(fid, "np_byte_order: '%c'\n", np->np_byte_order);

    fprintf(fid, "np_type: '%c'\n", np->np_type);

    fprintf(fid, "fortran_order: %d\n", np->fortran_order);
    fprintf(fid, "ndim: %d\n", np->ndim);
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
    fprintf(fid, "size of data: %zu x %d = %zu B\n",
            np->nel, np->np_bytes, np->data_size);
}


static int write_dictionary(FILE * fid, const int ndim, const int * shape,
                            const char * desc_str)
{

    if(ndim < 1)
    {
        return EXIT_FAILURE;
    }

    /// Create the shape string i.e. something like (2, 3,)
    char *shape_str;
    size_t size;
    FILE *stream = open_memstream(&shape_str, &size);
    if(stream == NULL)
    {
        return EXIT_FAILURE;
    }

    fprintf(stream, "(");
    for(int kk = 0; kk+1<ndim; kk++)
    {
        fprintf(stream, "%d, ", shape[kk]);
    }
    fprintf(stream, "%d,", shape[ndim-1]);
    fprintf(stream, ")");
    fclose(stream);
    //printf("Will write shape_str: %s\n", shape_str);

    /// Write the header
    // - header should end with '\n' and then be padded by '\x20'
    // so that
    // len(magic string) + 2 + len(length) + HEADER_LEN should be
    // divisible by 64
    // i.e. HEADER_LEN + 10 should be divisible by 64.

    char * dict = malloc(128+strlen(shape_str));
    if(dict == NULL)
    {
        return EXIT_FAILURE;
    }

    sprintf(dict,
            "{'descr': '%s', 'fortran_order': False, 'shape': %s, }",
            desc_str, shape_str);
    free(shape_str);

    // printf("dict to write: %s\n", dict);

    size_t _HEADER_LEN = strlen(dict);
    if(_HEADER_LEN > 2<<16)
    {
        goto fail;
    }
    uint16_t HEADER_LEN = (uint16_t) _HEADER_LEN;
    while( (10 + HEADER_LEN) % 64 != 63)
    {
        dict[HEADER_LEN++] = '\x20';
    }
    dict[HEADER_LEN++] = '\n';

    // write length of dictionary (including the padding)
    size_t nwritten = fwrite(&HEADER_LEN, 2, 1, fid);
    if(nwritten != 1)
    {
        goto fail;
    }

    // write the dictionary
    nwritten = fwrite(dict, HEADER_LEN, 1, fid);
    if(nwritten != 1)
    {
        goto fail;
    }

    //printf("nelements: %zu\n", nelements);
    free(dict);
    return EXIT_SUCCESS;

fail:
    free(dict);
    return EXIT_FAILURE;
}


static int npio_write_raw(const char * filename,
                           const int ndim, const int * shape,
                           const void * data,
                          size_t nelements, size_t element_size,
                          const char * desc_str)
{
    /// Write 6 byte magic header
    FILE * fid = fopen(filename, "w");
    if(fid == NULL)
    {
        fprintf(stderr, "npio_save_double: Unable to open %s for writing\n",
                filename);
        return EXIT_FAILURE;
    }
    char d[] = "\x93NUMPY";
    size_t nwritten = fwrite(d, 1, 6, fid);
    if(nwritten != 6)
    {
        goto fail;
    }

    /// Write major and minor version
    d[0] = 1;
    d[1] = 0;
    nwritten = fwrite(d, 1, 2, fid);
    if(nwritten != 2)
    {
        goto fail;
    }

    /// Write the part of the dictionary that describes the shape
    int status = write_dictionary(fid, ndim, shape, desc_str);
    if(status != EXIT_SUCCESS)
    {
        goto fail;
    }

    /// write the data
    //printf("Will write %zu elements\n", nelements);
    //printf("First element: %f\n", data[0]);
    nwritten = fwrite(data, element_size, nelements, fid);
    if(nwritten != nelements)
    {
        goto fail;
    }

    fclose(fid);
    return EXIT_SUCCESS;

fail:
    fclose(fid);
    return EXIT_FAILURE;
}


static int parse_shape_string(npio_t * npd,
                               const char * sstring, const int len)
{
    //printf("To parse shape string: %.*s\n", len, sstring);
    char * str = strndup(sstring, len);
    //printf("<%s>\n", str);
    str[len-1] = ' ';
    str[0] = ' ';
    //printf("<%s>\n", str);

    char * saveptr = NULL;
    char * str1 = str;
    int ndim = 0;
    size_t nel = 1;
    int * shape = malloc(len*sizeof(int)); // more than big enough
    for(str1 = str; ; str1 = NULL)
    {
        char * tok = strtok_r(str1, ",", &saveptr);
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

npio_t * npio_load(const char * filename)
{
    FILE * fid = fopen(filename, "r");
    if(fid == NULL)
    {
        fprintf(stderr,
                "npio: failed to open %s\n", filename);
        return NULL;
    }

    struct stat info;
    if(fstat(fileno(fid), &info) != 0)
    {
        fprintf(stderr,
                "npio: could not fstat %s\n", filename);
        fclose(fid);
        return NULL;
    }

    size_t filesize = info.st_size;

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
    npd->filename = strdup(filename);
    npd->fortran_order = 1;

    // Parse the dictionary
    dp_t dp = {0}; // dictionary parser
    dptok_t t[10]; // expect at most 9 tokens
    int r = dp_parse(dp, dict, strlen(dict), t, 10);
    assert(r < 10);

    for(int kk = 0; kk+1<r; kk = kk + 2)
    {
        if(dp_eq(dict, t+kk, "'descr'") == 0)
        {
            npd->descr = strndup(dict+t[kk+1].start,
                                 t[kk+1].end-t[kk+1].start);
            if(npd_parse_descr(npd))
            {
                goto fail;
            }

        }
        else if(dp_eq(dict, t+kk, "'fortran_order'") == 0)
        {
            if(dp_eq(dict, t+kk+1, "'False'"))
            {
                npd->fortran_order = 0;
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

    free(dict);

    if( npd->descr == NULL )
    {
        fprintf(stderr, "Failed to read the data type description\n");
        goto fail;
        }

    // Forward to the data
    long pos = ftell(fid);
    while(((size_t) pos % 64) != 0)
    {
        pos++;
    }
    r = fseek(fid, pos, SEEK_SET);
    if(r != EXIT_SUCCESS)
    {
        fprintf(stderr, "npio: fseek failed on line %d\n", __LINE__);
        goto fail;
    }

    /// Read the data
    size_t nBytes = (u64) npd->nel* (u64) npd->np_bytes;
    // Don't be tricked to allocate more memory than the actual
    // file size
    if(nBytes > filesize)
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
    fclose(fid);

    npd->data = data;
    return npd;

 fail1:
    free(data);
 fail:
    npio_free(npd);
    npd = NULL;
    fclose(fid);
    return NULL;
}

/* Save function signature */

#define NPIO_SAVE(x, y)                                  \
int npio_save_## x(const char * filename,                \
                     const int ndim, const int * shape,  \
                     const x * data)                     \
{                                                        \
    size_t nelements = 1;                                \
    for(int kk = 0; kk<ndim; kk++)                       \
    {                                                    \
        nelements *= shape[kk];                          \
    }                                                    \
    size_t element_size = sizeof(x);                     \
    return npio_write_raw(filename, ndim, shape,         \
                          (void *) data,                 \
                          nelements, element_size,       \
                          y);                            \
}                                                        \

NPIO_SAVE(f64, "<f8")
NPIO_SAVE(f32, "<f4")
NPIO_SAVE(i8, "|i1")
NPIO_SAVE(i16, "<i2")
NPIO_SAVE(i32, "<i4")
NPIO_SAVE(i64, "<i8")
NPIO_SAVE(u8, "|u1")
NPIO_SAVE(u16, "<u2")
NPIO_SAVE(u32, "<u4")
NPIO_SAVE(u64, "<u8")


/* DICTIONARY PARSER */

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

static int dp_parse(dp_t dp, const char * dict,
             const size_t dict_len, dptok_t * tok, const int tok_len)
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
                if(c == '(')
                {
                    nest ++;
                }

                state = 4;
                tok[dp.toknext].start = dp.pos;
                goto nextchar;
            }
        }
        if(state == 4) /* looking for end of value */
        {
            if( c == ')')
            {
                nest--;
            }
            if( c == '(')
            {
                nest++;
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
