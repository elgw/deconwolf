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

#include "ftab.h"

#include <assert.h>
#include <inttypes.h>
#define _USE_MATH_DEFINES
#include <math.h>


#include <stdlib.h>
#include <string.h>

#include <time.h>
#ifndef WINDOWS
#include <unistd.h>
#endif

typedef uint8_t u8;
typedef uint64_t u64;

/* Count the number of newlines in a file */
static size_t count_newlines(const char * fname)
{
    FILE * fid = fopen(fname, "r");
    assert(fid != NULL);
    size_t N = 0;
    int c = EOF;
    while( (c = fgetc(fid) ) != EOF )
    {
        if( c == '\n')
        {
            N++;
        }
    }
    fclose(fid);
    return N;
}

static void
trim_whitespace(char * str)
{
    if(str == NULL)
    {
        return;
    }
    size_t n = strlen(str);
    size_t first = 0;
    size_t last = n;
    if(n == 0)
    {
        return;
    }

    // Look for whitespaces from the beginning
    char * start = str;
    while(*start != '\0')
    {
        if(*start == ' ')
        {
            first++;
            start++;
        } else {
            break;
        }
    }

    // Look for whitespaces from the end
    char * end = str+n-1;
    while(end >= start)
    {
        char token = end[0];
        int ws = 0;
        switch(token)
        {
        case ' ':
            ws = 1;
            break;
        case '\n':
            ws = 1;
            break;
        case '\r':
            ws = 1;
            break;
        }
        if(ws)
        {
            last--;
            end--;
        } else {
            break;
        }
    }

    // copy from first to last
    //printf("first=%zu last=%zu\n", first, last);
    size_t wpos = 0;
    for(size_t kk = first; kk < last; kk++)
    {
        str[wpos++] = str[kk];
    }
    str[wpos] = '\0';
    return;
}

void ftab_head(ftab_t * T, int64_t n)
{

    if(n < 0)
    {
        return;
    }

    if(n > (int64_t) T->nrow)
    {
        n = T->nrow;
    }
    T->nrow = n;
}

void ftab_free(ftab_t * T)
{
    if(T == NULL)
    {
        return;
    }

    free(T->T);

    if(T->colnames != NULL)
    {
        for(size_t kk = 0; kk < T->ncol; kk++)
        {
            free(T->colnames[kk]);
        }
        free(T->colnames);
    }

    free(T);
    return;
}

int ftab_write_tsv(const ftab_t * T, const char * fname)
{
    FILE * fid = fopen(fname, "w");
    if(fid == NULL)
    {
        return EXIT_FAILURE;
    }
    int ret = ftab_print(fid, T, "\t");
    fclose(fid);
    return ret;
}

int ftab_write_csv(const ftab_t * T, const char * fname)
{
    FILE * fid = fopen(fname, "w");
    if(fid == NULL)
    {
        return EXIT_FAILURE;
    }
    int ret = ftab_print(fid, T, ",");
    fclose(fid);
    return ret;
}

int ftab_print(FILE * fid, const ftab_t * T, const char * sep)
{
    /* Write column names if they exist, otherwise col_1 etc */
    for(size_t cc = 0; cc<T->ncol; cc++)
    {
        int colname = 0;
        if(T->colnames != NULL)
        {
            if(T->colnames[cc] != NULL)
            {
                fprintf(fid, "%s", T->colnames[cc]);
                colname = 1;
            }
        }
        if(colname == 0)
        {
            fprintf(fid, "col_%zu", cc+1);
        }
        if(cc+1 != T->ncol)
        {
            fprintf(fid, "%s", sep);
        }
    }
    fprintf(fid, "\n");

    /* Write rows */
    for(size_t rr = 0; rr<T->nrow; rr++)
    {
        for(size_t cc = 0; cc<T->ncol; cc++)
        {
            fprintf(fid, "%f", T->T[rr*T->ncol + cc]);
            if(cc+1 != T->ncol)
            {
                fprintf(fid, "%s", sep);
            }
        }
        fprintf(fid, "\n");
    }
    return EXIT_SUCCESS;
}

ftab_t * ftab_new(int ncol)
{
    if(ncol < 1)
    {
        fprintf(stderr, "ftab_new requires at least 1 column\n");
        return NULL;
    }
    ftab_t * T = calloc(1, sizeof(ftab_t));
    assert(T != NULL);
    T->nrow = 0;
    T->ncol = ncol;
    T->nrow_alloc = 1024;
    T->T = calloc(T->ncol*T->nrow_alloc, sizeof(float));
    assert(T->T != NULL);
    return T;
}

static int
parse_col_names(ftab_t * T,
                const char * _line,
                const char * dlm)
{
    /* Figure out how many columns there are */
    if(strlen(_line) == 0)
    {
        return 0;
    }
    char * line = strdup(_line);
    assert(line != NULL);
    int ncol = 1;
    for(size_t kk = 0; kk<strlen(line); kk++)
    {
        if(line[kk] == dlm[0])
        {
            ncol++;
        }
    }
    T->ncol = ncol;

    /* Allocate memory */
    T->colnames = calloc(ncol, sizeof(char*));
    assert(T->colnames != NULL);

    char * line0 = line;
    /* Set columns */
    for(int kk = 0; kk<ncol; kk++)
    {
        // the tokens returned by strtok() are always nonempty strings
        char * token = strsep(&line, dlm);
        T->colnames[kk] = strdup(token);
        trim_whitespace(T->colnames[kk]);
        printf("read colnames: '%s'\n", T->colnames[kk]);
    }

    if(0){
        printf("Columns names:\n");
        for(int kk = 0; kk<ncol; kk++)
        {
            printf("%d '%s'\n", kk, T->colnames[kk]);
        }
        printf("\n");
    }
    //   exit(EXIT_FAILURE);
    free(line0);

    return ncol;
}

int ftab_get_col(const ftab_t * T, const char * name)
{
    int ret = -1;
    for(size_t kk = 0; kk<T->ncol; kk++)
    {
        if(T->colnames[kk] == NULL)
        {
            continue;
        }
        if(strcmp(T->colnames[kk], name) == 0)
        {
            if(ret != -1)
            {
                fprintf(stderr, "Warning multiple columns named '%s'", name);
            }
            ret = kk;
        }
    }
    //printf("%s is column %d\n", name, ret);
    return ret;
}

static int
parse_floats(char * l,
             float * row,
             int nval,
             const char * dlm)
{
    char * f = strsep(&l, dlm);
    if(f == NULL)
    {
        return 0;
    }
    //    printf("0/%d: %s\n", nval, f);
    row[0] = atof(f);
    for(int kk = 1; kk<nval; kk++)
    {
        f = strsep(&l, dlm);
        //      printf("%d: %s\n", kk, f);
        if(f == NULL)
        {
            return 0;
        }
        row[kk] = atof(f);
    }

    //   exit(EXIT_FAILURE);
    return 1;
}


static ftab_t *
ftab_from_dlm(const char * fname,
              const char * dlm)
{
    FILE * fid = fopen(fname, "r");
    if(fid == NULL)
    {
        fprintf(stderr, "Can not open %s\n", fname);
        return NULL;
    }
    ftab_t * T = calloc(1, sizeof(ftab_t));
    assert(T != NULL);
    // Get number of columns
    char * line = NULL;
    size_t len = 0;
    int read = getline(&line, &len, fid);
    if(strlen(line) == 0)
    {
        fprintf(stderr, "Empty header line\n");
        free(line);
        free(T);
        fclose(fid);
        return NULL;
    }

    int ncols = parse_col_names(T, line, dlm);

    // printf("%zu columns\n", T->ncol);
    size_t nrows = count_newlines(fname)+1;
    // printf("at most %zu lines\n", nrows);

    // Allocate memory
    T->T = calloc(nrows*ncols, sizeof(float));
    assert(T->T != NULL);
    T->ncol = ncols;
    T->nrow_alloc = nrows;
    // T->nrow = nrows; // set later

    // Read
    int row = 0;
    while ((read = getline(&line, &len, fid)) != -1)
    {
        if(strlen(line) == 0)
        {
            continue;
        }
        row += parse_floats(line, T->T + row*T->ncol, T->ncol, dlm);
    }
    T->nrow = row;
    //printf("Returning a %zux%zu table\n", T->nrow, T->ncol);
    free(line);
    fclose(fid);
    return T;
}


ftab_t * ftab_from_csv(const char * fname)
{
    return ftab_from_dlm(fname, ",");
}

ftab_t * ftab_from_tsv(const char * fname)
{
    return ftab_from_dlm(fname, "\t");
}

typedef struct{
    float value;
    size_t idx;
} ftab_sort_pair;

int ftab_sort_pair_cmp(const void * _A, const void * _B)
{
    ftab_sort_pair * A = (ftab_sort_pair*) _A;
    ftab_sort_pair * B = (ftab_sort_pair*) _B;

    //printf("idx: %zu, value: %f\n", A->idx, A->value);
    //printf("idx: %zu, value: %f\n", B->idx, B->value);
    //exit(1);
    float fa = A->value;
    float fb = B->value;
    if(fa == fb)
    {
        return 0;
    }
    if(fa < fb)
    {
        return 1;
    }
    return -1;
}

void ftab_sort(ftab_t * T, int col)
{
    if(col == -1)
    {
        fprintf(stderr, "ftab_sort: Can't use column %d for sorting\n", col);
        exit(EXIT_FAILURE);
    }
    ftab_sort_pair * P = calloc(T->nrow, sizeof(ftab_sort_pair));
    assert(P != NULL);
    /* Extract values from col %d and sort */

    for(size_t kk = 0; kk < T->nrow; kk++)
    {
        P[kk].idx = kk;
        P[kk].value = T->T[kk*T->ncol + col];
    }
    //printf("First: %f, Last: %f\n", P[0].value, P[T->nrow-1].value);
    //printf("Sorting dots %zu dots\n", T->nrow);
    qsort(P, T->nrow,
          sizeof(ftab_sort_pair),
          ftab_sort_pair_cmp);
    float * T2 = calloc(T->ncol*T->nrow_alloc, sizeof(float));
    assert(T2 != NULL);
    for(size_t kk = 0; kk< T->nrow; kk++)
    {
        size_t newpos = P[kk].idx*T->ncol;
        size_t oldpos = kk*T->ncol;
        memcpy(T2+oldpos,
               T->T + newpos,
               sizeof(float)*T->ncol);
    }
    free(P);
    free(T->T);
    T->T = T2;
    return;
}

void ftab_insert(ftab_t * T, float * row)
{
    assert(T != NULL);
    assert(row != NULL);
    if(T->nrow == T->nrow_alloc)
    {
        T->nrow_alloc += T->nrow_alloc*0.69;
        T->T = realloc(T->T, T->nrow_alloc*T->ncol*sizeof(float));
    }
    memcpy(T->T+T->ncol*T->nrow,
           row, T->ncol*sizeof(float));
    T->nrow++;
}




void ftab_set_colname(ftab_t * T, int col, const char * name)
{
    assert(name != NULL);
    // TODO check that only valid chars are used
    if(col < 0 || col >= (int) T->ncol)
    { return; }

    if(name == NULL)
    { return; }

    if(T->colnames == NULL)
    {
        T->colnames = calloc(T->ncol, sizeof(char*));
        assert(T->colnames != NULL);
    }

    free(T->colnames[col]);
    T->colnames[col] = strdup(name);
    return;
}


int ftab_set_coldata(ftab_t * T, int col, const float * data)
{
    if(T == NULL)
    {
        return EXIT_FAILURE;
    }
    if(col < 0 || (size_t) col >= T->ncol)
    {
        return EXIT_FAILURE;
    }

    float * C = T->T + col;
    for(size_t kk = 0; kk<T->nrow; kk++)
    {
        C[kk*T->ncol] = data[kk];
    }

    return EXIT_SUCCESS;
}


ftab_t * ftab_copy(const ftab_t * T)
{
    ftab_t * C = calloc(1, sizeof(ftab_t));
    assert(C != NULL);
    C->nrow = T->nrow;
    C->ncol = T->ncol;
    C->nrow_alloc = C->nrow;
    C->T = calloc(C->nrow*C->ncol, sizeof(float));
    if(C->T == NULL)
    {
        return NULL;
    }

    memcpy(C->T, T->T, C->nrow*C->ncol*sizeof(float));

    if(T->colnames != NULL)
    {
        C->colnames = calloc(C->ncol, sizeof(char*));
        assert(C->colnames != NULL);
        for(u64 kk = 0; kk < C->ncol; kk++)
        {
            if(T->colnames[kk] != NULL)
            {
                C->colnames[kk] = strdup(T->colnames[kk]);
                // TODO: gracefully tear down if NULL was returned.
            }
        }
    }
    return C;
}

ftab_t * ftab_concatenate_columns(const ftab_t * L , const ftab_t * R)
{
    if(L->nrow != R->nrow)
    {
        fprintf(stderr, "frab_concatenate_columns: "
                "ERROR: tables have different number of rows\n");
        return NULL;
    }
    size_t nrow = L->nrow;
    size_t ncol = L->ncol + R->ncol;

    ftab_t * T = ftab_new(ncol);
    assert(T != NULL);

    for(size_t kk = 0; kk < L->ncol; kk++)
    {
        ftab_set_colname(T, kk, L->colnames[kk]);
    }
    for(size_t kk = 0; kk < R->ncol; kk++)
    {
        ftab_set_colname(T, kk+L->ncol, R->colnames[kk]);
    }
    free(T->T);
    T->T = calloc(nrow*ncol, sizeof(float));
    assert(T->T != NULL);
    T->ncol = ncol;
    T->nrow = nrow;
    T->nrow_alloc = nrow;

    /* Insert data from L */
    for(size_t kk = 0; kk < L->nrow; kk++)
    {
        for(size_t ll = 0; ll < L->ncol; ll++)
        {
            T->T[kk*ncol + ll] = L->T[kk*L->ncol + ll];
        }
    }

    /* Insert data from R */
    for(size_t kk = 0; kk < R->nrow; kk++)
    {
        for(size_t ll = 0; ll < R->ncol; ll++)
        {
            T->T[kk*ncol + ll + L->ncol] = R->T[kk*R->ncol + ll];
        }
    }

    return T;
}

ftab_t *
ftab_concatenate_rows(const ftab_t * Top, const ftab_t * Down)
{

    if(Top == NULL && Down == NULL)
    {
        return NULL;
    }
    if(Top == NULL)
    {
        return ftab_copy(Down);
    }
    if(Down == NULL)
    {
        return ftab_copy(Top);
    }

    if(Top->ncol != Down->ncol)
    {
        fprintf(stderr,
                "ftab concatenate_rows: error: Tables does "
                "not have the same number of columns\n");
        return NULL;
    }
    // TODO: Check column names

    ftab_t * concat = calloc(1, sizeof(ftab_t));
    assert(concat != NULL);
    concat->ncol = Top->ncol;
    concat->nrow = Top->nrow + Down->nrow;
    concat->nrow_alloc = concat->nrow;
    concat->T = calloc(concat->nrow*concat->ncol, sizeof(float));
    if(concat->T == NULL)
    {
        free(concat);
        return NULL;
    }
    memcpy(concat->T,
           Top->T,
           Top->nrow*Top->ncol*sizeof(float));
    memcpy(concat->T+Top->nrow*Top->ncol,
           Down->T,
           Down->nrow*Down->ncol*sizeof(float));

    // TODO: Set column names
    return concat;
}

ftab_t *
ftab_new_from_data(int nrow, int ncol, const float * data)
{
    ftab_t * T = calloc(1, sizeof(ftab_t));
    if(T == NULL)
    {
        return T;
    }
    T->ncol = ncol;
    T->nrow = nrow;
    T->nrow_alloc = T->nrow;
    T->T = calloc(nrow*ncol, sizeof(float));
    if(T->T == NULL)
    {
        free(T);
        return NULL;
    }
    memcpy(T->T, data, nrow*ncol*sizeof(float));
    return T;
}

char * tempfilename()
{
    char * fname = calloc(1024, 1);
    assert(fname != NULL);
#ifdef WINDOWS
    printf("TODO: Windows functionality\n");
    sprintf(fname, "ftab-random");
#else
    sprintf(fname, "/tmp/ftab-XXXXXX");
    int filedes = mkstemp(fname);
    close(filedes);
#endif
    return fname;
}

int ftab_compare(const ftab_t * A, const ftab_t * B)
{
    if(A == NULL)
    {
        return 1;
    }
    if(B == NULL)
    {
        return 1;
    }
    if(A->ncol != B->ncol)
    {
        return 1;
    }
    if(A->nrow != B->nrow)
    {
        return 1;
    }

    int equal_colnames = 1;
    if(A->colnames == NULL)
    {
        if(B->colnames != NULL)
        {
            equal_colnames = 0;
        }

    }
    if(B->colnames == NULL)
    {
        if(A->colnames != NULL)
        {
            equal_colnames = 0;
        }
    }
    if(A->colnames != NULL && B->colnames != NULL)
    {
        for(size_t kk = 0; kk<A->ncol; kk++)
        {
             if(strcmp(A->colnames[kk], A->colnames[kk]))
             {
                 equal_colnames = 0;
            }
        }
    }
    if(equal_colnames == 0)
    {
        return 1;
    }

    if(A->nrow == 0 || A->nrow == 0)
    {
        return 0;
    }
    // Compare data
    return memcmp(A->T, B->T, A->ncol*A->nrow*sizeof(float));
}

int ftab_ut(int argc, char ** argv)
{
    if(argc == 1)
    {
        printf("Running some self tests.\n");
        printf("To test on a specific file, use:\n");
        printf("%s file.csv\n", argv[1]);
        printf("\n");
    }
    if(argc > 1)
    {
        printf("Reading %s as CSV\n", argv[1]);
        ftab_t * T = ftab_from_csv(argv[1]);
        printf("\n");
        if(T == NULL)
        {
            printf("Unable to read the file\n");
            return EXIT_FAILURE;
        }
        printf("Table size: %lu x %lu\n", T->nrow, T->ncol);
        if(T->colnames != NULL)
        {
            printf("Columns names\n");
            for(int kk = 0; kk < (int) T->ncol; kk++)
            {
                printf("%2d '%s'\n", kk+1, T->colnames[kk]);
            }
        }
        if(T->nrow > 0)
        {
            printf("First row:\n");
            for(int kk = 0; kk < (int) T->ncol; kk++)
            {
                printf("%f\t", T->T[kk]);
            }
            printf("\n");
        }

        char * fname = tempfilename();
        printf("Temporary file: %s\n", fname);
        ftab_write_tsv(T, fname);


        ftab_t * T2 = ftab_from_tsv(fname);
#ifndef WINDOWS
        unlink(fname);
#endif
        free(fname);

        int status = ftab_compare(T, T2);
        ftab_free(T);
        ftab_free(T2);
        if(status == 0)
        {
            printf("OK! File can be written and read back\n");
            return EXIT_SUCCESS;
        } else {
            printf("Failed to write and read back this file\n");
            return EXIT_FAILURE;
        }
    }

    ftab_t * T = ftab_new(4);
    assert(T != NULL);
    printf("T: %zu x %zu\n", T->nrow, T->ncol);
    ftab_set_colname(T, 0, "x");
    ftab_set_colname(T, 1, "y");
    ftab_set_colname(T, 2, "z");
    ftab_set_colname(T, 3, "value");
    ftab_print(stdout, T, "\t");

    float row[4] = {1, 2, 3, 1.23};
    ftab_insert(T, row);
    ftab_print(stdout, T, "\t");

    /* Save and read */
    char * fname = tempfilename();

    printf("Temporary file: %s\n", fname);
    ftab_write_tsv(T, fname);


    ftab_t * T2 = ftab_from_tsv(fname);
    int status = ftab_compare(T, T2);
    if( status != 0 )
    {
        printf("Test failed\n");
    }

#ifndef WINDOWS
    unlink(fname);
#endif
    free(fname);
    ftab_print(stdout, T, "\t");
    ftab_free(T);
    ftab_free(T2);

    return EXIT_SUCCESS;
}

void
ftab_subselect_rows(ftab_t * tab, const u8 * selection)
{
    u64 nsel = 0;
    for(u64 kk = 0; kk < tab->nrow; kk++)
    {
        if(selection[kk] > 0)
        {
            if(kk != nsel)
            {
                memcpy(tab->T + nsel*tab->ncol,
                       tab->T + kk*tab->ncol,
                       tab->ncol*sizeof(float));
            }
            nsel++;
        }
    }
    tab->nrow = nsel;
    return;
}
