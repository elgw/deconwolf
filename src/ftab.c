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

/* Forward declarations for non exported functions */
static int parse_floats(char * l, float * row, int nval);
static size_t count_newlines(const char * fname);

void ftab_free(ftab_t * T)
{
    if(T == NULL)
    {
        return;
    }
    if(T->T != NULL)
    {
        free(T->T);
    }

    if(T->colnames != NULL)
    {
        for(size_t kk = 0; kk<T->ncol; kk++)
        {
            if(T->colnames[kk] != NULL)
            {
                free(T->colnames[kk]);
            }
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
    int ret = ftab_print(fid, T);
    fclose(fid);
    return ret;
}

int ftab_print(FILE * fid, const ftab_t * T)
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
            fprintf(fid, "\t");
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
                fprintf(fid, "\t");
            }
        }
        fprintf(fid, "\n");
    }
    return EXIT_SUCCESS;
}

ftab_t * ftab_new(int ncol)
{
    ftab_t * T = malloc(sizeof(ftab_t));
    T->ncol = ncol;
    T->nrow = 0;
    T->nrow_alloc = 1024;
    T->T = malloc(T->ncol*T->nrow_alloc*sizeof(float));
    T->colnames = NULL;
    return T;
}

static int parse_col_names(ftab_t * T, const char * _line)
{
    /* Figure out how many columns there are */
    if(strlen(_line) == 0)
    {
        return 0;
    }
    char * line = strdup(_line);
    int ncol = 1;
    for(size_t kk = 0; kk<strlen(line); kk++)
    {
        if(line[kk] == '\t')
        {
            ncol++;
        }
    }
    T->ncol = ncol;

    /* Allocate memory */
    T->colnames = malloc(ncol*sizeof(char*));

    /* Set columns */
    char * f = strtok(line, "\t");
    assert(f != NULL); /* we already know that strlen > 0 */
    T->colnames[0] = strdup(f);
    for(int kk = 1; kk<ncol; kk++)
    {
        f = strtok(NULL, "\t");
        if(strlen(f) > 0)
        {
            if(f[strlen(f)-1] == '\n')
            {
                f[strlen(f)-1] = '\0';
            }
        }
        T->colnames[kk] = strdup(f);
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
    free(line);
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
    printf("%s is column %d\n", name, ret);
    return ret;
}

ftab_t * ftab_from_tsv(const char * fname)
{
    FILE * fid = fopen(fname, "r");
    if(fid == NULL)
    {
        fprintf(stderr, "Can not open %s\n", fname);
        return NULL;
    }
    ftab_t * T = malloc(sizeof(ftab_t));
    // Get number of columns
    char * line = NULL;
    size_t len = 0;
    int read = getline(&line, &len, fid);
    if(strlen(line) == 0)
    {
        fprintf(stderr, "Empty header line\n");
        if(line != NULL)
        {
            free(line);
        }
        return NULL;
    }

    int ncols = parse_col_names(T, line);

    printf("%zu columns\n", T->ncol);
    size_t nrows = count_newlines(fname)+1;
    printf("at most %zu lines\n", nrows);

    // Allocate memory
    T->T = malloc(nrows*ncols*sizeof(float));
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
        row += parse_floats(line, T->T + row*T->ncol, T->ncol);
    }
    T->nrow = row;
    //printf("Returning a %zux%zu table\n", T->nrow, T->ncol);
    free(line);
    fclose(fid);
    return T;
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
    ftab_sort_pair * P = malloc(T->nrow*sizeof(ftab_sort_pair));
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
    float * T2 = malloc(T->ncol*T->nrow_alloc*sizeof(float));
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

/* Count the number of newlines in a file */
    static size_t count_newlines(const char * fname)
{
    FILE * fid = fopen(fname, "r");
    assert(fid != NULL);
    size_t N = 0;
    char c = EOF;
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


static int parse_floats(char * l, float * row, int nval)
{
    char * f = strtok(l, "\t");
    if(f == NULL)
    {
        return 0;
    }
//    printf("0/%d: %s\n", nval, f);
    row[0] = atof(f);
    for(int kk = 1; kk<nval; kk++)
    {
        f = strtok(NULL, "\t");
        //      printf("%d: %s\n", kk, f);
        if(f == NULL)
        {
            return 0;
        }
        row[kk] = atof(f);
    }

    if(0){
        for(int kk = 0; kk<nval; kk++)
        {
            printf("%f ", row[kk]);
        }
        printf("\n");
    }
    //   exit(EXIT_FAILURE);
    return 1;
}

void ftab_set_colname(ftab_t * T, int col, const char * name)
{
    // TODO check that only valid chars are used
    if(col < 0 || col >= (int) T->ncol)
    {
        return;
    }

    if(name == NULL)
    {
        return;
    }

    if(T->colnames == NULL)
    {
        T->colnames = malloc(T->ncol*sizeof(char*));
        for(size_t kk = 0; kk<T->ncol; kk++)
        {
            T->colnames[kk] = NULL;
        }
    }

    if(T->colnames[col] != NULL)
    {
        free(T->colnames[col]);
    }
    T->colnames[col] = strdup(name);
}

int ftab_ut()
{
    ftab_t * T = ftab_new(4);
    printf("T: %zu x %zu\n", T->nrow, T->ncol);
    ftab_set_colname(T, 0, "x");
    ftab_set_colname(T, 1, "y");
    ftab_set_colname(T, 2, "z");
    ftab_set_colname(T, 3, "value");
    ftab_print(stdout, T);

    float row[4] = {1, 2, 3, 1.23};
    ftab_insert(T, row);
    ftab_print(stdout, T);

    /* Save and read */
    char * fname = malloc(1024);
    sprintf(fname, "/tmp/ftab-XXXXXX");
    int filedes = mkstemp(fname);
    close(filedes);

    printf("Temporary file: %s\n", fname);
    ftab_write_tsv(T, fname);


    ftab_t * T2 = ftab_from_tsv(fname);
    assert(T2 != NULL);
    assert(T2->ncol == T->ncol);
    assert(T2->nrow == T->nrow);
    for(size_t kk = 0; kk<T->ncol; kk++)
    {
        assert(strcmp(T->colnames[kk], T2->colnames[kk]) == 0);
    }


    unlink(fname);
    free(fname);
    ftab_print(stdout, T);
    ftab_free(T);
    ftab_free(T2);
    return EXIT_SUCCESS;
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

    float * C = T->T+col;
    for(size_t kk = 0; kk<T->nrow; kk++)
    {
        C[kk*T->nrow] = data[kk];
    }

    return EXIT_SUCCESS;
}
