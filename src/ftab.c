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
static size_t count_newlines(char * fname);

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
    free(T);
    return;
}

ftab_t * ftab_new(int ncol)
{
    ftab_t * T = malloc(sizeof(ftab_t));
    T->ncol = ncol;
    T->nrow = 0;
    T->nrow_alloc = 1024;
    T->T = malloc(T->ncol*T->nrow_alloc*sizeof(float));
    return T;
}

int parse_col_names(ftab_t * T, char * line)
{
    /* Figure out how many columns there are */
    if(strlen(line) == 0)
    {
        return 0;
    }
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
    return ncol;
}

int ftab_get_col(ftab_t * T, char * name)
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
    return ret;
}

ftab_t * ftab_from_tsv(char * fname)
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
        return NULL;
    }

    int ncols = parse_col_names(T, line);

    // printf("%d columns\n", ncols);
    size_t nrows = count_newlines(fname)+1;
    // printf("at most %zu lines\n", nrows);

    // Allocate memory
    T->T = malloc(nrows*ncols*sizeof(float));
    T->ncol = ncols;
    T->nrow_alloc = nrows;
    // T->nrow = nrows; // set later

    // Read
    int row = 0;
    while ((read = getline(&line, &len, fid)) != -1) {
        if(strlen(line) == 0 && line[0] == '#')
        {
            continue;
        }
        row += parse_floats(line, T->T + row*T->ncol, T->ncol);
    }
    T->nrow = row-1;
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
    static size_t count_newlines(char * fname)
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
