#include "kdtree.h"
#include <assert.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#include "ftab.h"
#include "dw_version.h"
#include "dw_util.h"

typedef struct{
    int verbose;
    /* Largest shift, in pixels, to consider between
     * the point clouds
     */
    double capture_distance;
    double sigma;
    int optind;
    /* Most number of points to use from each file.
     * A value of 0 is interpreted as "use all"
     */
    int npoint;
} opts;

static opts * opts_new()
{
    opts * s = calloc(1, sizeof(opts));
    assert(s != NULL);
    s->verbose = 1;
    s->capture_distance = 10;
    s->sigma = 1;
    return s;
}

static void opts_free(opts * s)
{
    free(s);
}

static void usage(void)
{
    printf("usage: dw align-dots [<options>] file1.tsv file2.tsv\n");
    printf("\n");
    printf("Options:\n");
    printf("--distance d, -d d\n"
           "\tSet the capture distance, i.e. largest expected shift in pixels\n");
    printf("--sigma s, -s s\n"
           "\tSet the size of the KDE used to identify the peak\n");
    printf("--npoint n, -n n\n"
           "\t set the maximum number of points to use from each file\n");
    printf("\n");
    return;
}

static void argparsing(int argc, char ** argv, opts * s)
{

    struct option longopts[] = {
        {"help", no_argument, NULL, 'h'},
        {"verbose", required_argument, NULL, 'v'},
        {"distance", required_argument, NULL, 'd'},
        {"npoint", required_argument, NULL, 'n'},
        {"sigma", required_argument, NULL, 's'},
        {NULL, 0, NULL, 0}};

    int ch;
    while( (ch = getopt_long(argc, argv,
                             "d:hn:s:v:",
                             longopts, NULL)) != -1)
    {
        switch(ch)
        {
        case 'd':
            s->capture_distance = atof(optarg);
            break;
        case 'h':
            usage();
            exit(EXIT_SUCCESS);
            break;
        case 'n':
            s->npoint = atol(optarg);
            break;
        case 's':
            s->sigma = atof(optarg);
            break;
        case 'v':
            s->verbose = atoi(optarg);
            break;
        default:
            exit(EXIT_FAILURE);
        }
    }
    s->optind = optind;
    return;
}

double * get_dots(opts * s, ftab_t * T)
{
    /* Extract the dots will look for fit_x, fit_y, ...
     * if not found, look for columns called x, y, z */

    int cx, cy, cz;
    cx = ftab_get_col(T, "fit_x");
    cy = ftab_get_col(T, "fit_y");
    cz = ftab_get_col(T, "fit_z");
    if(cx < 0)
    {
        printf("Warning: Could not find any column named 'fit_x'\n");
        cx = ftab_get_col(T, "x");
        cy = ftab_get_col(T, "y");
        cz = ftab_get_col(T, "z");
        if(cx < 0)
        {
            printf("Error: Could not find any column named 'x'\n");
        }
    }

    if(cx < 0 || cy < 0 || cz < 0)
    {
        printf("Error. Did not find the expected columns.\n");
        return NULL;
    }

    double * X = calloc(3*T->nrow, sizeof(double));
    assert(X != NULL);
    for(size_t kk = 0; kk< T->nrow; kk++)
    {
        float * row = T->T + kk*T->ncol;
        X[3*kk + 0] = row[cx];
        X[3*kk + 1] = row[cy];
        X[3*kk + 2] = row[cz];
    }
    return X;
}

/* Dynamic vector array */
struct dvarray {
    double * data;
    size_t n_used;
    size_t n_alloc;
};

static void dvarray_n_more(struct dvarray * A, size_t nmore)
{
    if(A->n_used + nmore >= A->n_alloc)
    {
        size_t new_size = A->n_alloc + nmore;
        if(new_size < 1.2 *A->n_alloc)
        {
            new_size = 1.2*A->n_alloc;
        }
        double * t = realloc(A->data, 3*new_size*sizeof(double));
        if(t == NULL) // This makes fanalyzer happy
        {
            free(A->data);
        }
        assert(t != NULL);
        A->data = t;


        assert(A->data != NULL);
        A->n_alloc = new_size;
    }
}

static void dvarray_insert_vector(struct dvarray * A, const double * X)
{
    dvarray_n_more(A, 1);
    memcpy(A->data + 3*A->n_used,
           X,
           3*sizeof(double));
    A->n_used++;
    return;
}

static void dvarray_free(struct dvarray * A)
{
    free(A->data);
    free(A);
    return;
}

struct dvarray * dvarray_new(size_t n)
{
    assert(n > 0);
    struct dvarray * A = calloc(1, sizeof(struct dvarray));
    assert(A != NULL);
    A->data = calloc(3*n, sizeof(double));
    assert(A->data != NULL);
    A->n_alloc = n;
    return A;
}


static void align_dots(const opts * s,
                       const double * XA, size_t nXA,
                       const double * XB, size_t nXB)
{
    kdtree_t * TA = kdtree_new(XA, nXA, 20);

    /* Detect all pairs of points within some distance */
    struct dvarray * arr = dvarray_new(1024);
    size_t nfound_total = 0;
    for(size_t kk = 0; kk < nXB; kk++)
    {
        const double * Q = XB + 3*kk; // Query point
        size_t nfound = 0;
        size_t * A_idx = kdtree_query_radius(TA, Q,
                                             s->capture_distance,
                                             &nfound);

        for(size_t ll = 0; ll < nfound; ll++)
        {
            const double * P = XA + 3*A_idx[ll];
            double D[3] = {0};
            D[0] = Q[0] - P[0];
            D[1] = Q[1] - P[1];
            D[2] = Q[2] - P[2];
            dvarray_insert_vector(arr, D);
        }
        nfound_total += nfound;
        free(A_idx);
    }
    printf("Found %zu pairs within the capture radius (%f)\n",
           nfound_total, s->capture_distance);
    kdtree_free(TA); TA = NULL;

    kdtree_t * TD = kdtree_new(arr->data, arr->n_used, 10);
    assert(TD != NULL);
//    kdtree_print_info(TD);
    dvarray_free(arr);
    arr = NULL;

    double maxkde = 0;
    double maxpos[3] = {0};
    printf("Grid search\n");
    double rs = s->capture_distance;
    for(double x = -rs; x <= rs; x+=0.5) {
        for(double y = -rs; y <= rs; y+=0.5) {
            for(double z = -rs; z <= rs; z+=0.5) {
                double P[] = {x, y, z};
                double v = kdtree_kde(TD, P, s->sigma, 0);
                if(v > maxkde)
                {
                    maxkde = v;
                    memcpy(maxpos, P, 3*sizeof(double));
                }
            }
        }
    }

    printf("     Grid search: (% .2f, % .2f, % .2f) (kde=%.1f)\n",
           maxpos[0], maxpos[1], maxpos[2], maxkde);

    /* Refinement over the grid search */
    rs = 2*s->sigma; // Region size
    while(rs > 1e-3)
    {
    double center[3];
    memcpy(center, maxpos, 3*sizeof(double));
    for(double x = -rs; x <= rs; x+= rs/5.0) {
        for(double y = -rs; y <= rs; y+= rs/5.0) {
            for(double z = -rs; z <= rs; z+= rs/5.0) {
                double P[] = {
                    x+center[0],
                    y+center[1],
                    z+center[2]};
                double v = kdtree_kde(TD, P, s->sigma, 0);
                //printf("%f, %f, %f -> %f\n", P[0], P[1], P[2], v);
                if(v > maxkde)
                {
                    maxkde = v;
                    memcpy(maxpos, P, 3*sizeof(double));
                }
            }
        }
    }
    rs /= 2.0;
    }

    printf("Refined position: (% .2f, % .2f, % .2f) (kde=%.1f)\n",
           maxpos[0],maxpos[1], maxpos[2], maxkde);

    double P[] = {0,0,0};
    printf("At origo, kde=%.2f\n", kdtree_kde(TD, P, s->sigma, 0));

    kdtree_free(TD);
    TD = NULL;
    return;

}

int dw_align_dots(int argc, char ** argv)
{

    opts * s = opts_new();
    argparsing(argc, argv, s);

    if(s->optind + 2 != argc)
    {
        printf("Please pass exactly two tsv files to process\n");
        opts_free(s);
        return EXIT_FAILURE;
    }

    /* Load the tables */
    if(s->verbose > 1)
    {
        printf("Reading %s\n", argv[optind]);
    }
    ftab_t * TA = ftab_from_tsv(argv[optind]);
    double * XA = get_dots(s, TA);
    size_t nXA = TA->nrow;
    ftab_free(TA); TA = NULL;
    if(XA == NULL)
    {
        return EXIT_FAILURE;
    }

    if(s->verbose > 1)
    {
        printf("Reading %s\n", argv[optind+1]);
    }
    ftab_t * TB = ftab_from_tsv(argv[optind+1]);
    double * XB = get_dots(s, TB);
    size_t nXB = TB->nrow;
    ftab_free(TB); TB = NULL;
    if(XB == NULL)
    {
        return EXIT_FAILURE;
    }

    if(s->verbose > 0)
    {
        printf("Will align %zu dots vs %zu\n", nXA, nXB);
    }

    align_dots(s, XA, nXA, XB, nXB);

    free(XA);
    free(XB);

    /* Present the result */

    opts_free(s);
    s = NULL;
    return EXIT_SUCCESS;
}
