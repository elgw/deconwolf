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
    size_t npoint;
    /* Name of file to write to (optional) */
    char * outfile;
    /* Overwrite outfile ?  */
    int overwrite;
    /* Append to the outfile ? */
    int append;

    /* For storing the result */
    double dx;
    double dy;
    double dz;
    double kde;
    /* Some measurement on how certain the result is */
    double goodness;
} opts;

// Golden ratio
const double PHI = 1.618033988749;

// The 12 vertices of an icosahedron;
const double icosahedron[36] = {0, PHI, 1,
                                0, -PHI, 1,
                                0, PHI, -1,
                                0, -PHI, -1,
                                1, 0, PHI,
                                -1, 0, PHI,
                                1, 0, -PHI,
                                -1, 0, -PHI,
                                PHI, 1, 0,
                                -PHI, 1, 0,
                                PHI, -1, 0,
                                -PHI, -1, 0};

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
    free(s->outfile);
    free(s);
}

static void opts_print(opts * s, FILE * fid)
{
    fprintf(fid, "capture distance: %f\n", s->capture_distance);
    fprintf(fid, "       KDE sigma: %f\n", s->sigma);
    return;
}

static void usage(void)
{
    printf("usage: dw align-dots [<options>] file1.tsv file2.tsv\n");
    printf("\n");
    printf("Options:\n");
    printf("--radius d, -d d\n"
           "\tSet the capture radius, i.e. largest expected shift in pixels\n");
    printf("--sigma s, -s s\n"
           "\tSet the size of the KDE used to identify the peak\n");
    printf("--npoint n, -n n\n"
           "\tset the maximum number of points to use from each file\n");
    printf("--out file, -o file\n"
           "\tWhere to write the results.\n");
    printf("--overwrite\n"
           "\tOverwrite the destination file if it exists\n");
    printf("--append, -a\n"
           "\tAppend to the output file\n");
    printf("\n");
    printf("Output columns:\n");
    printf(" 1. dx\n");
    printf(" 2. dy\n");
    printf(" 3. dz\n");
    printf(" 4. kde\n");
    printf(" 5. goodness\n");
    printf(" 6. reference file\n");
    printf(" 7. other file\n");
    printf("\n");
    return;
}

static void argparsing(int argc, char ** argv, opts * s)
{

    struct option longopts[] = {
        {"help", no_argument, NULL, 'h'},
        {"verbose", required_argument, NULL, 'v'},
        {"radius",   required_argument, NULL, 'd'},
        {"npoint", required_argument, NULL, 'n'},
        {"out",    required_argument, NULL, 'o'},
        {"append", no_argument, NULL, 'a'},
        {"sigma", required_argument, NULL, 's'},
        {"overwrite", no_argument, NULL, 'w'},
        {NULL, 0, NULL, 0}};

    int ch;
    while( (ch = getopt_long(argc, argv,
                             "ad:hn:o:s:v:w",
                             longopts, NULL)) != -1)
    {
        switch(ch)
        {
        case 'a':
            s->append = 1;
            break;
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
        case 'o':
            free(s->outfile);
            s->outfile = strdup(optarg);
            break;
        case 's':
            s->sigma = atof(optarg);
            break;
        case 'v':
            s->verbose = atoi(optarg);
            break;
        case 'w':
            s->overwrite = 1;
            break;
        default:
            exit(EXIT_FAILURE);
        }
    }

    if(s->outfile != NULL)
    {
        if(s->append == 0 && s->overwrite == 0)
        {
            if(dw_file_exist(s->outfile))
            {
                fprintf(stderr, "Error: The output file does already exists\n");
                fprintf(stderr, "Consider --overwrite or --append\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    s->optind = optind;
    return;
}

double * get_dots(opts * s, ftab_t * T)
{
    if(s->verbose > 1)
    {
        printf("Looking for columns f_x, f_y and f_z\n");
    }

    int cx, cy, cz;
    cx = ftab_get_col(T, "f_x");
    cy = ftab_get_col(T, "f_y");
    cz = ftab_get_col(T, "f_z");
    if(cx < 0)
    {
        printf("Warning: Could not find any column named 'f_x'\n");
        if(s->verbose > 1)
        {
            printf("Looking for columns x, y and z\n");
        }
        cx = ftab_get_col(T, "x");
        cy = ftab_get_col(T, "y");
        cz = ftab_get_col(T, "z");
        if(cx < 0)
        {
            fprintf(stderr,
                    "Error: Could not find any column named 'f_x' or x'\n");
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

static double v3_norm(const double * X)
{
    return sqrt( pow(X[0], 2) + pow(X[1], 2) + pow(X[2], 2));
}


static void align_dots(opts * s,
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
    if(s->verbose > 0)
    {
        printf("Found %zu pairs within the capture radius (%f)\n",
               nfound_total, s->capture_distance);
    }
    kdtree_free(TA); TA = NULL;

    kdtree_t * TD = kdtree_new(arr->data, arr->n_used, 10);
    assert(TD != NULL);
    //    kdtree_print_info(TD);
    dvarray_free(arr);
    arr = NULL;

    double maxkde = 0;
    double maxpos[3] = {0};
    if(s->verbose > 1)
    {
        printf("Grid search\n");
    }
    double rs = s->capture_distance;
    double grid_average = 0;
    double n_grid = 0;
    for(double x = -rs; x <= rs; x+=0.5) {
        for(double y = -rs; y <= rs; y+=0.5) {
            for(double z = -rs; z <= rs; z+=0.5) {
                double P[] = {x, y, z};
                if(v3_norm(P) > rs){
                    continue;
                }
                double v = kdtree_kde(TD, P, s->sigma, 0);
                grid_average += v;
                n_grid++;
                if(v > maxkde)
                {
                    maxkde = v;
                    memcpy(maxpos, P, 3*sizeof(double));
                }
            }
        }
    }
    grid_average /= n_grid;

    if(s->verbose > 1)
    {
        printf("     Grid search: (% .2f, % .2f, % .2f) (kde=%.1f)\n",
               maxpos[0], maxpos[1], maxpos[2], maxkde);
    }

    /* Refinement over the grid search */
    rs = 1.0*s->sigma; // Region size
    while(rs > 1e-4)
    {
        double center[3];
        memcpy(center, maxpos, 3*sizeof(double));
        for(double x = -rs; x <= rs; x += rs/7.0) {
            for(double y = -rs; y <= rs; y += rs/7.0) {
                for(double z = -rs; z <= rs; z += rs/7.0) {
                    // Only consider inside a sphere
                    double X[] = {x, y, z};
                    if(v3_norm(X) > rs){
                        continue;
                    }
                    // Shift by the best position
                    double P[] = {
                        x+center[0],
                        y+center[1],
                        z+center[2]};

                    double v = kdtree_kde(TD, P, s->sigma, 0);

                    if(v > maxkde)
                    {
                        //printf("%f, %f, %f -> %f\n", P[0], P[1], P[2], v);
                        maxkde = v;
                        memcpy(maxpos, P, 3*sizeof(double));
                    }
                }
            }
        }
        rs /= 2.0;
    }

    double sphere_radius = 3.0*s->sigma;
    double kde_sphere = 0;
    for(size_t kk = 0 ; kk < 12 ; kk++)
    {
        double P[3] = {icosahedron[3*kk],
                       icosahedron[3*kk+1],
                       icosahedron[3*kk+2]};
        // Give specific radius
        double nr = sphere_radius / v3_norm(P);
        P[0] *= nr;
        P[1] *= nr;
        P[2] *= nr;
        // Place at the maxpos
        P[0] += maxpos[0];
        P[1] += maxpos[1];
        P[2] += maxpos[2];
        double v = kdtree_kde(TD, P, s->sigma, 0);
        if(v > kde_sphere)
        {
            kde_sphere = v;
        }
    }
    if(kde_sphere > maxkde)
    {
        fprintf(stderr, "ERROR: Unexpected result\n");
    }

    double P[] = {0,0,0};
    double origo_kde = kdtree_kde(TD, P, s->sigma, 0);
    if(s->verbose > 0)
    {
        printf("At origo, kde=%.2f (without any correction)\n", origo_kde);
    }
    if(s->verbose > 0)
    {
        printf("Refined position: (% .2f, % .2f, % .2f) (kde=%.1f)\n",
               maxpos[0],maxpos[1], maxpos[2], maxkde);
    }
    if(maxkde < origo_kde)
    {
        fprintf(stderr, "Error: unexpected result!\n");
    } else {
        printf(" %.1fX better than in the surrounding sphere of radius %f\n",
               maxkde / kde_sphere, sphere_radius );
    }
    if(s->verbose > 1)
    {
        printf("Grid average: %f\n", grid_average);
    }


    s->dx = maxpos[0];
    s->dy = maxpos[1];
    s->dz = maxpos[2];
    s->kde = maxkde;
    s->goodness = maxkde / kde_sphere;
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

    if(s->verbose > 1)
    {
        opts_print(s, stdout);
    }

    /* Load the tables */
    if(s->verbose > 1)
    {
        printf("Reading %s\n", argv[optind]);
    }
    ftab_t * TA = ftab_from_tsv(argv[optind]);
    if(TA == NULL)
    {
        opts_free(s);
        exit(EXIT_FAILURE);
    }
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
    if(TB == NULL)
    {
        opts_free(s);
        free(XA);
        exit(EXIT_FAILURE);
    }
    double * XB = get_dots(s, TB);
    size_t nXB = TB->nrow;
    ftab_free(TB); TB = NULL;
    if(XB == NULL)
    {
        return EXIT_FAILURE;
    }

    /* Limit the number of dots to use? */
    if(s->npoint > 0)
    {
        nXA > s->npoint ? nXA = s->npoint : 0;
        nXB > s->npoint ? nXB = s->npoint : 0;
    }

    if(s->verbose > 0)
    {
        printf("Will align %zu dots vs %zu\n", nXA, nXB);
    }

    align_dots(s, XA, nXA, XB, nXB);

    free(XA);
    free(XB);

    /* Present the result */
    if(s->outfile != NULL)
    {
        FILE * fid = NULL;
        if(s->append)
        {
            fid = fopen(s->outfile, "a");
        } else {
            fid = fopen(s->outfile, "w");
        }
        fprintf(fid, "%f, %f, %f", s->dx, s->dy, s->dz);
        fprintf(fid, ", %f, %f", s->kde, s->goodness);
        fprintf(fid, ", '%s'", argv[optind]);
        fprintf(fid, ", '%s'", argv[optind+1]);
        fprintf(fid, "\n");
    }

    opts_free(s);
    s = NULL;
    return EXIT_SUCCESS;
}
