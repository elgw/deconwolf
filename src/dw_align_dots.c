#include "kdtree.h"
#include <assert.h>
#include <getopt.h>
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ftab.h"
#include "dw_version.h"
#include "dw_util.h"
#include "npio.h"
#include "qalign.h"
#include "quickselect.h"

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

    /* Magnification factor for the list of dots */
    int mag_is_set;
    double mag1;
    double mag2;

    /* Multiplier for the z-component to match the localization
       accuracy in the lateral plane and the axial direction */
    double weightz;

    /* Initial shifts for the dots */
    char * fname_d1;
    char * fname_d2;
    double delta1[3];
    double delta2[3];

    /* Rotation about the z-axis, i.e. in the plane */
    double rotz;

    /* For storing the result */
    union{
        double delta[3];
        struct{
            double dx;
            double dy;
            double dz;
        };
    };
    double kde;
    /* Some measurement on how certain the result is */
    double goodness;
} opts;


/* Poisson Probability Mass Function (PMF)
 * I.e. lambda^k*exp(-lambda)/k!
 */
static double poisson_pmf(double k, double lambda)
{
    return exp(k*log(lambda)-lambda-lgamma(k+1));
}

static double poisson_cdf(double k0, double lambda)
{
    double s = 0;
    for(int k = 0; k <= k0 ; k++)
    {
        s += poisson_pmf(k, lambda);
    }
    // printf("poisscdf(%e, %e) == %e\n", k0, lambda, s);
    return s;
}

/* Get the value of the linspace from a to b with n points at idx */
static double linspace(double a, double b, size_t n, size_t idx)
{
    return a + (b-a)*(double) idx/ (double) (n-1);
}

/* Return the smallest number of points, n, required for a
 * linspace [a, b] with dx being at most delta */
static size_t linspace_n(double a, double b, double delta)
{
    return (size_t) ceil((b-a) / delta) + 1;
}

static void linspace_ut()
{

    for(int kk = 1; kk<5; kk+=kk)
    {
        /* Middle point should be 0 */
        //printf("%f, %f, %zu, %zu, %f\n", -1.0, 1.0, 2*kk+1, kk, linspace(-1, 1, 2*kk+1, kk));
        assert( linspace(-1, 1, 2*kk+1, kk) == 0 );
        /* 0 -> a */
        assert( linspace(-1, 1, 2*kk+1, 0) == -1 );
        /* n-1 -> b */
        assert( linspace(-1, 1, 2*kk+1, 2*kk) == 1);
    }
    /* [0, 1] */
    assert(linspace_n(0, 1, 1) == 2);
    /* [-1, 0] */
    assert(linspace_n(-1,0, 1) == 2);
    assert(linspace_n(0, 5, 1) == 6);
}

static opts * opts_new()
{
    opts * s = calloc(1, sizeof(opts));
    assert(s != NULL);
    s->verbose = 1;
    s->capture_distance = 4;
    s->sigma = 0.4;
    s->mag1 = 1.0;
    s->mag2 = 1.0;
    s->npoint = 250;
    s->weightz = 0.5;
    return s;
}

static void opts_free(opts * s)
{
    free(s->fname_d1);
    free(s->fname_d2);
    free(s->outfile);
    free(s);
}

static void opts_print(opts * s, FILE * fid)
{
    fprintf(fid, "\n");
    fprintf(fid, "Settings:\n");
    fprintf(fid, "             verbosity: %d\n", s->verbose);
    fprintf(fid, "             KDE sigma: %f\n", s->sigma);
    fprintf(fid, "  Magnification is set: %s\n", dw_yes_no(s->mag_is_set));
    if(s->mag_is_set)
    {
        fprintf(fid, "                  mag1: %f\n", s->mag1);
        fprintf(fid, "                  mag2: %f\n", s->mag2);
    }
    fprintf(fid, "               weightz: %f\n", s->weightz);
    fprintf(fid, "                  rotz: %f\n", s->rotz);
    fprintf(fid, "             overwrite: %s\n", dw_yes_no(s->overwrite));
    fprintf(fid, "                append: %s\n", dw_yes_no(s->append));
    fprintf(fid, "      capture distance: %f\n", s->capture_distance);
    fprintf(fid, "\n");
    return;
}

static void load_delta(const char * fname, double * delta)
{
    npio_t * np = npio_load(fname);
    if(np == NULL)
    {
        fprintf(stderr, "Unable to load %s\n", fname);
        exit(EXIT_FAILURE);
    }
    if(np->nel != 3)
    {
        fprintf(stderr, "Wrong number of elements in %s\n", fname);
        exit(EXIT_FAILURE);
    }

    if(np->dtype == NPIO_F32)
    {
        float * X = (float*) np->data;
        for(int kk = 0; kk < 3; kk++)
        {
            delta[kk] = X[kk];
        }
        npio_free(np);
        return;
    }

    if(np->dtype == NPIO_F64)
    {
        double * X = (double*) np->data;
        for(int kk = 0; kk < 3; kk++)
        {
            delta[kk] = X[kk];
        }
        npio_free(np);
        return;
    }

    fprintf(stderr,
            "Unable to load data from %s\n"
            "Please check the data type and the number of elements\n",
            fname);
    exit(EXIT_FAILURE);
}

static void usage(void)
{
    opts * dopts = opts_new();
    printf(
           "Estimates a displacement vector d=[dx, dy, dz] so that\n"
           "A + d overlaps B, where A are dots from file1.tsv and\n"
           "B are dots from file2.tsv\n"
           "The input files should be generated with dw dots\n");
    printf("\n");
    printf("usage: dw align-dots [<options>] file1.tsv file2.tsv\n");
    printf("\n");
    printf("Options:\n");
    printf("  --out file, -o file\n"
           "\tWhere to write the results.\n");
    printf("  --overwrite\n"
           "\tOverwrite the destination file if it exists\n");
    printf("  --append, -a\n"
           "\tAppend to the output file\n");
    printf("  --sigma s\n"
           "\tSet the size of the KDE used to identify the peak\n"
           "\tShould be set approximately to the localization accurracy in the lateral plane\n"
           "\tDefault value: %.2f\n", dopts->sigma);
    printf("  --npoint n\n"
           "\tset the maximum number of points to use from each file"
           "\tDefault value: %zu\n", dopts->npoint);
    printf("  --mag1 f\n"
           "\tMultiplicative magnification factor for the 1st point set in the lateral plane\n"
           "\tDots coordinates will be multiplied with this factor directly after loading\n"
           "\tUnless either --mag1 and/or --mag2 is set the algorithm will be more restrictive\n"
           "\tlooking for correspondences\n");
    printf("  --mag2 f\n"
           "\tMultiplicative magnification factor for the 2nd point set\n");
    printf("  --rotz d\n"
           "\tRotation around the z-axis at (x,y)=(0,0), hence\n"
           "\tonly small rotations like +/- 0.0001 are meaningful\n");
    printf(" --weightz w\n"
           "\tSet the weight/scaling of the z-component of the dot coordinates\n"
           "\tDefault value: %.2f\n"
           "\tSet so this so that the localization accuracy in the lateral and axial planes conincide\n",
           dopts->weightz);
    printf("\n");
    printf("Depreciated options, that only applies to the fallback/brute force algorithm:\n");
    printf("  --radius d, -d d\n"
           "\tSet the capture radius, i.e. largest expected shift in pixels\n"
           "\tDefault value: %.1f\n", dopts->capture_distance);
    printf("  --delta1 delta1.npy\n"
           "\tInitial displacements for dots in file1\n");
    printf("  --delta2 delta2.npy\n"
           "\tInitial displacements for dots in file2\n");
    printf("\n"
           "Usage notes:\n"
           "\n"
           "If no correpondences are found, try increasing --sigma and/or --ndots\n");
    printf("\n");
    printf("Output columns, written to the --out file:\n");
    printf(" 1. dx\n");
    printf(" 2. dy\n");
    printf(" 3. dz\n");
    printf(" 4. kde, approximately the number of correspondences found\n");
    printf(" 5. goodness or confidence. 1=perfect\n");
    printf(" 6. file1.tsv \n");
    printf(" 7. file2.tsv\n");
    printf(" 8. magnification applied to the 1st set\n");
    printf(" 9. magnification for 2nd set\n");
    printf("10. sigma\n");
    printf("11. search radius (only relevant for the fallback algorithm)\n");
    printf("\n");
    opts_free(dopts);
    dopts = NULL;
    return;}


static void argparsing(int argc, char ** argv, opts * s)
{

    struct option longopts[] = {
        {"help", no_argument, NULL, 'h'},
        {"verbose", required_argument, NULL, 'v'},
        {"radius",   required_argument, NULL, 'd'},
        {"mag1",      required_argument, NULL, 'm'},
        {"mag2",      required_argument, NULL, 'M'},
        {"npoint", required_argument, NULL, 'n'},
        {"out",    required_argument, NULL, 'o'},
        {"append", no_argument, NULL, 'a'},
        {"rotz", required_argument, NULL, 'r'},
        {"sigma", required_argument, NULL, 's'},
        {"overwrite", no_argument, NULL, 'w'},
        {"weightz", required_argument, NULL, 'z'},
        {"delta1", required_argument, NULL, '1'},
        {"delta2", required_argument, NULL, '2'},
        {NULL, 0, NULL, 0}};

    int ch;
    while( (ch = getopt_long(argc, argv,
                             "ad:hm:M:n:o:r:s:v:wz:1:2:",
                             longopts, NULL)) != -1)
    {
        switch(ch)
        {
        case '1':
            free(s->fname_d1);
            s->fname_d1 = strdup(optarg);
            break;
        case '2':
            free(s->fname_d2);
            s->fname_d2 = strdup(optarg);
            break;
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
        case 'm':
            s->mag1 = atof(optarg);
            s->mag_is_set = 1;
            break;
        case 'M':
            s->mag2 = atof(optarg);
            s->mag_is_set = 1;
            break;
        case 'n':
            s->npoint = atol(optarg);
            break;
        case 'o':
            free(s->outfile);
            s->outfile = strdup(optarg);
            break;
        case 'r':
            s->rotz = atof(optarg);
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
        case 'z':
            s->weightz = atof(optarg);
            break;
        default:
            exit(EXIT_FAILURE);
        }
    }

    if(s->outfile != NULL)
    {
        if(s->append == 0 && s->overwrite == 0)
        {
            if(dw_isfile(s->outfile))
            {
                fprintf(stderr, "Error: The output file does already exists\n");
                fprintf(stderr, "Consider --overwrite or --append\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    s->optind = optind;

    if(s->fname_d1)
    {
        load_delta(s->fname_d1, s->delta1);
        if(s->verbose > 1)
        {
            printf("delta1=[%f, %f, %f]\n", s->delta1[0], s->delta1[1], s->delta1[2]);
        }
    }
    if(s->fname_d2)
    {
        load_delta(s->fname_d2, s->delta2);
        if(s->verbose > 1)
        {
            printf("delta2=[%f, %f, %f]\n", s->delta2[0], s->delta2[1], s->delta2[2]);
        }
    }
    return;
}

/* Multiply the x and y coordinates in X by mag */
static void
magnify_dots(double * X, size_t n, double mag)
{
    if(mag == 1.0)
    {
        return;
    }

    for(size_t kk = 0; kk<n ; kk++)
    {
        for(size_t dd = 0; dd < 2; dd++)
        {
            X[3*kk + dd] *= mag;
        }
    }
    return;
}

static void scale_z(double * X, i64 n, double weightz)
{
    if(weightz == 1.0)
    {
        return;
    }
    for(i64 kk = 0; kk < n; kk++)
    {
        X[3*kk + 2] *= weightz;
    }
    return;
}

static void rotz_dots(double * X, size_t n, double rot)
{
    if(rot == 0)
    {
        return;
    }
    double vcos = cos(rot);
    double vsin = sin(rot);
    for(size_t kk = 0; kk<n; kk++)
    {
        double x = X[3*kk];
        double y = X[3*kk+1];
        X[3*kk] = vcos*x - vsin*y;
        X[3*kk+1] = vsin*x + vcos*y;
    }
}

/* Extract the X,Y,Z values from the table */
static double * get_dots(opts * s, ftab_t * T)
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

    i64 nuse = T->nrow;

    int cv = ftab_get_col(T, "value");
    if(cv >= 0)
    {
        nuse = 0;
        for(i64 kk = 0; kk < (i64) T->nrow; kk++)
        {
            float * row = T->T + kk*T->ncol;
            float value = row[cv];
            if(value > 0)
            {
                nuse++;
            } else {
                break;
            }
        }
    }

    if(s->verbose > 1)
    {
        if(cv >= 0)
        {
            printf("Will load %ld / %ld dots where value > 0\n",
                   nuse, T->nrow);
        } else {
            printf("Will load %ld dots\n", nuse);
        }
    }

    double * X = calloc(3*nuse, sizeof(double));
    assert(X != NULL);
    for(i64 kk = 0; kk < nuse; kk++)
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

static double v3_norm_sq(const double * X)
{
    return pow(X[0], 2) + pow(X[1], 2) + pow(X[2], 2);
}

static double eudist3_sq(const double * A, const double * B)
{
    return pow(A[0]-B[0], 2.0) + pow(A[1]-B[1], 2.0) + pow(A[2]-B[2], 2.0);
}

static void align_dots_bf(opts * s,
                          const double * XA, size_t nXA,
                          const double * XB, size_t nXB)
{
    /* Brute Force (BF) alignment of dots from the two sets
     */

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

    if(nfound_total < 1)
    {
        printf("No points to work with, exiting\n");
        exit(EXIT_FAILURE);
    }

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
    double rs2 = pow(rs, 2.0);

    size_t n = linspace_n(-rs, rs, 0.5);

    linspace_ut();

    double * KDE = calloc(n*n*n, sizeof(double));
    assert(KDE != NULL);

    /* Initial grid search */
#pragma omp parallel for
    for(size_t iz = 0; iz < n; iz++) {
        double z = linspace(-rs, rs, n, iz);
        for(size_t iy = 0; iy < n; iy ++) {
            double y = linspace(-rs, rs, n, iy);
            for(size_t ix = 0; ix < n; ix++) {
                double x = linspace(-rs, rs, n, ix);

                double P[] = {x, y, z};
                if(v3_norm_sq(P) > rs2){
                    continue;
                }
                double v = kdtree_kde(TD, P, s->sigma, 0);
                KDE[n*n*iz + n*iy + ix] = v;
            }
        }
    }

    /* Search through the grid to determine the max pos and value */
    for(size_t iz = 0; iz < n; iz++) {
        for(size_t iy = 0; iy < n; iy ++) {
            for(size_t ix = 0; ix < n; ix++) {
                double v = KDE[n*n*iz + n*iy + ix];
                if( v > maxkde )
                {
                    double x = linspace(-rs, rs, n, ix);
                    double y = linspace(-rs, rs, n, iy);
                    double z = linspace(-rs, rs, n, iz);
                    maxkde = v;
                    maxpos[0] = x; maxpos[1] = y; maxpos[2] = z;
                }
            }
        }
    }

    if(s->verbose > 1)
    {
        printf("     Grid search: (% .2f, % .2f, % .2f) (kde=%.1f)\n",
               maxpos[0], maxpos[1], maxpos[2], maxkde);
    }

    /* Refinement over the grid search */
    rs = 1.0*s->sigma; // Region size
    rs2 = pow(rs, 2.0);
    while(rs > 1e-4)
    {
        double center[3];
        memcpy(center, maxpos, 3*sizeof(double));
        for(double x = -rs; x <= rs; x += rs/7.0) {
            for(double y = -rs; y <= rs; y += rs/7.0) {
                for(double z = -rs; z <= rs; z += rs/7.0) {
                    // Only consider inside a sphere
                    double X[] = {x, y, z};
                    if(v3_norm_sq(X) > rs2){
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


    /* 2nd grid search to figure out what the 2nd best position is */
    double maxkde2 = 0;
    double maxpos2[3] = {0};
    {
        /* We will exclude points withing 1 + s->sigma */
        double sigma2 = pow(2.0 + s->sigma, 2.0);
        if(s->verbose > 1)
        {
            printf("Grid search #2 -- searching for 2nd maxima\n");
        }
        double rs = s->capture_distance;
        double rs2 = pow(rs, 2.0);


        //for(double z = -rs; z <= rs; z+=0.5) {
        n = linspace_n(-rs, rs, 0.5);
        // double * KDE = calloc(n*n*n, sizeof(double)); // TODO

#pragma omp parallel for
        for(size_t iz = 0; iz < n; iz++) {
            double z = linspace(-rs, rs, n, iz);
            for(size_t iy = 0; iy < n; iy ++) {
                double y = linspace(-rs, rs, n, iy);
                for(size_t ix = 0; ix < n; ix++) {
                    double x = linspace(-rs, rs, n, ix);

                    double P[] = {x, y, z};
                    if(v3_norm_sq(P) > rs2){
                        continue;
                    }
                    /* Don't compare close the maxima that we found */
                    if(eudist3_sq(P, maxpos) < sigma2)
                    {
                        continue;
                    }
                    double v = KDE[n*n*iz + n*iy + ix];
                    if(v > maxkde2)
                    {
                        maxkde2 = v;
                        memcpy(maxpos2, P, 3*sizeof(double));
                    }
                }
            }
        }
    }

    free(KDE);
    KDE = NULL;

    double P[] = {0,0,0};
    double origo_kde = kdtree_kde(TD, P, s->sigma, 0);
    if(s->verbose > 1)
    {
        printf("At origo, kde=%.2f (without any correction)\n", origo_kde);
    }
    if(s->verbose > 1)
    {
        printf("Refined position: (% .2f, % .2f, % .2f) (kde=%.1f)\n",
               maxpos[0],maxpos[1], maxpos[2], maxkde);
    }

    if(maxkde < origo_kde)
    {
        fprintf(stderr, "Error: unexpected result! Please report this as a bug\n");
        fprintf(stderr, "       at %s %d\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
    }

    if(s->verbose > 1)
    {
        printf("secondary maxima value: %f at (%f, %f, %f)\n",
               maxkde2,
               maxpos2[0],
               maxpos2[1],
               maxpos2[2]);
    }

    s->dx = maxpos[0];
    s->dy = maxpos[1];
    s->dz = maxpos[2];

    /* Please note that these values will be overwritten by the
       function determine_alignment_quality */

    s->kde = maxkde;
    /* Gap / combined height */
    s->goodness = (maxkde-maxkde2) / (maxkde + maxkde2);

    if(s->verbose > 1)
    {
        printf("Estimated shift: [% .2f, % .2f, % .2f] pixels, G=%.2f, KDE=%.2f, 2nd KDE: %.2f mag1=%f mag2=%f\n",
               s->goodness,
               maxpos[0], maxpos[1], maxpos[2],
               maxkde, maxkde2, s->mag1, s->mag2);

        if(maxkde < 2)
        {
            printf("Warning: The KDE value is very low, should be at least 2\n");
        }

        if(maxkde < 2.0 * maxkde2)
        {
            printf("Warning: The KDE should be at least 2X of the secondary peak\n");
        }
    }


    kdtree_free(TD);
    TD = NULL;

    return;

}


/* Add a D=[dx, dy, dz] to each point in X */
static void shift_dots(double * X, double * D, i64 n)
{
    double s = fabs(D[0]) + fabs(D[1]) + fabs(D[2]);
    if(s == 0)
    {
        return;
    }

    for(i64 kk = 0; kk < n; kk++)
    {
        X[3*kk + 0] = X[3*kk + 0] + D[0];
        X[3*kk + 1] = X[3*kk + 1] + D[1];
        X[3*kk + 2] = X[3*kk + 2] + D[2];
    }
}


static void
get_displacement_from_qalign_result(opts * s,
                                    const double * qD, // 3xnqD vector
                                    i64 nqD,
                                    double * Q1) // return value
{
    /* Use mean shift to find a single displacement vector from the
       hypotheses returned from qalign */

    if(s->verbose > 1)
    {
        printf("Refining the result among %ld points\n", nqD);
    }
    kdtree_t * T = kdtree_new(qD, nqD, 20);
    i64 start_idx = 0;
    double best_kde = -1;
    for(i64 kk = 0; kk < nqD; kk++)
    {
        //printf("%f, %f, %f\n", qD[3*kk], qD[3*kk+1], qD[3*kk+2]);
        double kde = kdtree_kde(T, qD+3*kk, s->sigma, -1);
        if(kde > best_kde)
        {
            best_kde = kde;
            start_idx = kk;
        }
    }
    if(s->verbose > 1)
    {
        printf("Best start point at %ld (kde=%f)\n", start_idx, best_kde);
    }

    double Q0[3] = {qD[start_idx*3 + 0],
                    qD[start_idx*3 + 1],
                    qD[start_idx*3 + 2]};

    for(i64 kk = 0; kk < 100; kk++)
    {
        kdtree_kde_mean(T,
                        Q0,
                        s->sigma, // radius
                        -1, // cutoff (auto)
                        Q1);
        memcpy(Q0, Q1, 3*sizeof(double));
    }

    s->kde = kdtree_kde(T, Q1, s->sigma, -1);

    if(s->verbose > 1)
    {
        printf("qalign -> [%f, %f, %f] (kde=%.2f)\n", Q1[0], Q1[1], Q1[2], s->kde);
    }

    kdtree_free(T);
    return;
}

/* Volume of a point cloud calculated via the bounding box */
static double
bbx_volume(const double * X, i64 nX)
{
    double mincoord[3];
    double maxcoord[3];
    for(int kk = 0; kk < 3; kk++)
    {
        mincoord[kk] = X[kk];
        maxcoord[kk] = X[kk];
    }
    for(i64 kk = 0; kk < nX; kk++)
    {
        for(int ll = 0; ll < 3; ll++)
        {
            double c = X[3*kk+ll];
            c > maxcoord[ll] ? maxcoord[ll] = c : 0;
            c < mincoord[ll] ? mincoord[ll] = c : 0;
        }
    }
    return (maxcoord[0]-mincoord[0])
        *(maxcoord[1]-mincoord[1])
        *(maxcoord[2]-mincoord[2]);
}


static int determine_alignment_quality(opts * s,
                                       const double * XA, i64 nXA,
                                       const double * XB, i64 nXB)
{
    /* Determine the number of times that points overlap calculate the
     * probability for assuming that the dots are uniformly randomly
     * distributed and placed in bins so that the number of dots per
     * bin is distributed according to Poisson.  Finally calculate a
     * goodness score in [0, 1] such that a value of 1 indicates that
     * this could not happen if XA and XB were uncorrelated.
     */

    /* Number of times that a point from B overlaps a point in A after
       shift correction */
    double overlap_threshold = 0.7;
    kdtree_t * T = kdtree_new(XA, nXA, 20);
    i64 nalign = 0;
    double kde = 0;
    for(i64 kk = 0; kk < nXB; kk++)
    {
        double Q[3];
        for(int ii = 0; ii < 3; ii++)
        {
            Q[ii] = XB[3*kk + ii] - s->delta[ii];
        }
        double kde0 = kdtree_kde(T, Q, s->sigma, -1);
        kde += kde0;
        if(kde0 > overlap_threshold)
        {
            nalign++;
        }
    }
    printf("%ld correspondences were found after alignment\n", nalign);
    kdtree_free(T);
    s->kde = kde;


    /* If there was a 3D grid over the images, how many bins would it have?
     * We assume that the images have the same size. Would be more accurate
     * to calculate the bbx from all points directly */
    double ncell = bbx_volume(XA, nXA) / pow(s->sigma, 3);
    {
        double ncell2 = bbx_volume(XB, nXB) / pow(s->sigma, 3);
        ncell2 > ncell ? ncell = ncell2 : 0;
    }


    if(ncell < 1e-6 )
    {
        fprintf(stderr,
                "Failed to calculate the goodness value since the dot coordinates "
                "in combination with sigma does not make sense\n");
        s->goodness = -1;
        return EXIT_FAILURE;
    }

    /* We substract 1 since there will always be one aligned dot when
       successful and that should not count */
    double observed = nalign-1.0;
    if(observed < 1.0)
    {
        s->goodness = 0;
        return EXIT_FAILURE;
    }

    /* If the points in XA and XB are uniformly random and uncorrelated
     * how many times do we expect two dots in the same bin? */
    double lambda1 = nXA / ncell;
    double lambda2 = nXB / ncell;
    double prob1 = 1.0 - poisson_cdf(0, lambda1); // == 1 - poisson_pmf(0, lambda)
    double prob2 = 1.0 - poisson_cdf(0, lambda2);
    double expected = ncell*prob1*prob2;

    if(s->verbose > 1)
    {
        printf("ncell: %f\n", ncell);
        printf("lambda1: %e, lambda2: %e\n", lambda1, lambda2);
        printf("prob1: %e, prob2: %e\n", prob1, prob2);
        printf("observed: %f, expected %f\n", observed, expected);
    }

    s->goodness = 1.0 - expected / observed;
    s->goodness > 1.0 ? s->goodness = 1 : 0;
    s->goodness < 0.0 ? s->goodness = 0 : 0;

    return EXIT_SUCCESS;
}


static int
run_qalign(opts * s,
           const double * XA,
           i64 nXA,
           const double * XB,
           i64 nXB)
{
    if( (nXA < 3) || (nXB < 3) )
    {
        printf("Too few dots for qalign (%ld in set A, %ld in set B)\n", nXA, nXB);
        return EXIT_FAILURE;
    }

    struct timespec t0, t1;
    dw_gettime(&t0);
    float * fXA = malloc(3*nXA*sizeof(float));
    for(i64 kk = 0; kk < nXA; kk++)
    {
        fXA[3*kk] = XA[3*kk];
        fXA[3*kk+1] = XA[3*kk+1];
        fXA[3*kk+2] = XA[3*kk+2];
    }
    float * fXB = malloc(3*nXB*sizeof(float));
    for(i64 kk = 0; kk < nXB; kk++)
    {
        fXB[3*kk] = XB[3*kk];
        fXB[3*kk+1] = XB[3*kk+1];
        fXB[3*kk+2] = XB[3*kk+2];
    }

    double * qD = NULL;
    i64 nqD = 0;

    qalign_config * qconf = qalign_config_new();
    qconf->A = fXA;
    qconf->nA = nXA;
    qconf->B = fXB;
    qconf->nB = nXB;
    qconf->localication_sigma = s->sigma;
    qconf->verbose = s->verbose;
    if(s->mag_is_set)
    {
        qconf->rel_error = 0;
    }

    int res = qalign(qconf);
    free(fXA);
    free(fXB);

    if(res != EXIT_SUCCESS)
    {
        if(s->verbose > 1)
        {
            printf("qalign failed\n");
        }
        qalign_config_free(qconf);
        return EXIT_FAILURE;
    }


    float scaling = qselect_f32(qconf->S, qconf->nH, qconf->nH/2);
#if 0
    float scalingZ = qselect_f32(qconf->SZ, qconf->nH, qconf->nH/2);
    printf("S: %f, SZ: %f\n", scaling, scalingZ);
#endif

    if(scaling > (1.0+1e-4) || (1.0 / scaling) > (1.0+1e-4))
    {
        if(s->verbose > 0)
        {
            printf("\n");
            printf("!!! Lateral magnification mismatch: "
                   "median(mag A/ mag B) = %f\n", scaling);
            printf("    Either scale the input data and set --mag1 0\n"
                   "    or let use --mag1 %f OR --mag2 %f\n", 1.0/scaling, scaling);
            printf("\n");
        }
    }

    nqD = qconf->nH;
    qD = malloc(nqD*3*sizeof(double));
    for(i64 kk = 0; kk < nqD; kk++)
    {
        for(i64 ii = 0; ii < 3; ii++)
        {
            qD[3*kk + ii] = qconf->H[kk].delta[ii];
        }
    }


    qalign_config_free(qconf);
    dw_gettime(&t1);
    if(s->verbose > 1)
    {
        printf("qalign took %f s\n", timespec_diff(&t1, &t0));
    }

    if(nqD < 2)
    {
        if(s->verbose > 0)
        {
            printf("qalign could not find any correspondences\n");
        }
        free(qD);
        return EXIT_FAILURE;
    }

    double Q1[3] = {0};

    get_displacement_from_qalign_result(s,
                                        qD, nqD,
                                        Q1);
    free(qD);
    qD = NULL;
    for(int i = 0; i < 3; i++)
    {
        s->delta[i] = -Q1[i];
    }

    if(s->kde < 2)
    {
        if(s->verbose > 0)
        {
            printf("qalign could not find anything\n");
            return EXIT_FAILURE;
        }
    }

    return EXIT_SUCCESS;
}

static void print_results(FILE * fid, opts * s,
                          const char * file1,
                          const char * file2)
{

    fprintf(fid, "%f, %f, %f", s->dx, s->dy, s->dz);
    fprintf(fid, ", %f, %f", s->kde, s->goodness);
    fprintf(fid, ", '%s'", file1);
    fprintf(fid, ", '%s'", file2);
    fprintf(fid, ", %f", s->mag1);
    fprintf(fid, ", %f", s->mag2);
    fprintf(fid, ", %f", s->sigma);
    fprintf(fid, ", %f", s->capture_distance);
    fprintf(fid, "\n");
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

    if(!dw_isfile(argv[optind]))
    {
        fprintf(stderr, "%s does not exist or can not be opened\n",
                argv[optind]);
        opts_free(s);
        exit(EXIT_FAILURE);
    }

    if(!dw_isfile(argv[optind+1]))
    {
        fprintf(stderr, "%s does not exist or can not be opened\n",
                argv[optind+1]);
        opts_free(s);
        exit(EXIT_FAILURE);
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

    magnify_dots(XA, nXA, s->mag1);
    magnify_dots(XB, nXB, s->mag2);
    scale_z(XA, nXA, s->weightz);
    scale_z(XB, nXB, s->weightz);

    shift_dots(XA, s->delta1, nXA);
    shift_dots(XB, s->delta2, nXB);

    rotz_dots(XB, nXB, s->rotz);

    if(s->verbose > 1)
    {
        printf("Will align %zu dots vs %zu\n", nXA, nXB);
    }

    /* Only try the brute force algorithm if the quick does not find
       anything. It seems like the brute force is never better so it
       should be removed eventually. */
    if(run_qalign(s, XA, nXA, XB, nXB) != EXIT_SUCCESS)
    {
        if(s->verbose > 0)
        {
            printf("Trying the fallback algorithm\n");
        }
        align_dots_bf(s, XA, nXA, XB, nXB);
    }

    /* evaluate s->delta to set s->goodness and s->kde
     */
    determine_alignment_quality(s, XA, nXA, XB, nXB);

    /* Here we could make a higher order estimation from correspondence points */

    /* Scale back the z-component */
    s->dz /= s->weightz;

    free(XA);
    free(XB);

    /* Present the result */
    if(s->outfile != NULL)
    {
        FILE * fid;
        if(s->append)
        {
            fid = fopen(s->outfile, "a");
        } else {
            fid = fopen(s->outfile, "w");
        }
        print_results(fid, s, argv[optind], argv[optind+1]);
        fclose(fid);
    }
    if(s->verbose > 0)
    {
        printf("Shift: [%.2f, %.2f, %.2f]", s->dx, s->dy, s->dz);
        printf(", goodness = %.2e ", s->goodness);
        if(s->goodness > 0.999)
        {
            printf("(high confidence)");
        } else {
            printf("(low confidence)");
        }
        printf("\n");
    }

    opts_free(s);
    s = NULL;
    return EXIT_SUCCESS;
}

#ifdef STANDALONE
int main(int argc, char ** argv)
{
    return dw_align_dots(argc, argv);
}
#endif
