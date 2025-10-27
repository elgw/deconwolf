#pragma once

/* A space partitioning K-d tree, like described on
 * https://en.wikipedia.org/wiki/K-d_tree
 *
 * This particular implementation is hard coded for 3D data. You can
 * of course use 1D and 2D data as well by setting the coordinates of
 * the unused dimensions to 0.
 *
 * Usage:
 * Acquire a kdtree_t* by  kdtree_new and free it with kdtree_free.
 *
 * link with -kdtree -lm
 */


#include <stdint.h>
#include <stddef.h>

#define KDTREE_VERSION_MAJOR 0
#define KDTREE_VERSION_MINOR 1
#define KDTREE_VERSION_PATCH 0

/* Will only work when KDTREE_DIM=3 unless modified */
#define KDTREE_DIM 3

#ifdef WIN32
#define PUB __declspec(dllexport)
#else
#define PUB __attribute__((visibility("default")))
#endif

struct pqheap;

/* A tree or a leaf */
typedef struct {
    size_t id; // number; // node number TODO not needed

    /* Bounding box, [minx, maxx, miny, maxy, minz, maxz] */
    double bbx[2*KDTREE_DIM];

    /* Tells where in XID that the data for this node can be found */
    size_t data_idx;

    /* The number of points associated with the node
     * which are found in T->XID */
    uint32_t n_points;

    /* If not a leaf, this tells how this node was split */
    uint8_t split_dim; /* is set to KDTREE_DIM for leafs */
    double pivot;
} kdtree_node_t;

typedef struct{
    /* Node allocation  */
    kdtree_node_t * nodes; /* Array of nodes, nodes[0] is the root */
    size_t n_nodes_alloc; // Total number of nodes

    /* Storage for coordinates and their indexes, [(x, y, z), id] */
    double * XID;

    /* Maximum number of points per leaf (i.e. end node) */
    size_t max_leaf_size;
    size_t n_points; // Number of supplied points

    /* Temporary buffer used during tree construction */
    double * median_buffer;

    /* State variables for queries */
    struct pqheap * pq; // needs to be here?
    int direct_path;
    /* The latest query is stored internally to avoid an abundant
       number of malloc/free. Can of course be copied by the caller. */
    size_t * result; // KN for storing idx of K neighbours
    size_t result_alloc; /* number of elements allocated for result */
} kdtree_t;

/* Construct a new tree based on the N points stored in X
 *
 * binsize (or max_leaf_size) is an algorithmic parameter. According
 * to [1] bin size a of 4-32 elements is optimal regardless of the
 * number of dimensions
 */
PUB kdtree_t *
kdtree_new(const double * X,
           size_t N, int binsize);


/* Frees all resources associated with a tree */
PUB void kdtree_free(kdtree_t * T);

/* Query one point for its k nearest neighbours.  The returned array
 * contains the index of k points, sorted according to the distance of
 * the points, with the closest point first.
 *
 * Important: The returned array is owned by the tree and should not
 * be freed. It will be re-used with the next call to kdtree_query_*
 */
PUB size_t *
kdtree_query_knn(kdtree_t * T,
                 const double * Q,
                 size_t k);

/* TODO
 * Find all points within some radius of Q
 * Returns a newly allocated array of indexes of length nfound
 * On failure: Returns NULL and sets nfound to 0
 */
PUB size_t *
kdtree_query_radius(const kdtree_t * T,
                    const double * Q,
                    const double radius,
                    size_t * nfound);

/* Estimate the local density using non-normalized Gaussian symmetric
 * kernel with a fixed sigma.
 *
 * G(x, sigma) = exp(-x^2 / (2*sigma^2))
 *
 * For a custom kernel please use kdtree_query_radius and calculate
 * the KDE based on the found points.
 *
 * Cutoff: will use points up to sigma*cutoff away from the query
 * point Q. A default value will be used if cutoff == -1.
 */
PUB double
kdtree_kde(const kdtree_t * T,
           const double * Q,
           double sigma,
           double cutoff);

/* Calculate the weighted mean position of a Gaussian KDE at the point
 * Q as the following sum over all the neighbors
 *
 * m(x) = \frac{ \sum K(x_i-x) x_i }{\sum K(x_i -x)}
 *
 * returns: m -- the weighted mean position. Will be Q if no neighbours found.
 *
 * This is an ingredient of the mean shift algorithm. Please not that
 * using a Gaussian kernel a quite large radius contributes to the
 * kde, set cutoff to you liking. A cutoff of 0 means that it will be
 * set automagically.
 */

PUB void
kdtree_kde_mean(const kdtree_t *,
                const double * Q,
                double sigma,
                double cutoff,
                double * mean);

/* Wanted: Expectation Maximization (EM) with Gaussian Mixture Model (GMM)
void kdtree_emgmm(const kdtree_t * T,
const gaussian ** G0,
gaussian ** Gfinal);
 */



/* Find the index of the closest point */
PUB size_t kdtree_query_closest(kdtree_t * T, double * X);

PUB void node_print_bbx(const kdtree_node_t * N);

/* Make a shallow copy of a kd-tree for usage by another thread */
PUB kdtree_t * kdtree_copy_shallow(kdtree_t * );

/* Free a tree returned from kdtree_copy_shallow */
PUB void kdtree_free_shallow(kdtree_t * T);


/* Run some self-tests */
PUB void kdtree_validate(kdtree_t * T);

PUB void kdtree_print_info(kdtree_t * T);
