#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#ifdef GSL
#include <gsl/gsl_statistics_double.h>
#endif

#include "kdtree.h"
#include "pqheap.h"
#include "quickselect.h"

#define XID_STRIDE (KDTREE_DIM + 1)



/* Resolve the index of the right child based on the index of the
 * parent, following a Eytzinger scheme */
static size_t node_right_child_id(size_t node_id)
{
    return 2*node_id+2;
}

static size_t node_left_child_id(size_t node_id)
{
    return 2*node_id+1;
}


static void node_set_final(kdtree_node_t * node)
{
    node->split_dim = KDTREE_DIM;
    return;
}

static int node_is_final(const kdtree_node_t * node)
{
    return (node->split_dim == KDTREE_DIM);
}

/*  Hoare's partition scheme for vectors where the partitioning is
 *  performed based on a single dimension or index.  Adopted from
 *  arch/24/03/11_quickselect
 *
 * X : the data of size [XID_STRIDE x n] which will be partitioned.
 *
 * vim: the index or dimension of the pivot
 *
 * Returns:
 * nLow: the number of vectors where v[vdim] <= pivot
 * nHigh: the number of vectors where v[vdim] > pivot
 */
static void
partition_vectors(double * restrict X,
          const size_t n, /* Number of points */
          const size_t vdim, /* Dimension to take value from */
          const double pivot,
          size_t * nLow, size_t * nHigh)
{
    int64_t low = -1;
    int64_t high = n;
    int64_t n2 = n;

    while(1)
    {
        do { low++; } while ( (low < n2)  && X[low*XID_STRIDE+vdim] <= pivot );

        do { high--; } while ( (high > 0) && X[high*XID_STRIDE+vdim] > pivot );

        if(low >= high)
        { *nLow = low;  *nHigh = n-*nLow;
#ifndef NDEBUG
            assert(*nLow + *nHigh == n );
            for(int64_t kk = 0; kk < low; kk++)
            {
                //printf("Pivot = %f\n", pivot);
                //print_XID(X, low);
                assert(X[kk*XID_STRIDE+vdim] <= pivot);
            }
            for(int64_t kk = low; kk < n2; kk++)
            {
                assert(X[kk*XID_STRIDE+vdim] > pivot);
            }
#endif
            return;
        }

        /* Swap data */
        double t[XID_STRIDE];
        memcpy(t,
               X + low*XID_STRIDE,
               XID_STRIDE*sizeof(double));
        memcpy(X + low*XID_STRIDE,
               X + high*XID_STRIDE,
               XID_STRIDE*sizeof(double));
        memcpy(X + high*(KDTREE_DIM+1),
               t,
               XID_STRIDE*sizeof(double));
    }

    return;
}


void
node_print_bbx(const kdtree_node_t * N)
{

    printf("#%zu: ", N->id);
    for(size_t kk = 0; kk < KDTREE_DIM; kk++)
    {
        printf("[%f -- %f]", N->bbx[2*kk], N->bbx[2*kk+1]);
        if(kk + 1 < KDTREE_DIM)
        {
            printf("x");
        }
    }
    printf("\n");
    return;
}

void kdtree_free(kdtree_t * T)
{
    if(T == NULL)
    {
        return;
    }
    free(T->XID);
    T->XID = NULL;
    free(T->nodes);
    T->nodes = NULL;
    if(T->pq != NULL)
    {
        pqheap_free(&T->pq);
    }
    free(T->result);
    T->result = NULL;
    free(T);
    return;
}

kdtree_t * kdtree_copy_shallow(kdtree_t * _T)
{
    assert(_T != NULL);
    if(_T == NULL)
        return NULL;
    kdtree_t * T = calloc(1, sizeof(kdtree_t));
    assert(T != NULL);
    memcpy(T, _T, sizeof(kdtree_t));
    T->result = NULL;
    T->result_alloc = 0;
    T->pq = NULL;
    return T;
}

void kdtree_free_shallow(kdtree_t * T)
{
    pqheap_free(&T->pq);
    free(T->result);
    free(T);
    return;
}

/* Euclidean distance squared */
static double eudist_sq(const double * A, const double * B)
{
    double sum = 0;
    for(size_t ii = 0; ii<KDTREE_DIM; ii++)
    {
        sum+=pow(A[ii]-B[ii], 2);
    }
    return sum;
}


static double eudist(const double * A, const double * B)
{
    return sqrt(eudist_sq(A, B));
}


double get_median_from_strided(const double * X, // data
                               size_t N, // number of points
                               double * T, // temp buffer
                               size_t stride) // stride
{
    // T is a temporary buffer, should be N elements large
    // https://www.gnu.org/software/gsl/doc/html/statistics.html
    // quickselect
    assert(stride == KDTREE_DIM + 1);
    for(size_t kk = 0; kk < N; kk++)
    {
        T[kk] = X[stride*kk];
        //printf("(%f) ", T[kk]);
    }
    //printf("\n");
    //printf("N=%zu, N/2=%zu\n", N, N/2);
#ifdef GSL
    double median = gsl_stats_median(T, 1, N);
#else
    double median = quickselect(T, N, N/2);
#endif
    return median;
}


void bounding_box(const double * restrict X,
                  const size_t N, const size_t D,
                  double * restrict bbx)
{
    for(size_t dd = 0 ; dd < D; dd++)
    {
        bbx[2*dd] = X[dd]; // Min along dimension dd
        bbx[2*dd+1] = X[dd]; // Max along dimensions dd
    }
    for(size_t nn = 0; nn < N; nn++)
    {
        for(size_t dd = 0 ; dd < D; dd++)
        {
            X[D*nn + dd] < bbx[2*dd] ? bbx[2*dd] = X[D*nn + dd] : 0;
            X[D*nn + dd] > bbx[2*dd + 1] ? bbx[2*dd + 1] = X[D*nn + dd] : 0;
        }
    }
    return;
}


/* Recursive splitting  */
void
kdtree_split(kdtree_t * T,
             size_t node_id)
{
    kdtree_node_t * node = T->nodes + node_id;
    assert(node->data_idx % XID_STRIDE == 0);

    //node_print_bbx(node);

    /* Possible to append children without running out of nodes? */
    if(2*node_id+2 >= T->n_nodes_alloc) {
        goto final;
    }

    if(node->n_points < (size_t) T->max_leaf_size)
    {
    final: ; // Construct a "final" node without children
        node_set_final(node);
        return;
    }

    /* Decide along which dimension to split */
    size_t split_dim = 0; // dimension or variable to split on
    {
        double max_size = node->bbx[1] - node->bbx[0];
        for(size_t dd = 0; dd < KDTREE_DIM; dd++)
        {
            double t = node->bbx[2*dd+1] - node->bbx[2*dd];
            assert(t >= 0);
            if(t > max_size)
            {
                split_dim = dd;
                max_size = t;
            }
        }
    }
    node->split_dim = split_dim;

    double * XID = T->XID + node->data_idx;
    assert(node->data_idx + node->n_points <= T->n_points*XID_STRIDE);
    double pivot =
        get_median_from_strided(XID+split_dim,
                                node->n_points,
                                T->median_buffer,
                                XID_STRIDE);

    node->pivot = pivot;
    //printf("[%f   (pivot=%f)   %f]\n", node->bbx[2*split_dim], pivot, node->bbx[2*split_dim+1]);
    assert(node->bbx[2*split_dim] <= pivot);
    assert(node->bbx[2*split_dim+1] >= pivot);


    /* Avoid infinite recursion.
       This happens when points are flat in the splitting dimension */
    if(pivot == node->bbx[2*split_dim]
       || pivot == node->bbx[2*split_dim+1])
    {
        goto final;
    }

    /* Partition the data  */
    size_t nLow = 0;
    size_t nHigh = 0;


    assert(XID[KDTREE_DIM] < T->n_points);
    partition_vectors(XID, node->n_points, split_dim, pivot, &nLow, &nHigh);

    {
        size_t left_id = node_left_child_id(node_id);

        assert(left_id < T->n_nodes_alloc);
        kdtree_node_t * node_left = T->nodes+left_id;
        assert(node_left->id == 0);
        node_left->id = left_id;
        memcpy(node_left->bbx, node->bbx, 2*KDTREE_DIM*sizeof(size_t));
        node_left->bbx[2*split_dim + 1] = pivot;
        node_left->n_points = nLow;

        node_left->data_idx = node->data_idx;

        kdtree_split(T, left_id);
    }

    {
        size_t right_id = node_right_child_id(node_id);
        assert(right_id < T->n_nodes_alloc);
        kdtree_node_t * node_right = T->nodes+right_id;
        assert(node_right->id == 0); /* Unused? */
        node_right->id = right_id;
        memcpy(node_right->bbx, node->bbx, 2*KDTREE_DIM*sizeof(size_t));
        node_right->bbx[2*split_dim] = pivot;
        node_right->n_points = nHigh;

        node_right->data_idx = node->data_idx + nLow*XID_STRIDE;

        kdtree_split(T, right_id);
    }

    return;
}

kdtree_t *
kdtree_new(const double * X,
           size_t N,
           int max_leaf_size)
{

    if(max_leaf_size < 1)
    {
        printf("kdtree_new: invalid bin size\n");
        return NULL;
    }

    if(N < 1)
    {
        printf("kdtree_new: At least one data point needed\n");
        return NULL;
    }

#ifndef NDEBUG
    printf("kdtree warning: Not compiled with -DNDEBUG. Performance will be restrained. \n");
#endif

    /* Set up the tree and the basic settings*/
    kdtree_t * T = calloc(1, sizeof(kdtree_t));
    assert( T!= NULL);
    T->max_leaf_size = max_leaf_size;
    T->n_points = N;

    /* Allocate storage for the nodes. We allocate enough
     * nodes for a complete binary tree up to some depth.
     * I.e. we will have 1, 3, 7, 15, ... (2^(L+1)-1) nodes where L
     * is the number of leafs.
     */

    {
        /* If each leaf is 50% full we will have approximately */
        double n_leafs0 = 2.0* (double) N / (double) max_leaf_size;
        /* Since it has to be a power of two we pick */
        double n_leafs = pow(2.0, ceil(log2(n_leafs0)));
        /* Then the number of nodes needed is */
        T->n_nodes_alloc = n_leafs*2-1;
        T->n_nodes_alloc < 3 ? T->n_nodes_alloc = 3 : 0;
    }
    //T->n_nodes_alloc = 2*N;
    T->nodes = calloc(T->n_nodes_alloc, sizeof(kdtree_node_t));
    assert(T->nodes != NULL);
    if(T->nodes == NULL)
    {
        printf("kdtree_new: Memory allocation failed. Tried to allocate for %zu nodes\n"
               "            but couldn't get it from the system\n",
               T->n_nodes_alloc);
        return NULL;
    }

    /* Copy the data points and attach an index */
    T->XID = calloc((KDTREE_DIM+1)*N, sizeof(double));
    assert(T->XID != NULL);
    for(size_t kk = 0; kk<N; kk++)
    {
        memcpy(T->XID+kk*(KDTREE_DIM+1),
               X+kk*(KDTREE_DIM),
               KDTREE_DIM*sizeof(double));
        size_t * ID = (size_t *) T->XID + (KDTREE_DIM+1)*kk + KDTREE_DIM;
        *ID = kk;
    }

    //print_XID(T->XID, T->N);

    // Expand the region for correct ball-within-region checking
    //double dx = xmax-xmin;
    //double dy = ymax-ymin;
    //xmin -= 2*dx; xmax += 2*dx;
    //ymin -= 2*dy; ymax += 2*dy;

    T->median_buffer = calloc(N, sizeof(double));
    assert(T->median_buffer != NULL);

    kdtree_node_t * node = T->nodes;
    bounding_box(X, N, KDTREE_DIM, node->bbx);
    node->n_points = N;
    node->data_idx = 0;

    /* Recursive construction */
    kdtree_split(T, // Tree
                 0); // node_id (location in array)

    free(T->median_buffer);
    T->median_buffer = NULL;

    //print_XID(T->XID, T->N);
    return T;
}

size_t kdtree_query_closest(kdtree_t * T, double * X)
{

    kdtree_node_t * N = T->nodes;
    while( ! node_is_final(N) )
    {
        int split_dim = N->split_dim;
        if(X[split_dim] > N->pivot)
        {
            N = T->nodes + node_right_child_id(N->id);
        } else {
            N = T->nodes+ + node_left_child_id(N->id);
        }
    }

    double * NX = T->XID + N->data_idx;
    double dmin = eudist(X, NX);
    size_t imin = NX[KDTREE_DIM]; // N->idx[0];
    for(size_t kk = 0; kk<N->n_points; kk++)
    {
        double d = eudist(X, NX + kk*(KDTREE_DIM+1));
        if(d < dmin)
        {
            imin = NX[kk*(KDTREE_DIM+1) + KDTREE_DIM]; // ->idx[kk];
            dmin = d;
        }
    }
    assert(dmin < 1e-9);
    return imin;
}

int within_bounds(const kdtree_node_t * node, const double * Q, const double r)
{
    // Return 1 if the disk centered at Q
    // with radius r is inside the node
    // else 0
    const double x = Q[0];
    const double y = Q[1];
    const double z = Q[2];
    if(x + r > node->bbx[1] || x - r < node->bbx[0] ||
       y + r > node->bbx[3] || y - r < node->bbx[2] ||
       z + r > node->bbx[5] || z - r < node->bbx[4])
    {
        return 0;
    } else {
        //   printf("(%f, %f), r=%f is within [%f, %f, %f, %f]\n", x, y, r,
        //     node->xmin, node->xmax, node->ymin, node->ymax);
        return 1;
    }
}

static inline double min(double a, double b)
{
    if(a<b)
    {
        return a;
    }
    return b;
}

static int bounds_overlap_ball_raw(const double * bbx,
                                   const double * Q,
                                   const double r2)
{
    // find nearest x and nearest y
    double nx = 0;
    double ny = 0;
    double nz = 0;
    const double x = Q[0];
    const double y = Q[1];
    const double z = Q[2];

    if( x < bbx[0])
    {
        nx = bbx[0] - x;
    } else if (x > bbx[1])
    {
        nx = x - bbx[1];
    }

    if( y < bbx[2])
    {
        nx = bbx[2] - y;
    } else if (y > bbx[3])
    {
        ny = y - bbx[3];
    }

    if( z < bbx[4])
    {
        nz = bbx[4] - z;
    } else if (z > bbx[5])
    {
        nz = z - bbx[5];
    }

    if(pow(nx, 2) + pow(ny, 2) + pow(nz, 2) < r2)
    {
        return 1;
    }

    return 0;

}

static int
bounds_overlap_ball(const kdtree_t * T,
                    const kdtree_node_t * node,
                    const double * Q)
{

    const double rmax = pqheap_get_max_value(T->pq);
    return bounds_overlap_ball_raw(node->bbx, Q, rmax);
}

/* Recursive search until no more points can be found
   Return 1 if we are done
   Return 0 else
*/

static int kdtree_search(kdtree_t * T, const kdtree_node_t * node, const double * Q)
{
    pqheap_t * pq = T->pq;

    if(node_is_final(node))
    {
        T->direct_path = 0;
        const double * NX = T->XID + node->data_idx;
        const size_t * NID = (size_t *) T->XID + node->data_idx;

        // Add all points
        for(size_t kk = 0; kk<node->n_points; kk++)
        {
            double d2 = eudist_sq(NX + kk*XID_STRIDE, Q);
            //X += 2;
            //printf("Adding %f\n", sqrt(d2));
            //fprintf(T->log, "%f ", d2);

            pqheap_insert(pq, d2, NID[kk*XID_STRIDE + KDTREE_DIM]);
        }

        double rmax = sqrt(pqheap_get_max_value(pq));

        int done =  within_bounds(node, Q, rmax);
        //printf("rmax = %f, done = %d\n", rmax, done);
        return done;
    }

    // Descend depending on pivot
    // First take the path that gets us closer to the query point
    int split_dim = node->split_dim;
    int done = 0;
    if(Q[split_dim] > node->pivot)
    {
        // correct direction
        if(T->direct_path || bounds_overlap_ball(T, T->nodes + node_right_child_id(node->id), Q))
        {
            done = kdtree_search(T, T->nodes + node_right_child_id(node->id), Q);
            if(done == 1)
            {
                return done;
            }
        }
        // "wrong direction"
        if(bounds_overlap_ball(T, T->nodes + node_left_child_id(node->id), Q))
        {
            done = kdtree_search(T, T->nodes + node_left_child_id(node->id), Q);
            if(done == 1)
            {
                return done;
            }
        }
    } else {
        // "correct" direction
        if(T->direct_path || bounds_overlap_ball(T, T->nodes + node_left_child_id(node->id), Q))
        {
            done = kdtree_search(T, T->nodes + node_left_child_id(node->id), Q);
            if(done)
            {
                return 1;
            }
        }

        // "wrong" direction
        if(bounds_overlap_ball(T, T->nodes + node_right_child_id(node->id), Q))
        {
            done = kdtree_search(T, T->nodes + node_right_child_id(node->id), Q);
        }
        if(done)
        {
            return 1;
        }
    }
    // Now we have added all sub regions so we can check if the ball falls
    // within this non-end-node as well

    double rmax = sqrt(pqheap_get_max_value(pq));

    if(within_bounds(node, Q, rmax))
    {
        return 1;
    }
    return 0;
}

size_t * kdtree_query_knn(kdtree_t * T, const double * Q, size_t k)
{
    if(k > T->n_points)
    {
        fprintf(stderr,
                "kdtree_query_knn error: Impossible to call for %zu points when\n"
                "there are only %zu in the tree\n", k, T->n_points);
        return NULL;
    }
    //    printf("-> Q = (%f, %f)\n", Q[0], Q[1]);

    // If k changed from the last query, update:
    if(T->result_alloc != k)
    {
        if(T->pq != NULL)
        {
            pqheap_free(&T->pq);
        }
        T->pq = NULL;
        free(T->result);
        T->result = NULL;
        T->result_alloc = 0;
    }


    // Set up priority queue
    if(T->pq == NULL)
    {
        T->pq = pqheap_new(k);
    }
    pqheap_t * pq = T->pq;
    pq->n = 0;
    pqheap_insert(pq, 1e99, 0);


    if(T->result == NULL)
    {
        T->result = calloc(k, sizeof(size_t));
        assert(T->result != NULL);
    }


    // Traverse the tree
    T->direct_path = 1;
    kdtree_search(T, T->nodes, Q);

    // If we don't need an ordered answer we could just traverse
    // the pq and extract the elements as we go.

    for(size_t kk = 0; kk<k; kk++)
    {
        double val = 0;
        uint64_t idx = 0;
        pqheap_pop(pq, &val, &idx);
        //printf("Popped: %lu, d = %f\n", idx, val);
        T->result[k-kk-1] = idx;
    }

    return T->result;
}


void kdtree_validate(kdtree_t * T)
{
#ifdef NDEBUG
    printf("kdtree_validate does not work when NDEBUG is defined\n");
    if(T == NULL)
    {
        printf("T is null\n");
    }
#else
    printf("kdtree_validate()\n");
    assert(T != NULL);
    assert(sizeof(double) == sizeof(size_t));
    assert(T->n_nodes_alloc > 0);

    // Check that all points are within bounds
    for(size_t n = 0 ; n < T->n_nodes_alloc; n++)
    {
        kdtree_node_t * node = T->nodes + n;
        if(node->n_points > 0)
        {


            //printf("Node id: %zu (Left %d, Right %d), %zu points\n", n,
            //       node->node_left, node->node_right, node->n_points*(node->node_left == -1));
            // TODO: check XID that the points are within bounds ...
            // Looks like things are wrong. Too many points in some leafs.

            //node_print_bbx(node);

            double * XID = T->XID + node->data_idx;
            if(node_is_final(node))
            {
                //print_XID(XID, node->n_points);
                for(size_t pp = 0 ; pp < node->n_points; pp++)
                {
                    for(size_t dd = 0; dd < KDTREE_DIM; dd ++)
                    {
                        assert(XID[pp*XID_STRIDE + dd] >= node->bbx[2*dd]);
                        assert(XID[pp*XID_STRIDE + dd] <= node->bbx[2*dd+1]);
                    }
                }
            }
        }
    }
    printf("done\n");
#endif
}

struct darray {
    size_t * data;
    size_t n_used;
    size_t n_alloc;
};

static void darray_n_more(struct darray * A, size_t nmore)
{
    if(A->n_used + nmore >= A->n_alloc)
    {
        size_t new_size = A->n_alloc + nmore;
        if(new_size < 1.2 *A->n_alloc)
        {
            new_size = 1.2*A->n_alloc;
        }
        A->data = realloc(A->data, new_size*sizeof(size_t));
        assert(A->data != NULL);
        A->n_alloc = new_size;
    }
}

static void _kdtree_query_radius(const kdtree_t * T,
                                 const double * Q,
                                 size_t node_id,
                                 const double r,
                                 const double r2,
                                 struct darray * res)
{
    kdtree_node_t * node = T->nodes + node_id;
    if( ! bounds_overlap_ball_raw(node->bbx, Q, r2) )
    {
        return;
    }

    /* If we reached a leaf see what points match the criteria */
    if(node_is_final(node))
    {
        double * X = T->XID + node->data_idx;
        size_t * ID = (size_t * ) T->XID + node->data_idx;
        darray_n_more(res, node->n_points);
        for(size_t kk = 0; kk < node->n_points; kk++)
        {
            if(eudist_sq(X + kk*XID_STRIDE, Q) < r2)
            {
                res->data[res->n_used] = ID[kk*XID_STRIDE + KDTREE_DIM];
                res->n_used++;
            }
        }
        return;
    }
    /* If not in a leaf, we see what children it makes sense to traverse */
    // Some linear algebra: If Q + (mid-Q)/||mid-Q||*r crosses any of the 6 faces
    // we need to check. Simpler way to determine ?
    // if (x_intersect) if (y_intersect) if (z_intersect) then traverse ...

    _kdtree_query_radius(T, Q,
                         node_left_child_id(node_id),
                         r, r2, res);

    _kdtree_query_radius(T, Q,
                         node_right_child_id(node_id),
                         r, r2, res);

    return;
}

size_t *
kdtree_query_radius(const kdtree_t * T,
                    const double * Q,
                    const double radius,
                    size_t * nfound)
{
    struct darray * res = calloc(1, sizeof(struct darray));
    res->n_alloc = 100;
    res->data = calloc(res->n_alloc, sizeof(size_t));

    _kdtree_query_radius(T, Q, 0, radius, pow(radius, 2), res);


    size_t * result = res->data;
    *nfound = res->n_used;
    free(res);
    return result;
}

static double gaussian(double d2, double sigma22)
{
    // sigma22 = 2*sigma^2
    // d2 = d^2
    return exp(-d2/sigma22);
}

static double _kdtree_kde(const kdtree_t * T,
                          const double * Q,
                          size_t node_id,
                          const double r2,
                          const double sigma22)
{
    kdtree_node_t * node = T->nodes + node_id;
    if( ! bounds_overlap_ball_raw(node->bbx, Q, r2) )
    {
        return 0;
    }

    double kde = 0;
    /* If we reached a leaf see what points match the criteria */
    if(node_is_final(node))
    {
        double * X = T->XID + node->data_idx;

        for(size_t kk = 0; kk < node->n_points; kk++)
        {
            double d2 = eudist_sq(X + kk*XID_STRIDE, Q);
            kde += gaussian(d2, sigma22);
        }
        return kde;
    }


    /* If not in a leaf, we see what children it makes sense to traverse */
    // Some linear algebra: If Q + (mid-Q)/||mid-Q||*r crosses any of the 6 faces
    // we need to check. Simpler way to determine ?
    // if (x_intersect) if (y_intersect) if (z_intersect) then traverse ...

    kde += _kdtree_kde(T, Q,
                       node_left_child_id(node_id),
                       r2, sigma22);

    kde += _kdtree_kde(T, Q,
                       node_right_child_id(node_id),
                       r2, sigma22);

    return kde;
}


double kdtree_kde(const kdtree_t * T,
                  const double * Q,
                  const double sigma,
                  const double cutoff)
{
    /* How distant points are of interest for the given sigma value ?
     */
    double r = 2.5*sigma;
    if(cutoff > 0.0)
    {
        r = cutoff*sigma;
    }
    return _kdtree_kde(T, Q, 0, pow(r,2.0), 2.0*pow(sigma, 2.0));
}

void kdtree_print_info(kdtree_t * T)
{
    printf("kdtree info:\n");
    printf("n_points: %zu\n", T->n_points);
    printf("Root node bbx: ");
    node_print_bbx(T->nodes);
    return;
}
