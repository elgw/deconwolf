#include "sortbox.h"

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef WINDOWS
#include <unistd.h>
#endif
#include <omp.h>

/* Optimization opportunity:
 *
 * Currently the sys time is almost 3x the user time, i.e.
 * the alloc/free operations consume most time.
 *
 * If we use pointer arrays for feature, rank and class,
 * Then we can split each of those vectors without a new allocation,
 * as we do in qsort etc.
 *
 * We would only need one sortbox object and could do something like
 * sortbox_split(B)
 *
 * save_pointers(B, level);
 * sortbox_set_ns(B, nleft);
 * process(B) // left child
 *
 * restore_pointers(B, level);
 * sortbox_set_ns(B, nright)
 * sortbox_advance(B, nleft)
 * process(B) // right child
 *
 * Or we always have the same pointers in the sortbox, and at each
 * recursion level we store on the stack: left first, n left, first
 * right position, n right.  since the max stack size is limited it
 * might be better to write it without recursion... but that would be
 * a later modification.
 */

/* Features:
 * f1,1 f2,1 f2,2 ... fnf,1
 * f1,2
 * f1,3
 * ...
 * f1,np              f(nf,np)*/

typedef uint8_t u8;
typedef uint32_t u32;
typedef uint64_t u64;
typedef double f64;

struct _sortbox{
    /* Number of features */
    u64 nf;
    /* Number of samples */
    u64 ns;
    /* The features in column-major format */
    double * feature;
    /* Corresponding ranks for all features */
    u32 * rank;
    /* The classes of the elements */
    u32 * class;
    u32 level;
    u32 max_class;

    /* Buffers used by the split function */
    double * f64_buffer; // for feature
    u32 * u32_buffer; // for rank and class
    u8 * u8_buffer; // for left/right indexing

    /* After subsampling this maps the features_ids to the original
     * features ids, i.e. feature_id_map[id] = original_id */
    u32 * feature_id_map;
};


/* Set n/N elements to 1 in use, the other to 0 */
static void random_selection(size_t N, size_t n,
                             u8 * use)
{
    assert(use != NULL);

    /* Case 1: Select minority */
    if(n < N/2)
    {
        memset(use, 0, N);
        size_t sel = 0;
        while(sel < n)
        {
            size_t pos = rand() % N;
            if(use[pos] == 0)
            {
                use[pos] = 1;
                sel++;
            }
        }
        return;
    }

    /* Case 2: Select majority */
    {
        memset(use, 1, N);

        size_t sel = N;
        while(sel != n)
        {
            size_t pos = rand() % N;
            if(use[pos] == 1)
            {
                use[pos] = 0;
                sel--;
            }
        }
        return;
    }
}


struct vali {
    double value;
    u32 index;
};

static void rank_vector_update(u32 * R, u32 nR, u32 max, u32 * buffer)
{
    memset(buffer, 0, (max+1)*sizeof(u32));
    u32 * A = buffer;

    // Mark elements that are used:
    for(size_t kk = 0; kk < nR; kk++)
    {
        //printf("R: %u\n", R[kk]);
        assert(R[kk] < max);
        A[R[kk]] = 1;
    }

    // Cumsum
    size_t sum = 0;
    for(size_t kk = 0; kk <= max; kk++)
    {
        sum += A[kk];
        A[kk] = sum - 1;
    }

    for(size_t kk = 0; kk < nR; kk++)
    {
        R[kk] = A[R[kk]];
        assert(R[kk] < max);
    }
    return;
}

/* Compare value and index struct */
static int
cmp_vali(const void * _a,
         const void * _b,
         __attribute__((unused)) void * p)
{
    struct vali * a = (struct vali *) _a;
    struct vali * b = (struct vali *) _b;
    if(a->value > b->value)
    {
        return 1;
    } else if(a->value < b->value)
    {
        return -1;
    }
    return 0;
}

/* Based on the values in F, set correponding ranks
 * in R */
static void
set_rank(u32 * restrict R,
         const double * restrict F,
         size_t n)
 {
     struct vali * VI = calloc(n, sizeof(struct vali));
     assert(VI != NULL);
    for(size_t kk = 0; kk < n; kk++)
    {
        VI[kk].value = F[kk];
        VI[kk].index = kk;
    }

    _quicksort(VI, n, sizeof(struct vali), cmp_vali, NULL);

    for(size_t kk = 0; kk < n; kk++)
    {
        R[ VI[kk].index ] = kk;
    }
    free(VI);
}

static void sortbox_init_buffers(sortbox * S)
{
    assert(S != NULL);
    const size_t nsample = S->ns;
    assert(nsample > 0);

    S->f64_buffer = calloc(nsample, sizeof(double));
    assert(S->f64_buffer != NULL);
    S->u32_buffer = calloc(nsample+1, sizeof(u32));
    assert(S->u32_buffer != NULL);
    S->u8_buffer = calloc(nsample, sizeof(u8));
    assert(S->u8_buffer != NULL);
}

sortbox * sortbox_init(const double * F,
                       const u32 * class,
                       size_t nsample,
                       size_t nfeature)
{
    u32 max_class = 0;
    for(size_t kk = 0; kk < nsample; kk++)
    {
        class[kk] > max_class ? max_class = class[kk] : 0;
    }
    if(max_class == 0)
    {
        /* Is this really a problem? Alternative: just continue... */
        printf("All samples have the same class\n");
        return NULL;
    }

    sortbox * S = calloc(1, sizeof(sortbox));
    assert(S != NULL);
    S->max_class = max_class;

    S->nf = nfeature;
    S->ns = nsample;

    S->feature = calloc(nsample*nfeature, sizeof(double));
    assert(S->feature != NULL);
    memcpy(S->feature, F, nsample*nfeature*sizeof(double));

    S->class = calloc(nsample, sizeof(u32));
    assert(S->class != NULL);
    memcpy(S->class, class, nsample*sizeof(u32));

    S->rank = calloc(nsample*nfeature, sizeof(u32));
    assert(S->rank != NULL);



// The memory grows alot here unless we limit the number
// of arenas used by malloc...
// The obvious work-around is to allocate all the buffers
// before going into the parallel region. Only having 1 arena
// most likely means that there will be a mutex for each malloc/free...


#pragma omp parallel for
    for(size_t ff = 0; ff < nfeature; ff++)
    {
        set_rank(S->rank + ff*nsample,
                 S->feature + ff*nsample,
                 nsample);
    }

    sortbox_init_buffers(S);

    return S;
}

sortbox *
sortbox_clone(sortbox * B)
{
    sortbox * C = calloc(1, sizeof(sortbox));
    assert(C != NULL);
    C->nf = B->nf;
    C->ns = B->ns;

    C->feature = calloc(C->nf*C->ns, sizeof(double));
    assert(C->feature != NULL);
    memcpy(C->feature, B->feature, C->nf*C->ns*sizeof(double));

    C->rank = calloc(C->nf*C->ns, sizeof(u32));
    assert(C->rank != NULL);
    memcpy(C->rank, B->rank, C->nf*C->ns*sizeof(u32));

    C->class = calloc(C->ns, sizeof(u32));
    assert(C->class != NULL);
    memcpy(C->class, B->class, C->nf*sizeof(u32));

    sortbox_init_buffers(C);

    return C;
}

static void split_double_vector(double * restrict V,
                                const u8 * restrict to_right,
                                const u32 nleft, const u32 nright,
                                double * restrict Rbuffer)
{
    memset(Rbuffer, 0, nright*sizeof(double));
    size_t wL = 0;
    size_t wR = 0;
    for(u32 kk = 0; kk < nleft+nright; kk++)
    {
        if(to_right[kk])
        {
            Rbuffer[wR++] = V[kk];
        } else {
            V[wL++] = V[kk];
        }
    }
    memcpy(V+nleft, Rbuffer, nright*sizeof(double));
}

static void split_u32_vector(u32 * restrict V,
                             const u8 * restrict to_right,
                             const u32 nleft, const u32 nright,
                             u32 * restrict Rbuffer)
{
    memset(Rbuffer, 0, nright*sizeof(u32));
    size_t wL = 0;
    size_t wR = 0;
    for(u32 kk = 0; kk < nleft+nright; kk++)
    {
        if(to_right[kk])
        {
            Rbuffer[wR++] = V[kk];
        } else {
            V[wL++] = V[kk];
        }
    }
    memcpy(V+nleft, Rbuffer, nright*sizeof(u32));
}


void sortbox_split(sortbox * P,
                   const int feature_id,
                   const double th,
                   // Index of first element
                   u32 start,
                   // Number of elements
                   u32 len,
                   // Return values:
                   u32 * _nleft, u32 * _nright)
{

    assert(feature_id >= 0);
    assert((u32) feature_id < P->nf);

    const size_t ns0 = P->ns;
    assert(ns0 > 0);


    // 1. Figure out what goes left and what goes right
    // and store in binary indicator

    u8 * to_right = P->u8_buffer;
    memset(to_right, 0, len);
    u32 nright = 0;

    {
        double * feature = P->feature + feature_id*P->ns + start;

        for(size_t kk = 0; kk < len; kk++)
        {
            if(feature[kk] >= th)
            {
                to_right[kk] = 1;
                nright++;
            }
        }
    }

    u32 nleft = len - nright;

    if(nleft == 0 || nright == 0)
    {
        *_nleft = 0;
        *_nright = 0;
        return;
    }

    *_nleft = nleft;
    *_nright = nright;

    // 2. Go through the data:
    //    - Split features.
    //    - Split order
    //    - Split class



    // Split class array
    split_u32_vector(P->class+start, to_right, nleft, nright, P->u32_buffer);


    for(size_t ff = 0; ff < P->nf; ff++)
    {
        // Feature ff
        split_double_vector(P->feature + ff*ns0 + start,
                            to_right,
                            nleft, nright,
                            P->f64_buffer);
        // Rank ff
        split_u32_vector(P->rank + ff*ns0 + start,
                         to_right,
                         nleft, nright,
                         P->u32_buffer);
        // Left rank
        rank_vector_update(P->rank + ff*ns0 + start,
                           nleft,
                           nleft+nright,
                           P->u32_buffer);
        // Right rank
        rank_vector_update(P->rank + ff*ns0 + start+ nleft,
                           nright,
                           nright+nleft,
                           P->u32_buffer);
    }

    return;
}

void sortbox_split_n(sortbox * P,
                     const int feature_id,
                     u32 nleft,
                     // Index of first element
                     u32 start,
                     // Number of elements
                     u32 len,
                     // Return values:
                     u32 * _nleft, u32 * _nright)
{

    assert(feature_id >= 0);
    assert((u32) feature_id < P->nf);

    const size_t ns0 = P->ns;
    assert(ns0 > 0);


    // 1. Figure out what goes left and what goes right
    // and store in binary indicator

    u8 * to_right = P->u8_buffer;
    memset(to_right, 0, len);
    u32 nright = 0;

    u32 * frank = P->rank + feature_id*P->ns + start;
    for(size_t kk = 0; kk < len; kk++)
    {
        if(frank[kk] >= nleft)
        {
            to_right[ kk ] = 1;
            nright++;
        }
    }

    if(0)
    {
    printf("split_n, feature_id=%u\n", feature_id);
    for(size_t kk = 0 ; kk < len; kk++)
    {
        printf("%u ", to_right[kk]);
    }

    printf("\n");
    printf("nleft = %u nright = %u len = %u\n", nleft, nright, len);
    }
    assert(nleft + nright == len);

    if(nleft == 0 || nright == 0)
    {
        *_nleft = 0;
        *_nright = 0;
        return;
    }

    *_nleft = nleft;
    *_nright = nright;

    // 2. Go through the data:
    //    - Split features.
    //    - Split order
    //    - Split class

    // Split class array
    split_u32_vector(P->class+start,
                     to_right,
                     nleft, nright,
                     P->u32_buffer);

    for(size_t ff = 0; ff < P->nf; ff++)
    {
        // Feature ff
        split_double_vector(P->feature + ff*ns0 + start,
                            to_right,
                            nleft, nright,
                            P->f64_buffer);
        // Rank ff
        split_u32_vector(P->rank + ff*ns0 + start,
                         to_right,
                         nleft, nright,
                         P->u32_buffer);
        // Left rank
        rank_vector_update(P->rank + ff*ns0 + start,
                           nleft,
                           nleft+nright,
                           P->u32_buffer);
        // Right rank
        rank_vector_update(P->rank + ff*ns0 + start+ nleft,
                           nright,
                           nright+nleft,
                           P->u32_buffer);
    }

    return;
}


void sortbox_free(sortbox * B)
{
    if(B == NULL)
    {
        return;
    }
    free(B->rank);
    free(B->class);
    free(B->feature);
    free(B->u8_buffer);
    free(B->f64_buffer);
    free(B->u32_buffer);
    free(B->feature_id_map);
    free(B);
}

/* The returned arrays are owned by the sortbox and will be overwritten
 * by other sortbox_commands */
void
sortbox_get_feature(const sortbox * B,
                    const u32 feature_id,
                    const u32 start,
                    const u32 len,
                    f64 ** _ret_feature,
                    u32 ** _ret_class)
{
    // printf("get_feature(%u)\n", feature_id);
    assert(B != NULL);
    assert(feature_id < B->nf);

    f64 * ret_feature = B->f64_buffer;
    u32 * ret_class = B->u32_buffer;
    *_ret_feature = ret_feature;
    *_ret_class = ret_class;

    assert(B->feature != NULL);
    double * feature = B->feature + (size_t) feature_id*B->ns + start;

    assert(B->rank != NULL);
    u32 * R = B->rank + (size_t) feature_id * B->ns + start;

    u32 * C = B->class + start;

    for(size_t kk = 0; kk < len; kk++)
    {
        //printf("kk = %zu f = %f r = %d\n", kk, feature[kk], R[kk]);
        assert(R[kk] < B->ns);
        ret_feature[ R[kk] ] = feature[kk];
        ret_class[ R[kk] ] = C[kk];
    }

    return;
}

size_t sortbox_get_nsample(const sortbox * B)
{
    return B->ns;
}

size_t sortbox_get_nfeature(const sortbox * B)
{
    return B->nf;
}

/* This is essential when constructing random forests, for each tree
   it is required that both the data points and the features are
   sampled. */

sortbox *
sortbox_subsample(const sortbox * B,
                  const size_t n_sample,
                  const size_t n_feature)
{
    assert(B != NULL);
    if(n_sample == 0)
    {
        return NULL;
    }
    if(n_feature == 0)
    {
        return NULL;
    }
    assert(n_sample <= B->ns);
    assert(n_feature <= B->nf);
    sortbox * S = calloc(1, sizeof(sortbox));
    assert(S != NULL);
    S->nf = n_feature;
    S->ns = B->ns;
    sortbox_init_buffers(S);
    S->ns = n_sample;
    S->max_class = B->max_class;

    // Initialize all buffers
    S->feature = calloc(n_sample*n_feature, sizeof(double));
    assert(S->feature != NULL);
    S->rank = calloc(n_sample*n_feature, sizeof(u32));
    assert(S->rank != NULL);
    S->class = calloc(n_sample, sizeof(u32));
    assert(S->class != NULL);

    // Select features and samples randomly
    u8 * sel_sample = S->u8_buffer;
    u8 * sel_feature = (u8*) S->u32_buffer;

    random_selection(B->ns, n_sample, sel_sample);
    random_selection(B->nf, n_feature, sel_feature);

    S->feature_id_map = calloc(n_feature, sizeof(u32));
    assert(S->feature_id_map != NULL);

    {
        size_t wpos = 0;
        for(size_t kk = 0; kk < B->nf; kk++)
        {
            if(sel_feature[kk])
            {
                assert(wpos < n_feature);
                S->feature_id_map[wpos++] = kk;
            }
        }
    }


    // Copy data
    size_t wpos = 0;
    for(size_t kk = 0; kk < B->ns; kk++)
    {
        if(sel_sample[kk])
        {
            S->class[wpos++] = B->class[kk];
        }
    }


    // Copy ranks and features

    u32 f_to = 0; // feature index to write to
    for(u32 f = 0; f < B->nf; f++)
    {
        if(sel_feature[f])
        {
            double * F_from = B->feature + f*B->ns;
            double * F_to = S->feature + f_to*S->ns;

            u32 * R_from = B->rank + f*B->ns;
            u32 * R_to = S->rank + f_to*S->ns;

            size_t writepos = 0;
            for(size_t kk = 0; kk < B->ns; kk++)
            {
                if(sel_sample[kk])
                {
                    F_to[writepos] = F_from[kk];
                    R_to[writepos] = R_from[kk];
                    writepos++;
                }
            }
            f_to++;
        }
    }

    // Update ranks
    for(u32 f = 0; f < S->nf; f++)
    {
        u32 * R = S->rank + f*S->ns;
        rank_vector_update(R, S->ns, B->ns, S->u32_buffer);
    }

    return S;
}

u32
sortbox_get_nclass(const sortbox * B)
{
    return B->max_class+1;
}

u32
sortbox_map_feature(const sortbox * B, u32 id)
{
    assert(B->feature_id_map != NULL);
    return B->feature_id_map[id];
}

const u32 *
sortbox_get_class_array_unsorted(const sortbox * B)
{
    return B->class;
}
