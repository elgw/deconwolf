#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifndef WINDOWS
#include <unistd.h>
#endif

#include <omp.h>

#include "entropy.h"
#include "gini.h"
#include "sortbox.h"
#include "trafo.h"
#include "trafo_util.h"

/* Compile with -fvisibility=hidden
 * then only functions decorated with this macro will be visible
 * in the .so file, check with
 * nm libtrafo.so | grep ' T '
 */
#if HAVE___ATTRIBUTE__VISIBILITY_HIDDEN
#define EXPORT __attribute__ ((visibility("default")))
#else
#define EXPORT
#endif

static void
dw_gettime(struct timespec * t)
{
#ifdef WINDOWS
    timespec_get(t, TIME_UTC); // since C11
#else
    clock_gettime(CLOCK_REALTIME, t);
#endif
    return;
}


enum trafo_criterion {
    trafo_gini, trafo_entropy
};

typedef struct  {
    u32 left_node; // below threshold, refering to the index of the node
    u32 right_node; // above threshold (or CLASS)
    i32 var; // Variable, set to -1 if not decided yet
    f32 th; // Threshold
} tnode;


/* For storing the table for a tree */
typedef struct {
    tnode * nodes;
    u32 nnode;
    size_t nalloc;
    size_t nf;
    size_t maxclass;
    double * importance;
} ttable;


/* For storing the settings, more or less a superset of the CLI
 * configuration  */
struct _trf {
    /*
     * Input data
     */
    u32 n_sample;
    u32 n_feature;

    /* Feature data, set one of these */
    const double * F_col_major;
    const double * F_row_major;

    /* Label data */
    const u32 * label;

    /*
     * Required settings
     */

    /* Number of trees to construct */
    u32 n_tree;

    u32 max_label; // The maximum label + 1;

    /* What to split by, 0 = gini, 1 = entropy */
    enum trafo_criterion criterion;

    /*
     * Optional settings
     */

    u8 verbose;

    /* Fraction of samples to use per tree. "bagging" if set to 0
     * the default will be used*/
    float tree_f_sample;
    /* Number of features to use per tree
     * if set to zero sqrt(n_features) will be used
     */
    u32 tree_n_feature;

    u32 min_samples_leaf;

    ttable * trees;
};




/* Transpose the M x N matrix into N x M. M is the size of the first,
 * non-strided dimension */
static double *
transpose_f64(const double * D, const size_t M, const size_t N)
{
    double * T = calloc(M*N, sizeof(f64));
    assert(T != NULL);
    for(size_t nn = 0; nn < N; nn++)
    {
        for(size_t mm = 0; mm < M; mm++)
        {
            T[nn + mm*N] = D[mm + nn*M];
        }
    }
    return T;
}

static u32
argmax_u32(const u32 * x, const size_t n)
{
    u32 am = 0;
    float vmax = x[0];
    for(size_t kk = 1; kk < n; kk++)
    {
        if(x[kk] > vmax)
        {
            am = kk;
            vmax = x[kk];
        }
    }
    return am;
}

static u32
most_commomax_label(const u32 * C, const size_t n, const size_t max_label_id)
{
    size_t nc = max_label_id+1;
    u32 H[nc];
    memset(H, 0, nc*sizeof(u32));
    for(size_t kk = 0; kk < n; kk++)
    {
        H[C[kk]]++;
    }
    return argmax_u32(H, nc);
}


static u32
ttable_classify_vector(const ttable * restrict ttab,
                       size_t ntree,
                       const double * restrict V,
                       u32 max_label_id)
{

    u32 n = max_label_id+1;
    u32 H[n];
    memset(H, 0, n*sizeof(u32));
    assert(H[0] == 0);

    for(size_t tt = 0; tt < ntree; tt++)
    {
        const ttable * T = &ttab[tt];
        tnode * node = &T->nodes[0];
        while(node->left_node != 0)
        {
            if(V[node->var] < node->th)
            {
                assert(node->left_node < ttab->nnode);
                node = &T->nodes[node->left_node];
            } else {
                assert(node->right_node < ttab->nnode);
                node = &T->nodes[node->right_node];
            }
        }
        H[node->right_node]++;
    }
    //printf("%f, %f, %f, %f, %u\n", V[0], V[1], V[2], V[3], argmax_u32(H, n));
    return argmax_u32(H, n);
}


static void
ttable_grow(ttable * T)
{
    size_t new_size = T->nalloc*1.2;
    tnode * old_location = T->nodes;
    #ifdef WINDOWS
    T->nodes = realloc(T->nodes, new_size*sizeof(tnode));
    #else
    T->nodes = reallocarray(T->nodes, new_size, sizeof(tnode));
    #endif
    assert(T->nodes != NULL);
    /* "If the new size is larger than the old size, the
       added memory will not be initialized".

       So we set it to 0 */
    if(T->nodes != old_location)
    {
        memset(T->nodes + T->nalloc, 0, (new_size-T->nalloc)*sizeof(tnode));
    }
    T->nalloc = new_size;
    return;
}

static void
recurse_tree(sortbox * B,
             ttable * T,
             const uint32_t minsize,
             const size_t level,
             const size_t start,
             const size_t npoints,
             const enum trafo_criterion criterion,
             double node_disorder,
             const size_t tnode_id)
{
    tnode * node =  &T->nodes[tnode_id];

    if(level == 0)
    {
        const u32 * cl_vector = sortbox_get_class_array_unsorted(B);
        if(criterion == trafo_gini)
        {
            node_disorder =
                gini_evaluate(cl_vector,
                              npoints,
                              sortbox_get_nclass(B));
            assert(node_disorder <= 1);
        } else {
            node_disorder =
                entropy_evaluate(cl_vector,
                                 npoints,
                                 sortbox_get_nclass(B));
        }
        assert(node_disorder >= 0);
    }

    if(node_disorder == 0 && level > 0)
    {
        assert(level > 0);
        // Finalize the node. All elements have the same class
        const u32 * cl = sortbox_get_class_array_unsorted(B);
        // tnode_id=95 node->left_node=3930079040
        // printf("tnode_id=%zu node->left_node=%u\n", tnode_id, node->left_node);
        fflush(stdout);

        node->left_node = 0;
        node->right_node = cl[start];

        for(size_t kk = start; kk < start+npoints; kk++)
        {
            assert(cl[kk] == node->right_node);
        }

        return;
    }

    //    printf("Level: %zu, start = %zu, npoints=%zu\n", level, start, npoints);
    if(npoints <= minsize)
    {
        //printf("Can't split any more\n");
        const u32 * cl_vector = sortbox_get_class_array_unsorted(B);
        node->left_node = 0;
        node->right_node =
            most_commomax_label(cl_vector+start, npoints,
                                sortbox_get_nclass(B));
        //printf("leaf: minsizes: class:%u\n", node->right_node);
        return;
    }

    int sel_feature = -1; // What feature to use
    double mindisorder = 1e99; // An upper bound on the disorder

    double th = 0;
    double disorder_left = 1;
    double disorder_right = 1;
    u32 best_nleft = 0;

    // TODO: In some cases we want to keep all the features but sample a fraction
    // of them at each split. Here we are cheating -- subsampling features at the start ...

    for(u32 kk = 0; kk < sortbox_get_nfeature(B); kk++)
    {
        //int ff = rand() % sortbox_get_nfeature(B);
        int ff = kk;
        //printf("ff = %u\n", ff);

        double * feature;
        uint32_t * class;

        sortbox_get_feature(B, ff, start, npoints,
                            // Output:
                            &feature, &class);
        if(0){
            printf("C=[");
            for(size_t tt = 0; tt < npoints; tt++)
            {
                printf("%u", class[tt]);
                if(tt + 1 != npoints)
                {
                    printf(", ");
                }
            }
            printf("]\n");
        }
        for(size_t kk = 0 ; kk+1 < npoints; kk++)
        {
            assert(feature[kk] <= feature[kk+1]);
        }


        u32 nleft = 0;
        u32 nright;

        // if criterion == gini
        f64 d_left, d_right;
        double disorder;

        if(criterion == trafo_gini)
        {
            disorder =
                gini_split(class, feature,
                           npoints, sortbox_get_nclass(B),
                           &nleft, &nright,
                           &d_left, &d_right);
            assert(disorder <= 1.0);
        } else {
            disorder =
                entropy_split(class, feature,
                              npoints, sortbox_get_nclass(B),
                              &nleft, &nright,
                              &d_left, &d_right);
            //printf("ff=%u E: %f -> %f (%f, %f) d=%f\n",ff,  node_disorder, disorder, d_left, d_right, disorder-node_disorder);
        }

        assert(disorder >= 0.0);

        if(nleft == 0)
        {
            /* This can happen if there is no threshold that can split
               the data. It is ok. */
        }

        if(disorder < mindisorder && nleft > 0)
        {
            mindisorder = disorder;
            best_nleft = nleft;
            sel_feature = ff;
            th = 0.5*(feature[nleft-1] + feature[nleft]);

            disorder_left = d_left;
            disorder_right = d_right;
            if(0)
            {
                printf("]\n");
                printf("sel_feature = %u\n", sel_feature);
                printf("nleft = %u, nright = %u\n", nleft, nright);
                printf("disorder = %f, dleft=%f, dright=%f\n", disorder,
                       disorder_left, disorder_right);
            }
        }
    }

    if(sel_feature == -1)
    {
        /* Unable to pick a feature to split on. This can happen if
           features are constant and typically happens more often when
           constructing trees from sub sampled data.

           In this case the tree will not be perfect, i.e. it will not
           be able to classify all training data correctly.
        */

        const u32 * cl = sortbox_get_class_array_unsorted(B);
        node->left_node = 0;
        node->right_node =
            most_commomax_label(cl+start, npoints, sortbox_get_nclass(B));
        return;
    }


    assert(sortbox_map_feature(B, sel_feature) >= (u32) sel_feature);
    node->th = th;
    /* Note: the variable is mapped back to the original variable id
     * in case that we are processing a subselection of variabes */
    node->var = sortbox_map_feature(B, sel_feature);

    //printf("\n --- will split over feature %u\n", sel_feature);

    uint32_t nleft, nright;

    /* Split the sortbox with n elements to the left */
    sortbox_split_n(B,
                    sel_feature, best_nleft,
                    start, npoints,
                    &nleft, &nright);


    assert(nleft == best_nleft);

    if(nleft == 0 || nright == 0)
    {
        assert(0); // this should not happen
        // printf("Can't split any more\n");
        return;
    }
    assert(nleft + nright == npoints);

    /* Grow size if needed */
    if(T->nnode + 3 >= T->nalloc)
    {
        ttable_grow(T);
        /* the address of T->nodes might be changed */
        node = &T->nodes[tnode_id];
    }

    /* Reserve slots for the children nodes in the table */
    size_t left_id = ++T->nnode;
    size_t right_id = ++T->nnode;
    node->left_node = left_id;
    node->right_node = right_id;

    /* Record variable and decrease of disorder as a proxy for importance */

    if(0){
    printf("l=%zu v=%u %f (%f, %f), L=%u, R=%zu\n", level, node->var,
           node_disorder, disorder_left, disorder_right, best_nleft, (npoints - best_nleft));
    }


    if(criterion == trafo_gini)
    {
        T->importance[node->var]
            += ( (double) npoints*node_disorder -
                 ( (double) best_nleft*disorder_left + (double) (npoints-best_nleft)*disorder_right ) );
    } else if (criterion == trafo_entropy){
        T->importance[node->var]
            += (node_disorder - mindisorder);
    }




    /* Or possibly just the number of affected samples  */
    //T->importance[sortbox_map_feature(B,sel_feature)] += npoints;

    //printf("Running Left\n");
    recurse_tree(B, T,
                 minsize, level+1, start, nleft,
                 criterion,
                 disorder_left, left_id);
    //printf("Running Right\n");
    recurse_tree(B, T,
                 minsize, level+1, start+nleft, nright,
                 criterion,
                 disorder_right, right_id);
    return;
}

EXPORT void
trafo_free(trf * s)
{
    if(s == NULL)
    {
        return;
    }
    if(s->trees != NULL)
    {
        for(size_t kk = 0; kk < s->n_tree; kk++)
        {
            free(s->trees[kk].nodes);
            free(s->trees[kk].importance);
        }
        free(s->trees);
    }

    free(s);
    return;
}

EXPORT void
trafo_print(FILE * fid, const trf * s)
{
    int got_features = 0;
    if(s->F_col_major != 0)
    {
        fprintf(fid, "Features provided in column major format\n");
        got_features++;
    }
    if(s->F_row_major != 0)
    {
        fprintf(fid, "Features provided in row major format (to be transposed)\n");
        got_features++;
    }
    if(got_features == 0)
    {
        fprintf(fid, "Warning: No features provided\n");
    }
    if(got_features > 1)
    {
        fprintf(fid, "ERROR: Ambiguous feaures provided\n");
    }
    if(s->label != NULL)
    {
        fprintf(fid, "Label array provided\n");
    } else {
        fprintf(fid, "Warning: Class array missing\n");
    }

    fprintf(fid, "Number of features: %u\n", s->n_feature);
    fprintf(fid, "Number of samples: %u\n", s->n_sample);
    fprintf(fid, "Number of trees: %u\n", s->n_tree);
    fprintf(fid, "Fraction of samples per tree: %.2f\n", s->tree_f_sample);
    fprintf(fid, "Features per tree: %u\n", s->tree_n_feature);
    fprintf(fid, "min_samples_leaf: %u\n", s->min_samples_leaf);
    fprintf(fid, "Largest label id: %u\n", s->max_label);
    fprintf(fid, "Splitting criterion: ");
    switch (s->criterion)
    {
    case trafo_gini:
        fprintf(fid, "Gini Impurity\n");
        break;
    case trafo_entropy:
        fprintf(fid, "Entropy\n");
        break;
    default:
        fprintf(fid, "INVALID!\n");
        break;
    }
    return;
}


static int
trafo_check(trf * s)
{
    if(s->tree_f_sample <= 0)
    {
        s->tree_f_sample = 0.632;
        if(s->tree_f_sample > 1)
        {
            s->tree_f_sample = 1;
        }
    }
    if(s->tree_n_feature == 0)
    {
        s->tree_n_feature = ceil(sqrt(s->n_feature));
    }
    if(s->min_samples_leaf == 0)
    {
        s->min_samples_leaf = 1;
    }
    if(s->n_tree == 0)
    {
        printf("Error: Invalid number of trees\n");
        return -1;
    }
    if(s->n_feature == 0)
    {
        printf("Error: Number of features has to be at least one\n");
        return -1;
    }
    if(s->n_sample == 0)
    {
        printf("Error: Number of samples has to be at least one\n");
        return -1;
    }
    if(s->tree_n_feature > s->n_feature)
    {
        printf("Error: Invalid number of features per tree\n");
    }
    if(s->tree_f_sample > 1 || s->tree_f_sample < 0)
    {
        printf("Error: Invalid samples per tree fraction\n");
    }

    return 0;
}


/* Predict the class for n_point points in X. L is optional. If
   provided the function will print out the number of correctly
   predicted points. */

EXPORT u32 *
trafo_predict(trf * s,
              const f64 * X_cm,
              const f64 * X_rm,
              u64 n_point)
{
    if( (X_rm == NULL) & (X_cm == NULL) )
    {
        return NULL;
    }
    if( (X_rm != NULL) & (X_cm != NULL) )
    {
        return NULL;
    }
    struct timespec t0, t1;
    dw_gettime(&t0);

    u32 * P = calloc(n_point, sizeof(u32));
    assert(P != NULL);

    if(s->verbose > 0)
    {
        printf("Classifying using %u tables/trees\n", s->n_tree);
    }


    if(X_cm != NULL)
    {
        //printf("CM predictions. s->n_feature = %u\n", s->n_feature);
        const double * X = X_cm;

#pragma omp parallel for
        for(size_t ss = 0; ss < n_point; ss ++)
        {
            double V[s->n_feature];
            for(u32 ff = 0; ff < s->n_feature; ff++)
            {
                V[ff] = X[ff*n_point + ss];
            }
            u32 H[s->max_label+1];
            memset(H, 0, (s->max_label+1)*sizeof(u32));
            for(u32 tt = 0; tt < s->n_tree; tt++)
            {
                u32 res = ttable_classify_vector(s->trees + tt, 1,
                                                 V, s->max_label+1);
                H[res] ++ ;
            }
            P[ss] = argmax_u32(H, s->max_label+1);
        }
    } else {
        /* If ntree == 1 we don't need H etc, that could be a separate case  */
#pragma omp parallel for
        for(size_t ss = 0; ss < n_point; ss ++)
        {
            u32 H[s->max_label+1];
            memset(H, 0, (s->max_label+1)*sizeof(u32));
            i64 class = -1;
            for(u32 tt = 0; tt < s->n_tree; tt++)
            {
                u32 res = ttable_classify_vector(s->trees + tt, 1,
                                                 X_rm + ss*s->n_feature,
                                                 s->max_label+1);
                H[res] ++ ;

                /* If more than 50% of the votes are in a single bin,
                   the final result can be predicted and there is no
                   need to continue */
                if(2*H[res] > s->n_tree){
                    class = res;
                    break;
                }
            }
            if(class == -1)
            {
                class = argmax_u32(H, s->max_label+1);
            }
            assert(class >= 0);
            P[ss] = class;
        }
    }
    dw_gettime(&t1);
    if(s->verbose > 0)
    {
        printf("Prediction took %f s\n", timespec_diff(&t1, &t0));
    }
    return P;
}


EXPORT trf *
trafo_fit(trafo_settings * conf)
{
    trf * s = calloc(1, sizeof(trf));
    s->F_col_major = conf->F_col_major;
    s->F_row_major = conf->F_row_major;
    s->label = conf->label;
    s->n_sample = conf->n_sample;
    s->n_tree = conf->n_tree;
    s->n_feature = conf->n_feature;
    s->tree_n_feature = conf->tree_n_feature;
    s->tree_f_sample = conf->tree_f_sample;
    s->verbose = conf->verbose;
    s->min_samples_leaf = conf->min_samples_leaf;
    if(conf->entropy)
    {
        s->criterion = trafo_entropy;
    } else {
        s->criterion = trafo_gini;
    }

    if(trafo_check(s))
    {
        printf("Invalid settings, unable to continue\n");
        return NULL;
    }

    const double * X = s->F_col_major;
    double * XT = NULL;
    int free_X = 0;
    if(X == NULL)
    {
        free_X = 1;
        XT = transpose_f64(s->F_row_major,  s->n_feature, s->n_sample);
        X = XT;
    }

    assert(X != NULL);

    {
        size_t mc = 0;
        for(size_t kk = 0; kk < s->n_sample; kk++)
        {
            s->label[kk] > mc ? mc = s->label[kk] : 0;
        }
        s->max_label = mc;
    }

    if(s->verbose > 0)
    {
        trafo_print(stdout, s);
    }


    struct timespec tictoc_start, tictoc_end;

    dw_gettime(&tictoc_start);

    sortbox * B = sortbox_init(X, s->label, s->n_sample, s->n_feature);

    if(B == NULL)
    {
        printf("Unable to initialize the data structures\n");
        trafo_free(s);
        return NULL;
    }
    if(free_X)
    {
        free(XT);
        X = NULL;
    }

    /* Initialize the tables */
    u32 tree_n_sample = s->tree_f_sample*s->n_sample;
    tree_n_sample < 1 ? tree_n_sample = 1 : 0;
    tree_n_sample > s->n_sample ? tree_n_sample = s->n_sample : 0;
    ttable * ttab = calloc(s->n_tree, sizeof(ttable));
    assert(ttab != NULL);
    for(size_t kk = 0 ; kk < s->n_tree; kk++)
    {
        ttab[kk].nalloc = tree_n_sample; // Will grow as needed
        ttab[kk].nodes = calloc(ttab[kk].nalloc, sizeof(tnode));
        assert(ttab[kk].nodes != NULL);
        ttab[kk].importance = calloc(s->n_feature, sizeof(double));
        assert(ttab[kk].importance != NULL);
    }


#pragma omp parallel for
    for(size_t tt = 0; tt < s->n_tree; tt++)
    {
        // printf("Tree = %zu\n", tt); fflush(stdout);
        //sortbox * BC = sortbox_clone(B);

        sortbox * BC = sortbox_subsample(B,
                                         tree_n_sample,
                                         s->tree_n_feature);
        if(BC == 0)
        {
            continue;
        }
        assert(sortbox_get_nsample(BC) == tree_n_sample);
        assert(sortbox_get_nfeature(BC) == s->tree_n_feature);
        recurse_tree(BC,
                     &ttab[tt],
                     s->min_samples_leaf, // leaf termination condition
                     0, // level
                     0, // start index
                     tree_n_sample, // number of selected samples
                     s->criterion,
                     1.0, // node_disorder (does not matter for the root)
                     0); // node_id

        ttab[tt].nnode++;
        sortbox_free(BC);
    }


    sortbox_free(B);
    dw_gettime(&tictoc_end);
    if(conf->verbose > 1)
    {
        printf(" Forest training took %f s\n",
               timespec_diff(&tictoc_end, &tictoc_start));
    }
    s->trees = ttab;
    return s;
}

static int
test_transpose_f64(void)
{

    size_t M = 10 + rand() % 10;
    size_t N = 10 + rand() % 10;
    double * X = calloc(M*N, sizeof(f64));
    assert(X != NULL);
    for(size_t kk = 0; kk < M*N; kk++)
    {
        X[kk] = (double) rand() / (double) RAND_MAX;
    }

    double * T = transpose_f64(X, M, N);
    double * TT = transpose_f64(T, N, M);

    for(size_t kk = 0; kk < M*N; kk++)
    {
        assert(X[kk] == TT[kk]);
    }
    free(X);
    free(T);
    free(TT);
    return 0;
}

EXPORT int
trafo_ut(void)
{
    test_transpose_f64();
    return 0;
}


/* Input / Output */

const u32 trafo_forest_magic = 0x018A08FA;
const u32 trafo_tree_magic = 0x3ADEFA12;

static int ttable_to_file(ttable * T, FILE * fid)
{
    size_t nwritten;
    nwritten = fwrite(&trafo_tree_magic, sizeof(u32), 1, fid);
    if(nwritten != 1)
    {
        printf("Error writing to disk (magic)\n");
        return 1;
    }
    nwritten = fwrite(&(T->nnode), sizeof(u32), 1, fid);
    if(nwritten != 1)
    {
        printf("Error writing to disk (n_nodes)\n");
        return 1;
    }
    #ifndef NDEBUG
    printf("Writing %u nodes to disk\n", T->nnode);
    #endif
    nwritten = fwrite(T->nodes, sizeof(tnode),
                      T->nnode, fid);
    if(nwritten != T->nnode)
    {
        printf("Error writing to disk (nodes)\n");
        return 1;
    }
    return 0;
}


static int ttable_from_file(ttable * T, FILE * fid)
{
    u32 magic;
    size_t nread = fread(&magic, sizeof(u32), 1, fid);
    if(nread != 1)
    {
        printf("Error reading from file (magic)\n");
        return -1;
    }
    if(magic != trafo_tree_magic)
    {
        printf("Error reading from file: wrong magic number\n");
    }

    u32 n_nodes;
    nread = fread(&n_nodes, sizeof(u32), 1, fid);
    if(nread != 1)
    {
        printf("Error reading from file (n_nodes)\n");
        return -1;
    }

#ifndef NDEBUG
    printf("Reading %u nodes from disk\n", n_nodes);
#endif

    T->nodes = calloc(n_nodes, sizeof(tnode));
    assert(T->nodes != NULL);
    T->nnode = n_nodes;

    nread = fread(T->nodes, sizeof(tnode), n_nodes, fid);
    if(nread != n_nodes)
    {
        fprintf(stderr, "Error reading nodes from disk (got %zu, expected %u)\n",
                nread, n_nodes);
        free(T->nodes);
        free(T);
        return -1;
    }

    return 0;
}

EXPORT int
trafo_save(trf * F,
               const char * filename)
{
    assert(F != NULL);
    FILE * fid = fopen(filename, "wb");
    if(fid == NULL)
    { goto fail2; }

    size_t nwritten;
    nwritten = fwrite(&trafo_forest_magic, sizeof(u32), 1, fid);
    if(nwritten != 1)
    { goto fail1; }

    nwritten = fwrite(&(F->n_tree), sizeof(u32), 1, fid);
    if(nwritten != 1)
    { goto fail1; }

    nwritten = fwrite(&F->n_feature, sizeof(u32), 1, fid);
    if(nwritten != 1)
    { goto fail1; }

    nwritten = fwrite(&F->max_label, sizeof(u32), 1, fid);
    if(nwritten != 1)
    { goto fail1; }

    for(size_t kk = 0; kk < F->n_tree; kk++)
    {
        if(ttable_to_file(&F->trees[kk], fid))
        { goto fail1; }
    }

    fclose(fid);

    return 0;

 fail1:
    printf("Error writing to %s\n", filename);
    fclose(fid);

 fail2:
    printf("Error opening %s for writing\n", filename);
    return -1;
}

EXPORT trf *
trafo_load(const char * filename)
{
    FILE * fid = fopen(filename, "rb");
    if(fid == NULL)
    {
        printf("Unable to open %s\n", filename);
        return NULL;
    }
    u32 magic;
    size_t nread = fread(&magic, sizeof(u32), 1, fid);
    if(nread != 1)
    {
        printf("Error reading from file (magic)\n");
        return NULL;
    }
    if(magic != trafo_forest_magic)
    {
        printf("Error reading from file: wrong magic number\n");
        fclose(fid);
        return NULL;
    }

    u32 n_tree;
    nread = fread(&n_tree, sizeof(u32), 1, fid);
    if(nread != 1)
    {
        printf("Error reading from file (n_ntree)\n");
        return NULL;
    }

    u32 n_feature;
    nread = fread(&n_feature, sizeof(u32), 1, fid);
    if(nread != 1)
    {
        printf("Error reading from file (n_feature)\n");
        return NULL;
    }

    u32 max_label;
    nread = fread(&max_label, sizeof(u32), 1, fid);
    if(nread != 1)
    {
        printf("Error reading from file (max_label)\n");
        return NULL;
    }

    trf * F = calloc(1, sizeof(trf));
    assert(F != NULL);
    F->n_tree = n_tree;
    F->n_feature = n_feature;
    F->trees = calloc(F->n_tree, sizeof(ttable));
    F->max_label = max_label;
    assert(F->trees != NULL);

    for(size_t kk = 0; kk < n_tree; kk++)
    {
        if(ttable_from_file(&F->trees[kk], fid))
        {
            printf("Error reading tree %zu\n", kk);
            assert(0);
        }
    }

    fclose(fid);
    return F;
}


EXPORT double *
trafo_importance(trf * F)
{
    if(F == NULL)
    {
        return NULL;
    }

    if(F->n_tree == 0)
    {
        return NULL;
    }

    if(F->trees == NULL)
    {
        return NULL;
    }

    if(F->trees[0].importance == NULL)
    {
        return NULL;
    }

    double * imp = calloc(F->n_feature, sizeof(double));
    assert(imp != NULL);
    memset(imp, 0, F->n_feature*sizeof(double));

    for(u32 kk = 0; kk < F->n_tree; kk++)
    {
        for(u32 ff = 0; ff < F->n_feature; ff++)
        {
            imp[ff] += F->trees[kk].importance[ff];
        }
    }
    double sum = 0;
    for(u32 ff = 0; ff < F->n_feature; ff++)
    {
        sum+=imp[ff];
    }
    for(u32 ff = 0; ff < F->n_feature; ff++)
    {
        imp[ff] /= sum;
    }

    return imp;
}
