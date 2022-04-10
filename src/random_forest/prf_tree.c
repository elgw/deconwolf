#include "prf_tree.h"

//
// Internal functions / Forward declarations
//

PrfNode * prf_node_new();


// Gini from a vector of classes
float gini_from_C(const float * C, size_t N, int nclass);

float gini_from_VC(float * VC, size_t N, float * g, size_t * nL, size_t * nR, float * lgini, float * rgini);
float gini_from_H(size_t nL, size_t * HL,
                   size_t nR, size_t * HR,
                   size_t nclass,
                   float * lgini, float * rgini);
int most_common_class(const float * C, const size_t N, const int nclass);
float squaresum(size_t * V, size_t T, size_t N);


PrfTree * prf_tree_new()
{
    PrfTree * T = malloc(sizeof(PrfTree));
    prf_tree_set_defaults(T);
    return T;
}

void prf_tree_set_defaults(PrfTree *T)
{
    T->maxClass = PRF_AUTOMATIC;
    T->min_size = 10;
    T->verbose = 0;
    T->max_features = PRF_AUTOMATIC;
    // A list where the variables to use are
    // randomized
    T->nodes_root = NULL;
    T->node_pointer = NULL; // Where to insert the next node
    T->varlist = NULL;
    T->best_idx = NULL;
    T->VC = NULL;
    T->feature_map = NULL;
    T->n_nodes = 0;
}

void prf_tree_show(FILE * f, PrfTree * T)
{
    fprintf(f, "Largest integer number used for a class\n");
    fprintf(f, "maxClass: ");
    if(T->maxClass == PRF_AUTOMATIC)
    {
        printf(" automatic \n");
    } else {
        printf("%d \n", T->maxClass);
    }
    fprintf(f, "Don't split nodes with fewer than min_size samples\n");
    fprintf(f, "min_size: %d\n", T->min_size);
    fprintf(f, "verbose: %d\n", T->verbose);
    fprintf(f, "The number of features to consider when splitting a node\n");
    fprintf(f, "max_features: ");
    if(T->max_features == PRF_AUTOMATIC)
    {
        printf(" automatic \n");
    } else {
        printf("%d\n", T->max_features);
    }
    return;
}

void prf_tree_free(PrfTree * T)
{
    if(T == NULL)
    {
        return;
    }
    if(T->best_idx != NULL)
    {
        free(T->best_idx);
    }
    if(T->VC != NULL)
    {
        free(T->VC);
    }
    free(T->varlist);
    if(T->nodes_root != NULL)
    {
        free(T->nodes_root);
    }
    if(T->randperm_buff != NULL)
    {
        free(T->randperm_buff);
    }
    if(T->feature_map != NULL)
    {
        free(T->feature_map);
    }
    free(T);

}

size_t prf_node_enumerate_from(PrfNode * T, size_t N)
{
    if(T == NULL)
    {
        return N;
    }
    T->id = N++;
    N = prf_node_enumerate_from(T->left, N);
    N = prf_node_enumerate_from(T->right, N);
    return N;
}

size_t prf_tree_enumerate(PrfTree * T)
{
    return prf_node_enumerate_from(T->nodes_root, 0);
}

float prf_tree_classify(const PrfTree * TT, const float * X)
{
    PrfNode * T = TT->nodes_root;

    while(T->class == PRF_NODE_NO_CLASS)
    {
        // Index of the variable to consider
        int var = T->var;
        // Threshold
        float th = T->th;
        if( X[var] < th )
        {
            assert(T->left != NULL);
            T = T->left;
        } else {
            assert(T->right != NULL);
            T = T->right;
        }
    }
    return T->class;
}

int float_cmp(const void * a, const void * b, __attribute__((unused)) void * p)
{

    if (*(float*)a > *(float*)b)
    {
        return 1;
    }
    else if (*(float*)a < *(float*)b)
    {
        return -1;
    }
    else
    {
        return 0;
    }

}

void prf_node_set_defaults(PrfNode * N)
{
    N->left = 0;
    N->right = 0;
    N->th = NAN;
    N->var = -1;
    N->class = PRF_NODE_NO_CLASS;
    N->parent = NULL;
    N->left = NULL;
    N->right = NULL;
    return;
}

PrfNode * prf_node_new()
{
    PrfNode * N = malloc(sizeof(PrfNode));
    prf_node_set_defaults(N);
    return N;
}

void prf_node_free(PrfNode * N)
{
    if(N->left != NULL)
    {
        prf_node_free(N->left);
    }
    if(N->right != NULL)
    {
        prf_node_free(N->right);
    }
    free(N);
    return;
}

int most_common_class(const float * C, const size_t N, const int nclass)
{
    // Figure out the most common class of C
    // C is supposed to contain integers
    // in the range 0, ..., nclass
    size_t * H = calloc(nclass+1, sizeof(size_t));
    for(size_t kk = 0; kk < N; kk++)
    {
        assert( (int) C[kk] > -1);
        assert( (int) C[kk] <= nclass);

        H[ (size_t) C[kk] ]++;
    }
    size_t max = 0;
    size_t maxclass = 0;
    for(int kk = 0; kk<= nclass; kk++)
    {
        if(H[kk] > max)
        {
            max = H[kk];
            maxclass = kk;
        }
    }
    free(H);
    return maxclass;
}


float vector_max(const float * V, const size_t N)
{
    float max = V[0];
    for(size_t kk = 0; kk<N; kk++)
    {
        if(V[kk] > max)
        {
            max = V[kk];
        }
    }
    return max;
}


void prf_tree_init_varlist(PrfTree * T)
{
    //printf("Init varlist\n");
    assert(T->n_features > 0);
    if(T->varlist == NULL)
    {
        //printf("Allocating for %zu features\n", c->n_features);
        assert(T->maxClass > 0);
        assert(T->varlist == NULL);
        T->varlist = malloc((T->n_features)*sizeof(int));
        for(size_t kk = 0; kk < T->n_features; kk++)
        {
            T->varlist[kk] = kk;
        }
    }
}

int cmp_int(const void * a, const void * b, __attribute__((unused)) void * p)
{
    int * A = (int *) a;
    int * B = (int *) b;

    if(A[0] > B[0])
    {
        return 1;
    } else {
        if(A[0] < B[0])
        {
            return -1;
        } else {
            return 0;
        }
    }
}

void randperm_int(int * X, size_t N, int * X2)
{
    assert(malloc_usable_size(X2) >= 2*N);
    //int * X2 = malloc(2*N*sizeof(int));
    for(size_t kk = 0; kk<N; kk++)
    {
        X2[2*kk] = rand() % INT_MAX; // Not truly random even if rand() would be
        X2[2*kk+1] = X[kk];
    }
    assert(N < 100);
    //printf("%zu>", N); fflush(stdout);
    _quicksort(X2, N, 2*sizeof(int), cmp_int, NULL);
    //printf("<"); fflush(stdout);
    for(size_t kk = 0; kk<N; kk++)
    {
        X[kk] = X2[2*kk+1];
    }
    //free(X2);
}

void prf_tree_varlist_randperm(PrfTree *T)
{
    randperm_int(T->varlist, T->n_features, T->randperm_buff);
}

void prf_tree_train(PrfTree * T, const float * X)
{
    T->nodes_root = malloc(2*T->n_samples*sizeof(PrfNode));
    T->n_nodes_alloc = 2*T->n_samples;
    T->node_pointer = T->nodes_root;
    T->randperm_buff = malloc(2*T->n_features*sizeof(int));
    prf_tree_train_with_parent(T,
                               NULL, // parent node
                               T->n_samples, X, -1);
}


void split_X(const float * V, const float * X,
             float ** P_Ltab, float ** P_Rtab,
             size_t N, size_t M,
             size_t nL, size_t nR, float th)
{
// TODO: This might break if there are duplicate x-values
// I guess that the number of points below the threshold should be counted...
float * Ltab = malloc(nL * M * sizeof(float));
float * Rtab = malloc(nR * M * sizeof(float));
size_t ll = 0; // Row to write to in Ltab
size_t rr = 0; // Row to write to in Rtab

for(size_t rpos = 0; rpos < N; rpos++)
 {
     // printf("-> %zu: x=%f\n", rpos, V[rpos]);
     if(V[rpos] < th)
     {
         assert(ll < nL);
         for(size_t mm = 0; mm < M; mm++)
         {
             Ltab[ll+mm*nL] = X[rpos + mm*N];
         }
         ll++;
     } else {
         assert(rr < nR);
         for(size_t mm = 0; mm < M; mm++)
         {
             Rtab[rr+mm*nR] = X[rpos + mm*N];
         }
         rr++;
     }
 }
 P_Ltab[0] = Ltab;
 P_Rtab[0] = Rtab;

assert(rr == nR);
assert(ll == nL);
}



PrfNode * prf_tree_train_with_parent(PrfTree * T, PrfNode * parent,
                                     const size_t N,
                                     const float * X,
                                     float gini)
{


    const size_t M = T->n_features+1;
    if(T->verbose > 3)
    {
        printf("N= %zu\n", N);
    }

    const float * CL = X+(M-1)*N; // class vector

    if(T->maxClass == PRF_AUTOMATIC)
    {
        T->maxClass = (int) vector_max(CL, N) + 1;
    }

    T->n_nodes++;
    assert(T->n_nodes < T->n_nodes_alloc);
    PrfNode * root = T->node_pointer++;
    prf_node_set_defaults(root);
    root->N = N;

    if(parent == NULL)
    {
        // Some initialization if this is the root node of the tree.
        root->gini = gini_from_C(CL, N, T->maxClass);
        if(T->max_features == -1)
        {
            T->max_features = floor(sqrt(T->n_features));
            if(T->max_features < 1)
            {
                T->max_features = 1; // At least one features
            }
        }
        prf_tree_init_varlist(T);
        T->best_idx = malloc(N*sizeof(uint32_t));
        T->VC = malloc(3*N*sizeof(float));
    } else {
        //printf("Parent told: %f, I think %f\n", gini, gini_from_C(C, N, T->maxClass));
        root->gini = gini;
    }

    // Determine if this node is final, i.e.,
    // has a decision and no children.

    int final = 0;
    if((int) N <= T->min_size)
    {
        final = 1;
    }
    if(root->gini == 0)
    {
        final = 1;
    }
    if(parent != NULL)
    {
        if(parent->gini == 0)
        {
            final = 1;
        }

        if((int) parent->N <= 2*T->min_size)
        {
            final = 1;
        }
    }
    if(final)
    {
        root->class = most_common_class(CL, N, T->maxClass);
        assert(root->class <= T->maxClass);
        assert(root->left == NULL);
        assert(root->right == NULL);
//        printf("Final node with class: %d\n", root->class);
        return root;
    }


    /* Make an interlaced copy of the variable and class vectors
     * This is done in order to be able to sort the data efficiently with qsort
     */

    float min_gini = INFINITY;
    int best_var = 0;
    float best_th = -1;
    float best_lgini = -1;
    float best_rgini = -1;
    uint32_t * best_idx = T->best_idx;

    // Randomize the list of variables to consider
    prf_tree_varlist_randperm(T);

    // value, class, idx
    float * VC = T->VC;

    assert(T->max_features <= (int) T->n_features);

    // Check which features is best to split over
    for(int vv = 0; vv < T->max_features; vv++)
    {
        int var = T->varlist[vv];
        const float * V = X+N*var; // variable vector


    for(size_t kk = 0; kk<N; kk++)
    {
        VC[3*kk] = V[kk];
        VC[3*kk+1] = CL[kk];
        VC[3*kk+2] = kk;
    }

    //qsort(VC, N, 3*sizeof(float), float_cmp);
    //printf("("); fflush(stdout);
    _quicksort(VC, N, 3*sizeof(float), float_cmp, NULL);
    //printf(")"); fflush(stdout);

    for(size_t kk = 0; kk<(N-1); kk++)
    {
        assert( VC[3*kk] <= VC[3*(kk+1)] );
    }

    // TODO: For each variable
    /* Get the threshold that minimizes Gini's impurity */
    float g = 0;
    size_t nL = 0, nR = 0;
    float lgini = 0;
    float rgini = 0;
    float th = gini_from_VC(VC, N, &g, &nL, &nR, &lgini, &rgini);
    if(g < min_gini)
    {
        min_gini = g;
        best_var = var;
        best_th = th;
        best_lgini = lgini;
        best_rgini = rgini;
        for(size_t ii = 0; ii < N; ii++)
        {
            best_idx[ii] = VC[3*ii+2];
        }
    }
    }

    if(T->verbose > 1)
    {
        printf("min_gini: %f\n", min_gini);
        printf("best_var: %d\n", best_var);
    }


// If it is impossible to split by any variable, make this a final node.
// This happens when all features under consideration are constant
    if(min_gini > 1)
    {
        root->class = most_common_class(CL, N, T->maxClass);
        assert(root->class <= T->maxClass);
        assert(root->left == NULL);
        assert(root->right == NULL);
        return root;
    }

    /* To further avoid issues with multiple identical feature
     * values we count the number of vectors left and right of the threshold
     */

    const float * V = X+N*best_var; // variable vector
    size_t _nL = 0;
    for(size_t kk = 0; kk<N; kk++)
    {
        V[kk] < best_th ? _nL ++ : 0;
    }

    if(_nL == 0 || _nL == N)
    {
        root->class = most_common_class(CL, N, T->maxClass);
        assert(root->class <= T->maxClass);
        assert(root->left == NULL);
        assert(root->right == NULL);
        return root;
    }

    size_t nL = _nL;
    size_t nR = N-_nL;



    if(T->feature_map == NULL)
    {
        root->var = best_var;
    } else {
        /* If we are working with a table with reduced number of
         * features, map this feature to the correct feature in the full set */
        root->var = T->feature_map[best_var];
    }



//    size_t nL = best_nL;
//    size_t nR = best_nR;

    root->gini = min_gini;
    root->th = best_th;
    float th = best_th;

    float lgini = best_lgini;
    float rgini = best_rgini;




    // If one size has zero elements we run into infinite recursion.
    if(nL == 0 || nR == 0)
    {
        fprintf(stderr, "Void split.... This is a bug. N = %zu, nL = %zu, nR = %zu\n",
                N, nL, nR);
        fprintf(stderr, "Quitting to avoid infinite recursion\n");
        exit(EXIT_FAILURE);
        getchar();
    }

    assert(nL + nR == N);


    /*
     * Split the table into two sub-tables,
     * one for each child node.
     */

    float * Ltab = NULL;
    float * Rtab = NULL;
    split_X(V, X, &Ltab, &Rtab, N, M, nL, nR, th);
    root->left = prf_tree_train_with_parent(T, root, nL, Ltab, lgini);
    free(Ltab);
    root->right = prf_tree_train_with_parent(T, root, nR, Rtab, rgini);
    free(Rtab);

    return root;
}

float squaresum(size_t * V, size_t T, size_t N)
{
    if(T == 0)
        return -1;
    float s = 0;
    for(size_t kk = 0; kk<N; kk++)
    {
        s+= pow((float) V[kk]/ (float) T, 2);
    }
    return s;
}

float * histogram(const float * C, const size_t N, const int nclass)
{
    float * H = calloc(nclass, sizeof(float));
    for(size_t kk = 0; kk<N; kk++)
    {
        H[(int) C[kk]]++;
    }
    return H;
}

float gini_from_C(const float * C, const size_t N, int const nclass)
{
    if(N < 2)
    {
        return 0;
    }
    float * H = histogram(C, N, nclass);
    float gini = 1;
    for(int cc = 0; cc<nclass; cc++)
    {
        //printf("%f ", H[cc]);
        gini -= pow(H[cc]/(float) N, 2);
    }
    //getchar();
    //printf("\n");
    free(H);
    return gini;
}

// Gini from two histograms
float gini_from_H(size_t nL, size_t * HL,
                   size_t nR, size_t * HR,
                   size_t nclass,
                   float * lgini, float * rgini)
{
    float gL = 1 - squaresum(HL, nL, nclass);
    float gR = 1 - squaresum(HR, nR, nclass);
    float N = nL + nR;
    lgini[0] = gL;
    rgini[0] = gR;
    return  (float) nL/N * gL + (float) nR/N * gR;
}

float gini_from_VC(float * VC, size_t N,
                    float * g, size_t * nLeft, size_t * nRight,
                    float * lgini, float * rgini)
{

    const int verbose = 0;

    if(verbose > 2)
    {
        for(size_t kk = 0 ; kk< N; kk++)
        {
            printf("g: %zu: %f, %f\n", kk, VC[3*kk], VC[3*kk+1]);
        }
    }
    assert(VC[0] <= VC[3]);


    /* Sweep threshold from left to right,
     * between points that does not have the same
     * value.
     */

    const size_t nclass = 3;
    size_t * HL = calloc(nclass, sizeof(size_t));
    size_t * HR = calloc(nclass, sizeof(size_t));

    for(size_t kk = 1; kk<3*N; kk+=3)
    {
        HR[(int) VC[kk]]++;
    }
    // The number of elements to the left and right of the threshold
    size_t nL = 0, nR = 0;
    for(size_t kk = 0; kk<nclass; kk++)
    {
        nR += HR[kk];
    }
    assert(nR == N);

    nLeft[0] = 0;
    nRight[0] = N;

    float minth = -INFINITY;
    float mingini = 2;
    float best_lgini = -1;
    float best_rgini = -1;

    assert(N > 1);

    float prev = VC[0];
    float next = VC[0];

    // Between element kk and element kk+1
    //printf("N = %zu\n", N);
    #ifndef NDEBUG
    size_t bestkk = 0;
    #endif

    for(size_t kk = 0; kk+1< N; kk++)
    {
        //printf("%zu ", kk);
        next = prev + 1e-6;
        if((kk+1) < N)
        {
            next = VC[3*(kk+1)];
        }
        int prevclass = VC[3*kk + 1];
        assert(prevclass >= 0);

        // printf("(%zu--%zu) prev: %f, next: %f\n", kk, kk+1, prev, next);

        nL++;
        HL[prevclass]++;
        nR--;
        HR[prevclass]--;

        float clgini = 0;
        float crgini = 0;
        float g = gini_from_H(nL, HL, nR, HR, nclass, &clgini, &crgini);
        //printf("(%f)", g);
        if(next > prev)
        {
            if(g < mingini)
            {
                #ifndef NDEBUG
                bestkk = kk;
                #endif
                mingini = g;
                minth = 0.5*(next + prev);
                nLeft[0] = nL;
                nRight[0] = nR;
                best_lgini = clgini;
                best_rgini = crgini;
            }
        }

        prev = next;
    }
//    printf("\n");
    lgini[0] = best_lgini;
    rgini[0] = best_rgini;
    #ifndef NDEBUG
    if(verbose > 2)
    {
        printf("nLeft: %zu, nRight: %zu, minth: %f (from kk = %zu)\n", nLeft[0], nRight[0], minth, bestkk);
    }
    #endif

    if(0 && mingini > 1)
    {
        printf("Unable to split\n");
        printf("mingini: %f\n", mingini);
        printf("nL: %zu, nR: %zu\n", nL, nR);
        getchar();
    }

    g[0] = mingini;

    free(HL);
    free(HR);

    return minth;
}
