#include "prf_forest.h"

void prf_forest_alloc_trees(PrfForest * F)
{
    assert(F != NULL);
    if(F->ntrees == 0)
    {
        F->trees = NULL;
        return;
    }
    assert(F->trees == NULL);

    F->trees = malloc(F->ntrees*sizeof(PrfTree*));
    for(size_t kk = 0; kk < F->ntrees; kk++)
    {
        F->trees[kk] = NULL;
    }
    return;
}

PrfForest * prf_forest_new(size_t ntrees)
{
    PrfForest * F = malloc(sizeof(PrfForest));
    F->trees = NULL;
    F->ntrees = ntrees;
    prf_forest_alloc_trees(F);

    F->samples_per_tree = PRF_AUTOMATIC;
    F->features_per_tree = PRF_AUTOMATIC;
    F->nthreads = 16;
    F->min_node_size = 200;
    F->verbose = 1;
// Largest integer used for a class, set by prf_forest_train
    F->maxClass = -1;
    F->H = NULL;
    F->nH_allocated = 0;
    return F;
}

void prf_forest_free_trees(PrfForest * F)
{
    if(F->trees != NULL)
    {
        for(size_t tt = 0; tt < F->ntrees; tt++)
        {
            prf_tree_free(F->trees[tt]);
        }
        free(F->trees);
    }
    return;
}

void prf_forest_free(PrfForest * F)
{
    if(F == NULL)
    {
        print_warning("Trying to free NULL PrfForest");
        return;
    }

    prf_forest_free_H(F);
    prf_forest_free_trees(F);

    free(F);
}

float * select_m_from_table(const float * T,
                             const size_t N, const size_t M, const size_t m,
                             uint32_t ** feature_map)
{
    /* Select m random features from an NxN table with M-1 features */

    float * TS = malloc(N*(m+1)*sizeof(float));
    if(m + 1 == M)
    {
        memcpy(TS, T, M*N*sizeof(float));
        feature_map[0] = NULL;
        return TS;
    }

    uint32_t * fmap = malloc(m*sizeof(uint32_t));
    feature_map[0] = fmap;

    assert(m < M);
    uint8_t * use = get_random_selection(M-1, m);
    use[M-1] = 1;

    size_t wpos = 0;
    size_t fpos = 0;
    for(size_t mm = 0; mm < M; mm++)
    {
        if(use[mm] && mm+1 < M)
        {
            fmap[fpos++] = mm;
        }
        if(use[mm])
        {
            for(size_t nn = 0; nn < N; nn++)
            {
                TS[wpos++] = T[nn + mm*N];
            }
        }
    }
    free(use);
    return TS;
}

float * select_n_from_table(const float * T, const size_t N, const size_t M, const size_t n)
{
    // Col-Major
    // Select n out of N samples from the table T
    // with size NxM
    // I.e., return a nxM table

    assert(n<=N);
    if(n == N)
    {
        float * TS = malloc(M*N*sizeof(float));
        memcpy(TS, T, M*N*sizeof(size_t));
        return TS;
    }

    float * TS = malloc(n*M*sizeof(float));
    uint8_t * use = get_random_selection(N, n);

    size_t ww = 0; // write to this row
    for(size_t rr = 0; rr<N; rr++)
    {
        if(use[rr])
        {
            for(size_t mm = 0; mm < M; mm++)
            {
                assert(rr + mm*N < M*N);
                TS[ww + mm*n] = T[rr + mm*N];
            }
            ww++;
        }
    }
    free(use);
    return TS;
}


void prf_forest_train_tree(PrfForest * F, const float * X,
                           const size_t N, const size_t M, int Tid)
{
    uint32_t * feature_map = NULL;

    size_t MS = F->features_per_tree + 1;
    if(F->verbose > 1)
    {
        printf("Subselecting %d / %d features\n", F->features_per_tree, (int) M-1);
    }
    // Subselect features
    float * XSf = select_m_from_table(X, N, M, F->features_per_tree,
                                       &feature_map);

    // Subselect samples
    float * XSfSn = select_n_from_table(XSf, N, MS, F->samples_per_tree);
    free(XSf);

    PrfTree * T = prf_tree_new();
    assert(F->trees != NULL);
    F->trees[Tid] = T;
    T->n_samples = F->samples_per_tree;
    T->n_features = F->features_per_tree;
    T->min_size = F->min_node_size;
    T->feature_map = feature_map;
    prf_tree_train(T, XSfSn);

    free(XSfSn);
}

struct prf_forest_threaddata {
    PrfForest * F;
    const float * X;
    size_t N;
    size_t M;
    int thread;
    // For classification
    int * C;
};

void * prf_forest_train_thread(void * p)
{
    struct prf_forest_threaddata * D = (struct prf_forest_threaddata*) p;
    for(size_t tt = D->thread; tt < D->F->ntrees; tt += D->F->nthreads)
    {
        prf_forest_train_tree(D->F, D->X, D->N, D->M, tt);
    }
    return NULL;
}

void prf_forest_print(FILE * f, PrfForest * F)
{
    fprintf(f, "-> PrfForest with settings:\n");
    fprintf(f, "\t# trees: %zu\n", F->ntrees);
    fprintf(f, "\tthreads: %d\n", F->nthreads);
    fprintf(f, "\tMin node size: %d\n", F->min_node_size);
    fprintf(f, "\tFeatures per tree: ");
    if(F->features_per_tree == PRF_AUTOMATIC)
    {
        printf("automatic\n");
    } else {
        printf("%d\n", F->features_per_tree);
    }
    fprintf(f, "\tSamples per tree: ");
    if(F->samples_per_tree == PRF_AUTOMATIC)
    {
        printf("automatic\n");
    } else {
        printf("%d\n", F->samples_per_tree);
    }

    return;
}

/* X is of size NxM i.e. M includes the labels */
int prf_forest_train(PrfForest * F, const float * X, const size_t N, const size_t M)
{
    printf("Training a forest with %zu vectors of length %zu\n", N, M);
    assert(N > M);

    if(F->samples_per_tree == PRF_AUTOMATIC)
    {
        // Also known as bootstrapping
        F->samples_per_tree = round( 0.5 * (float) N);
        // F->samples_per_tree = round( (1.0 - 1.0/exp(1.0)) * (float) N);
    }
    if(F->features_per_tree == PRF_AUTOMATIC)
    {
        F->features_per_tree = ceil(sqrt(M-1));
    }

    if(F->trees == NULL)
    {
        F->trees = malloc(F->ntrees * sizeof(PrfTree));
    }

    if(F->verbose > 0)
    {
        prf_forest_print(stdout, F);
    }


    const float * L = X + (M-1)*N;
    float maxClass = L[0];
    float minClass = L[0];
    for(size_t kk = 0; kk<N; kk++)
    {
        if(L[kk] != (int) L[kk])
        {
            fprintf(stderr, "prf_forest_train: Error -- class of element "
                    "%zu is %f which isn't an integer\n", kk, L[kk]);
            return EXIT_FAILURE;
        }
        L[kk] > maxClass ? maxClass = L[kk] : 0;
        L[kk] < minClass ? minClass = L[kk] : 0;
    }
    printf("Classes: [%f, %f]\n", minClass, maxClass);
    F->maxClass = (int64_t) maxClass;

    if(F->nthreads > 1)
    {
        pthread_t * threads = malloc(F->nthreads*sizeof(pthread_t));
        struct prf_forest_threaddata * threaddata
            = malloc(F->nthreads*sizeof(struct prf_forest_threaddata));

        for(int tt = 0; tt < F->nthreads; tt++)
        {
            threaddata[tt].X = X;
            threaddata[tt].F = F;
            threaddata[tt].N = N;
            threaddata[tt].M = M;
            threaddata[tt].thread = tt;
            pthread_create(threads+tt, NULL,
                           prf_forest_train_thread, threaddata+tt);
        }

        for(int tt = 0; tt < F->nthreads; tt++)
        {
            pthread_join(threads[tt], NULL);
        }
        free(threaddata);
        free(threads);
    } else
    {
        if(F->verbose > 0)
        {
            printf("\n");
        }
        for(size_t tt = 0; tt<F->ntrees; tt++)
        {
            if(F->verbose > 0)
            {
                printf("\r                                      ");
                printf("\rTraining tree %zu / %zu ", tt+1, F->ntrees);
                fflush(stdout);
            }

            prf_forest_train_tree(F, X, N, M, tt);
        }
        printf("\n");
    }
    if(F->verbose > 0)
    {
        printf("Done training\n");
        fflush(stdout);
    }
    return EXIT_SUCCESS;
}

void prf_forest_copy_settings(PrfForest * F_target, const PrfForest * F_source)
{
    if(F_target->trees != NULL)
    {
        for(size_t kk = 0; kk < F_target->ntrees; kk++)
        {
            free(F_target->trees[kk]);
        }
        free(F_target->trees);
    }

    memcpy(F_target, F_source, sizeof(PrfForest));
    F_target->trees = NULL;
}

void prf_forest_free_H(PrfForest * F)
{
    if(F->H != NULL)
    {
        for(size_t hh = 0; hh < F->nH_allocated; hh++)
        {
            free(F->H[hh]);
        }
        free(F->H);
        F->H = NULL;
    }
}

void prf_forest_init_H(PrfForest * F)
{
    if((int) F->nH_allocated == F->nthreads)
    {
        return;
    }

    // Free if already existing
    prf_forest_free_H(F);

    assert(F->maxClass > 0);
    F->H = malloc(F->nthreads * sizeof(size_t*));
    for(int hh = 0; hh < F->nthreads; hh++)
    {
        F->H[hh] = malloc((F->maxClass + 1)*sizeof(size_t));
    }
    F->nH_allocated = F->nthreads;
}

int prf_forest_classify_thread(PrfForest * F, const float * X, int thread)
{
    const size_t ntrees = F->ntrees;
    prf_forest_init_H(F);

    const size_t maxlabel = F->maxClass;
    assert(maxlabel > 0);

    // Single thread:
    size_t * C = F->H[thread];
    for(size_t kk = 0; kk < maxlabel+1; kk++)
        C[kk] = 0;

    for(size_t tt = 0; tt < ntrees; tt++)
    {
        int c = prf_tree_classify(F->trees[tt], X);
        assert(c <= (int) maxlabel);
        C[c]++;
    }
    size_t class = vector_size_t_argmax(C, maxlabel+1);

    return (int) class;
}

int prf_forest_classify(PrfForest * F, const float * X)
{
    int thread = 0;
    return prf_forest_classify_thread(F, X, thread);
}

void * prf_forest_classify_table_thread(void * p)
{
    struct prf_forest_threaddata * D = (struct prf_forest_threaddata*) p;
    int * C = D->C;
    size_t M = D->M;
    size_t N = D->N;
    const float * X = D->X;
    PrfForest * F = D->F;
    int thread = D->thread;

    float * XF = malloc(M*sizeof(float));
    for(size_t kk = D->thread; kk < D->N; kk += D->F->nthreads)
    {
        copy_row(XF, 1, X+kk, N, M);
        C[kk] = prf_forest_classify_thread(F, XF, thread);
    }
    free(XF);
    return NULL;
}

int * prf_forest_classify_table(PrfForest * F, float * X, size_t N, size_t M)
{
    int * C = malloc(N*sizeof(int));

    if(F->nthreads == 1)
    {
        float * XF = malloc(M*sizeof(float));
        for(size_t kk = 0; kk<N; kk++)
        {
            copy_row(XF, 1, X+kk, N, M);
            C[kk] = prf_forest_classify(F, XF);
        }
        free(XF);
        return C;
    } else {
        prf_forest_init_H(F);
        pthread_t * threads = malloc(F->nthreads*sizeof(pthread_t));
        struct prf_forest_threaddata * threaddata
            = malloc(F->nthreads*sizeof(struct prf_forest_threaddata));

        for(int tt = 0; tt < F->nthreads; tt++)
        {
            threaddata[tt].X = X;
            threaddata[tt].F = F;
            threaddata[tt].N = N;
            threaddata[tt].M = M;
            threaddata[tt].C = C;
            threaddata[tt].thread = tt;
            pthread_create(threads+tt, NULL,
                           prf_forest_classify_table_thread, threaddata+tt);
        }

        for(int tt = 0; tt < F->nthreads; tt++)
        {
            pthread_join(threads[tt], NULL);
        }
        free(threaddata);
        free(threads);

        return C;
    }
}

float prf_forest_cross_validate_k(const PrfForest * F_template,
                                   const float * X, const size_t N, const size_t M,
                                   const int K)
{
    int verbose = 1;
    int ntrees = 10;
    PrfForest * F0 = prf_forest_new(ntrees);
    if(F_template != NULL)
    {
        prf_forest_copy_settings(F0, F_template);
        verbose = F_template->verbose;
    } else {
        printf("No settings supplied to prf_forest_cross_validate, using:\n");
        prf_forest_print(stdout, F0);
    }

    size_t total = 0;
    size_t total_correct = 0;

    if(verbose > 0)
    {
        printf("\n");
    }

    // K-fold cross-validation
    for(int k = 0; k < K; k++)
    {
        // Define training set
        float * Tr = NULL;
        size_t nTr = 0;
        // Define validation set
        float * Va = NULL;
        size_t nVa = 0;

        get_subset_k(X, N, M, K, k, &Tr, &nTr, &Va, &nVa);

        PrfForest * F = prf_forest_new(ntrees);

        prf_forest_copy_settings(F, F0);
        F->nthreads = 8;
        F->verbose = 0;

        prf_forest_train(F, Tr, nTr, M);

        // Validate
        float * XF = malloc(M*sizeof(float));
        size_t correct = 0;
        for(size_t kk = 0; kk<nVa; kk++)
        {
            copy_row(XF, 1, Va+kk, nVa, M);
            int class = prf_forest_classify(F, XF);
            int true_class = Va[(M-1)*nVa + kk];
            if(class == true_class)
            {
                correct++;
            }
        }
        free(XF);

        if(verbose > 0)
        {
            printf("partition %2d/%2d: %zu / %zu correct\n", k+1, K, correct, nVa);
        }
        total += nVa;
        total_correct += correct;

        prf_forest_free(F);
        free(Va);
        free(Tr);
    }
    prf_forest_free(F0);
    float pcorrect = 100.0*(float) total_correct / (float) total;
    if(verbose > 0)
    {
        printf("\n");
        printf("\tAverage classification:\n");
        printf("\t%zu / %zu correct, i.e., %.2f %%\n", total_correct, total,
               pcorrect);
        printf("\n");
    }
    return pcorrect;
}


void prf_forest_feature_importance(const PrfForest * _F,
                                   const float * X, const size_t N, const size_t M,
                                   float ** _FI)
{
    float * FI = malloc(M*sizeof(float));
    if(_FI != NULL)
    {
        _FI[0] = FI;
    }

    PrfForest * F = prf_forest_new(50);
    if(_F == NULL)
    {
        printf("No settings supplied to prf_forest_feature_importance, using:\n");
        prf_forest_print(stdout, F);
        F->verbose = 0;
    } else {
        prf_forest_copy_settings(F, _F);
    }

    F->verbose = 0;
    printf("\n");
    FI[M-1] = prf_forest_cross_validate_k(F, X, N, M, 10);
    for(int ff = 0; ff+1 < (int) M; ff++)
    {
        printf("\rTesting without feature %d / %ld...", ff+1, M-1);
        fflush(stdout);
        float * XD = malloc(M*N*sizeof(float));
        memcpy(XD, X, N*M*sizeof(float));
        scramble_feature(XD, N, M, ff);
        FI[ff] = prf_forest_cross_validate_k(F, XD, N, M, 10);
        free(XD);
    }
    printf("\n");
    printf("Feature importance:\n");
    printf("Excluded\tClassification\tImportance\n");
    printf("   none  \t%.2f        \t% 5.2f %%\n", FI[M-1], 100.0);
    for(int ff = 0; ff+1 < (int) M; ff++)
    {
        printf("%5d     \t%.2f        \t% 5.2f %%\n", ff+1, FI[ff], FI[M-1]-FI[ff]);
    }
    printf("\n");


    if(_FI == NULL)
    {
        free(FI);
    }

    prf_forest_free(F);

    return;
}
