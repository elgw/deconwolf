#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>


#include "prf_tree.h"
#include "prf_forest.h"



void test_forest(float * X, const size_t N, const size_t M)
{
    struct timespec tictoc_start, tictoc_end;
    size_t ntrees = 200;
    printf("Creating a forest with %zu trees\n", ntrees);
    clock_gettime(CLOCK_REALTIME, &tictoc_start);
    PrfForest * F = prf_forest_new(ntrees);
    F->nthreads = 8;

    prf_forest_train(F, X, N, M);
    clock_gettime(CLOCK_REALTIME, &tictoc_end);
    printf(" Forest training took %f s\n", timespec_diff(&tictoc_end, &tictoc_start));

    print_todo("validate multi-threaded implementation");

    clock_gettime(CLOCK_REALTIME, &tictoc_start);

    int * class = prf_forest_classify_table(F, X, N, M);

    size_t ncorrect = 0;
    for(size_t kk = 0; kk<N; kk++)
    {
        if( class[kk] == X[kk + (M-1)*N])
        {
            ncorrect ++;
        }
    }

    printf("%zu / %zu correctly classified\n", ncorrect, N);

    clock_gettime(CLOCK_REALTIME, &tictoc_end);
    printf(" Forest classification took %f s\n", timespec_diff(&tictoc_end, &tictoc_start));

    printf("Validating multi-threaded implementation... ");
    F->nthreads = 1;
    int * class_single = prf_forest_classify_table(F, X, N, M);
    int err = 0;
    for(size_t kk = 0; kk<N; kk++)
    {

        if(class[kk] != class_single[kk])
        {
            err = 1;
            printf("%zu: %d %d \n", kk, class_single[kk], class[kk]);
        }

    }
    if(err)
    {
        print_error("Results differs between single and multiple-core implementation");
        exit(1);
    }
    print_ok();

    free(class_single);
    free(class);
    prf_forest_free(F);
}

float * read_raw(char * file, size_t N, size_t M)
{
    printf("Reading %s\n", file);
    printf("N = %zu, M = %zu\n", N, M);
    FILE * f = fopen(file, "r");
    if(f == NULL)
    {
        printf("Cant open %s\n", file);
        exit(1);
    }
    float * X = malloc(M*N*sizeof(float));
    size_t nread = fread(X, sizeof(float), M*N, f);
    if(nread != M*N)
    {
        printf("Didn't read the expected amount of data from %s\n", file);
        assert(nread == M*N);
        exit(1);
    }
    fclose(f);
    return X;
}

void preview_data(float * X, size_t N, size_t M)
{
    printf("Input data:\n");
    size_t nrows = 10;
    if(nrows > N)
    {
        nrows = N;
    }
    for(int kk = 0; kk<10; kk++)
    {
        printf(" %d :  ", kk);
        for(size_t rr = 0; rr<M; rr++)
        {
            printf("% f ", X[kk + rr*N]);
        }
        printf("\n");
    }
    printf(" ...\n");
}


    void test_one_tree(float * X, size_t N, size_t M)
    {

        int verbose = 0;
        struct timespec tictoc_start, tictoc_end;

    PrfTree * T = prf_tree_new();
    T->verbose = 1;
    T->n_features = M-1;
    T->n_samples = N;
    T->min_size = 1; // 1 for classification tree, 10 for random forest?
    printf("Settings:\n");
    prf_tree_show(stdout, T);
    clock_gettime(CLOCK_REALTIME, &tictoc_start);
    prf_tree_train(T, X);


    clock_gettime(CLOCK_REALTIME, &tictoc_end);
    printf(" Training took %f s\n", timespec_diff(&tictoc_end, &tictoc_start));

    clock_gettime(CLOCK_REALTIME, &tictoc_start);


    size_t ncorrect = 0;
    float * v = malloc((M-1)*sizeof(float));
    for(size_t kk = 0; kk<N; kk++)
    {
        for(size_t ii = 0; ii < (M-1); ii++)
        {
            v[ii] = X[kk + N*ii];
        }

        int class = prf_tree_classify(T, v);
        if(verbose > 0)
        {
        if(kk<10){
            printf("X = "); vector_show(v, (M-1));
        printf("row %zu : x = %f class = %f classification = %d\n",
               kk, X[kk], X[kk+(M-1)*N], class);
        }
        }
        if(class == X[kk+(M-1)*N])
        {
            ncorrect++;
        }

    }
    free(v);
    printf("%zu/%zu correctly classified\n", ncorrect, N);
    clock_gettime(CLOCK_REALTIME, &tictoc_end);
    printf(" Classification took %f s\n", timespec_diff(&tictoc_end, &tictoc_start));


    printf("-> Classifying using PrfTreeTable\n");
    PrfTreeTable * TT = prf_tree_to_tree_table(T);
    clock_gettime(CLOCK_REALTIME, &tictoc_start);
    ncorrect = 0;
    v = malloc((M-1)*sizeof(float));
    for(size_t kk = 0; kk<N; kk++)
    {
        for(size_t ii = 0; ii < (M-1); ii++)
        {
            v[ii] = X[kk + N*ii];
        }

        int class = prf_tree_table_classify(TT, v);
        if(verbose > 0)
        {
            if(kk<10){
                printf("X = "); vector_show(v, (M-1));
                printf("row %zu : x = %f class = %f classification = %d\n",
                       kk, X[kk], X[kk+(M-1)*N], class);
            }
        }
        if(class == X[kk+(M-1)*N])
        {
            ncorrect++;
        }

    }
    free(v);

    printf("%zu/%zu correctly classified\n", ncorrect, N);
    clock_gettime(CLOCK_REALTIME, &tictoc_end);
    printf(" Classification took %f s\n", timespec_diff(&tictoc_end, &tictoc_start));

    prf_tree_table_free(TT);

    prf_tree_free(T);
    }

int main(int argc, char ** argv)
{
    struct timespec tictoc_start, tictoc_end;

    float * X = NULL;
    size_t M = 7;
    size_t N = 1000000;
    if(argc > 1)
    {
        N = atoi(argv[1]);
    }

    /* Generate some random test data */
    if(M > 0)
    {
        X = malloc(M*N*sizeof(float));
        for(size_t kk = 0; kk<N*(M-1); kk++)
        {
            X[kk] = (float) rand() / (float) RAND_MAX;
        }
        /* Two classes */
        for(size_t kk = N*(M-1); kk< N*M; kk++)
        {
            X[kk] = 1 + (kk % 2);
        }
    }

    if(X == NULL)
    {
        char file[] = "../prototype/X150001x14.float";
        N = 150001;
        M = 14;
        X = read_raw(file, N, M);
    }


    print_section("Input data");

    preview_data(X, N, M);

    print_section("testing PrfTree");
    test_one_tree(X, N, M);
    printf("Peak Memory %zu kb\n", get_peakMemoryKB());

    if(argc > 2 && atoi(argv[2]) == 1)
    {
        free(X);
        return EXIT_SUCCESS;
    }

    print_todo("Cross-validation");

    clock_gettime(CLOCK_REALTIME, &tictoc_start);

    print_section("testing PrfForest");
    test_forest(X, N, M);
    printf("Peak Memory %zu kb\n", get_peakMemoryKB());

    clock_gettime(CLOCK_REALTIME, &tictoc_end);
    printf(" Forest test took %f s\n", timespec_diff(&tictoc_end, &tictoc_start));

    print_section("PrfForest: Cross validation test");
    prf_forest_cross_validate_k(NULL, X, N, M, 10);

    print_section("PrfForest: Feature importance test:");
    prf_forest_feature_importance(NULL, X, N, M, NULL);

    free(X);

    return 0;
}
