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
    int nrows = 10;
    if(nrows > N)
    {
        nrows = N;
    }
    for(int kk = 0; kk<10; kk++)
    {
        printf(" %d, ", kk);
        for(int rr = 0; rr<M; rr++)
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
        for(int ii = 0; ii < (M-1); ii++)
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

    if(verbose > 0)
    {
        //size_t nnodes = prf_tree_enumerate(C);
        //printf("The tree has %zu nodes\n", nnodes);
    }

    prf_tree_free(T);
    }

int main(int argc, char ** argv)
{
    struct timespec tictoc_start, tictoc_end;

    char file[] = "../prototype/X150001x14.float"; size_t N = 150001; size_t M = 14;

    //char file[] = "../prototype/X15001x3.float"; size_t N = 15001; size_t M = 3;
    //char file[] = "../prototype/X1500001x15.float"; size_t N = 1500001; size_t M = 15;
    //char file [] = "/home/erikw/code/pixelClassifier/ieg728_60x_50_iter_test/X22521x51.float"; size_t N = 22521; size_t M = 51;
    //char file [] = "/home/erikw/code/pixelClassifier/ieg728_60x_50_iter_test/X2816x10.float"; size_t N = 2816; size_t M = 10;

    print_section("Input data");
    float * X = read_raw(file, N, M);
    preview_data(X, N, M);

    print_section("testing PrfTree");
    test_one_tree(X, N, M);
    printf("Peak Memory %zu kb\n", get_peakMemoryKB());

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
