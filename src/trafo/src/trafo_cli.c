#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <getopt.h>

#include "iris.h"
#include "trafo.h"
#include "trafo_util.h"
#include "ftab.h"

typedef uint32_t u32;
typedef double f64;

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


static int is_csv_file_name(const char * name)
{
    assert(name != NULL);
    size_t l = strlen(name);
    if(l < 4)
    {
        return 0;
    }

    if(strcmp(name + l-4, ".csv") == 0)
    {
        return 1;
    }
    return 0;
}

double * read_raw(const char * file, size_t N, size_t M)
{
    printf("Reading %s\n", file);
    printf("N = %zu, M = %zu\n", N, M);
    FILE * f = fopen(file, "rb");
    if(f == NULL)
    {
        printf("Can't open %s\n", file);
        exit(1);
    }
    double * X = calloc(M*N, sizeof(double));
    assert(X != NULL);
    size_t nread = fread(X, sizeof(double), M*N, f);
    if(nread != M*N)
    {
        printf("Didn't read the expected amount of data from %s\n", file);
        assert(nread == M*N);
        exit(1);
    }
    fclose(f);
    return X;
}

void compare_results(const u32 * R,
                     const u32 * P,
                     const size_t n)
{
    size_t neq = 0;
    for(size_t kk = 0; kk < n; kk++)
    {
        if(R[kk] == P[kk])
        {
            neq++;
        }
    }
    printf("%.2f %% correctly classified (%zu / %zu)\n",
           100.0 * (double) neq / (double) n,
           neq,
           n);
    return;
}

void test_iris_tree(int entropy)
{
    /* IRIS dataset -- builtin */
    if(entropy == 0)
    {
        print_section("IRIS -- single tree -- Gini Impurity");
    } else {
        print_section("IRIS -- single tree -- Entropy");
    }
    trafo_settings conf = {0};
    conf.F_col_major = iris_F;
    conf.label = iris_L;
    conf.n_tree = 1;
    conf.n_sample = iris_ns;
    conf.n_feature = iris_nf;
    conf.tree_f_sample = 1.0;
    conf.tree_n_feature = iris_nf;
    conf.min_samples_leaf = 1;
    conf.verbose = 1;
    conf.entropy = entropy;

    trf * F = trafo_fit(&conf);
    assert(F != NULL);

    /* Predict the training data and compare to training labels */
    u32 * P = trafo_predict(F, iris_F, NULL, iris_ns);
    compare_results(P, iris_L, iris_ns);

    printf("Feature importance*:\n");
    double * imp = trafo_importance(F);
    for(u32 ff = 0; ff < conf.n_feature; ff++)
    {
        printf("#%3d : %3.1f %%\n", ff, imp[ff]*100.0);
    }
    free(imp);

    printf("\n-> Saving to disk, reading from disk and comparing\n");
    trafo_save(F, "iris1.trf");
    trafo_free(F);
    trf * F2 = trafo_load("iris1.trf");
    if(F2 == NULL)
    {
        printf("Error reading from iris1.trf!\n");
    }

    u32 * P2 = trafo_predict(F2, iris_F, NULL, iris_ns);

    u32 neq = 0;
    for(size_t kk = 0; kk < iris_ns; kk++)
    {
        if(P[kk] == P2[kk])
        {
            neq++;
        }
    }

    printf("%u / %u predictions are equal\n", neq, iris_ns);
    if(neq != iris_ns)
    {
        assert(0);
    }

    free(P2);

    trafo_free(F2);

    free(P);

    return;
}

void test_iris_forest(int entropy)
{
    if(entropy == 0)
    {
        print_section("IRIS -- Forest -- Gini Impurity");
    } else {
        print_section("IRIS -- Forest -- Entropy");
    }

    /* IRIS dataset -- builtin */
    trafo_settings conf = {0};
    conf.F_col_major = iris_F;
    conf.label = iris_L;
    conf.n_tree = 20;
    conf.n_sample = iris_ns;
    conf.n_feature = iris_nf;
    conf.verbose = 1;
    conf.entropy = entropy;

    trf * F = trafo_fit(&conf);
    u32 * P = trafo_predict(F, iris_F, NULL, iris_ns);
    compare_results(P, iris_L, iris_ns);

    printf("Feature importance*:\n");
    double * imp = trafo_importance(F);
    for(u32 ff = 0; ff < conf.n_feature; ff++)
    {
        printf("#%3d : %3.1f %%\n", ff, imp[ff]*100.0);
    }
    free(imp);

    fflush(stdout);
    printf("\n-> Saving to disk, reading from disk and comparing\n");
    trafo_save(F, "iris20.trf");
    trafo_free(F);

    trf * F2 = trafo_load("iris20.trf");
    if(F2 == NULL)
    {
        printf("Error reading from iris20.trf!\n");
    }

    u32 * P2 = trafo_predict(F2, iris_F, NULL, iris_ns);
    trafo_free(F2);

    u32 neq = 0;
    for(size_t kk = 0; kk < iris_ns; kk++)
    {
        if(P[kk] == P2[kk])
        {
            neq++;
        }
    }
    free(P);
    free(P2);

    printf("%u / %u predictions are equal\n", neq, iris_ns);
    if(neq != iris_ns)
    {
        assert(0);
    }
}

typedef struct {
    char * file_train;
    char * file_predict;
    char * classifier_out;
    char * classifier_in;
    char * classcol;
    int classcol_id;
    u32 n_tree;
    u32 xfold;
    int builtin_tests;
    u32 min_leaf_size;
    float tree_f_sample;
    u32 tree_features;
    u32 verbose;
    int entropy;
} trafo_cli_settings;

void trafo_cli_settings_free(trafo_cli_settings * s)
{
    if(s == NULL)
    {
        return;
    }
    free(s->file_train);
    free(s->file_predict);
    free(s->classifier_out);
    free(s->classifier_in);
    free(s->classcol);
    free(s);
}

static void
trafo_cli_settings_show(trafo_cli_settings * s)
{
    printf("OMP_NUM_THREADS: %s\n", getenv("OMP_NUM_THREADS"));
    printf("file_train: %s\n", s->file_train);
    printf("file_predict: %s\n", s->file_predict);
    printf("classifier_out: %s\n", s->classifier_out);
    printf("classifier_in: %s\n", s->classifier_in);
    printf("classcol: %s\n", s->classcol);
    printf("n_tree=%u\n", s->n_tree);
    return;
}

static void usage(void)
{
    printf("Usage:\n");
    printf("--train file.tsv\n\t"
           "table to train on\n");
    printf("--cout file.trf\n\t"
           "Write the classifier to disk\n");
    printf("--ntree n\n\t"
           "number of trees in the forest\n");
    printf("--predict file.tsv\n\t"
           "table of point to classify\n");
    printf("--model file.trf\n\t"
           "classifer to use\n");
    printf("--classcol name\n\t"
           "specify the name of the column that contain the class/label\n");
    printf("--tree_samples n\n\t"
           "Fraction of samples per tree (to override default)\n");
    printf("--tree_features n\n\t"
           "Number of features per tree\n");
    printf("--min_leaf_size\n\t"
           "How small a node be before it is automatically turned into a leaf\n");
    printf("--verbose n\n\t"
           "Set verbosity level\n");
    printf("--entropy\n\t"
           "Split on entropy instead of Gini impurity\n");
    printf("--xfold n\n\t"
           "Perform n-fold cross validataion\n");
    printf("\n");
    printf("Example: 10-fold cross validation\n");
    printf("$ trafo --xfold 10 --train file.csv\n");
    return;
}

trafo_cli_settings * parse_cli(int argc, char ** argv)
{
    trafo_cli_settings * conf = calloc(1, sizeof(trafo_cli_settings));
    assert(conf != NULL);

    conf->classcol = strdup("class");
    conf->classcol_id = -1;
    conf->verbose = 1;
    conf->n_tree  = 1;

    struct option longopts[] = {
        {"train",         required_argument, NULL, 'c'},
        {"help",          no_argument,       NULL, 'h'},
        {"model",         required_argument, NULL, 'm'},
        {"classcol",      required_argument, NULL, 'n'},
        {"classcolID",    required_argument, NULL, 'i'},
        {"entropy",       no_argument,       NULL, 'e'},
        {"cout",          required_argument, NULL, 'w'},
        {"ntree",         required_argument, NULL, 't'},
        {"predict",       required_argument, NULL, 'p'},
        {"tree_samples",  required_argument, NULL, 's'},
        {"tree_features", required_argument, NULL, 'f'},
        {"min_leaf_size", required_argument, NULL, 'l'},
        {"verbose",       required_argument, NULL, 'v'},
        {"version",       no_argument,       NULL, 'V'},
        {"xfold",         required_argument, NULL, 'x'},
        {NULL, 0, NULL, 0}};

    int ch;
    while((ch = getopt_long(argc, argv,
                            "c:ef:hi:l:m:n:p:s:t:v:Vw:x:",
                            longopts, NULL)) != -1)
    {
        switch(ch){
        case 'c':
            free(conf->file_train);
            conf->file_train = strdup(optarg);
            break;
        case 'e':
            conf->entropy = 1;
            break;
        case 'f':
            conf->tree_features = atol(optarg);
            break;
        case 'h':
            usage();
            exit(EXIT_SUCCESS);
        case 'i':
            conf->classcol_id = atol(optarg);
            break;
        case 'l':
            conf->min_leaf_size = atol(optarg);
            break;
        case 'm':
            free(conf->classifier_in);
            conf->classifier_in = strdup(optarg);
            break;
        case 'n':
            free(conf->classcol);
            conf->classcol = strdup(optarg);
            break;
        case 'p':
            free(conf->file_predict);
            conf->file_predict = strdup(optarg);
            break;
        case 's':
            conf->tree_f_sample = atof(optarg);
            break;
        case 't':
            conf->n_tree = atol(optarg);
            break;
        case 'v':
            conf->verbose = atol(optarg);
            break;
        case 'V':
            printf("trafo_cli version %d.%d.%d\n",
                   TRAFO_VERSION_MAJOR, TRAFO_VERSION_MINOR, TRAFO_VERSION_PATCH);
            exit(EXIT_SUCCESS);
        case 'w':
            free(conf->classifier_out);
            conf->classifier_out = strdup(optarg);
            break;
        case 'x':
            conf->xfold = atol(optarg);
            break;

        default:
            usage();
            trafo_cli_settings_free(conf);
            return NULL;
        }
    }

    if( (conf->file_predict == NULL) && (conf->file_train == NULL) )
    {
        conf->builtin_tests = 1;
    }
    return conf;
}

ftab_t * ftab_from_dlm(const char * filename)
{
    const int verbose = 0;
    ftab_t * tab = NULL;
    if(is_csv_file_name(filename))
    {
        if(verbose > 1)
        {
            printf("Reading as csv file\n");
        }
        tab = ftab_from_csv(filename);
    } else {
        if(verbose > 1)
        {
            printf("Reading as tsv file\n");
        }
        tab = ftab_from_tsv(filename);
    }
    return tab;
}

void predict_file(trafo_cli_settings * conf)
{
    printf("Reading model from %s\n", conf->classifier_in);
    trf * T = trafo_load(conf->classifier_in);
    if(T == NULL)
    {
        fprintf(stderr, "Unable to load model\n");
        exit(EXIT_FAILURE);
    }

    printf("Reading data from %s\n", conf->file_predict);
    ftab_t * tab = ftab_from_dlm(conf->file_predict);
    if(tab == NULL)
    {
        fprintf(stderr, "Unable to load dataset\n");
        exit(EXIT_FAILURE);
    }
    printf("ncol: %zu, nrow: %zu\n", tab->ncol, tab->nrow);

    u32 n_feature = tab->ncol;

    int c_col = conf->classcol_id;
    if(c_col < 0)
    {
        c_col = ftab_get_col(tab, conf->classcol);
    }
    if(c_col < 0)
    {
        printf("Unable to find any column with the name \"%s\"\n",
               conf->classcol);
    } else {
        printf("Found feature column \"%s\"\n",
               conf->classcol);
        n_feature--;
    }

    u32 n_sample = tab->nrow;

    double * F = calloc(tab->ncol*tab->nrow, sizeof(double));
    assert(F != NULL);

    u32 * L = calloc(tab->nrow, sizeof(u32));
    assert(L != NULL);

    for(int row = 0; row < (int) tab->nrow; row ++)
    {
        u32 wcol = 0;
        for(int col = 0; col < (int) tab->ncol; col++)
        {
            if(col == c_col)
            {
                continue;
            }
            F[row*n_feature + wcol] = tab->T[row*tab->ncol + col];
            wcol++;
        }
    }

    if(c_col >= 0)
    {
        for(size_t row = 0; row < tab->nrow; row++)
        {
            L[row] = tab->T[row*tab->ncol + c_col];
        }
    }
    ftab_free(tab);


    struct timespec t0, t1;

    dw_gettime(&t0);
    u32 * P = trafo_predict(T, NULL, F, n_sample);
    dw_gettime(&t1);
    if(conf->verbose > 1)
    {
        printf("trafo: Prediction took %.6f s\n",
               timespec_diff(&t1, &t0));
    }


    if(L == 0)
    {
        printf("No ground truth to compare with\n");
    } else {
        u32 neq = 0;
        int max_label = 0;
        for(size_t kk = 0; kk < n_sample; kk++)
        {
            (int) L[kk] > max_label ? max_label = L[kk] : 0 ;
            (int) P[kk] > max_label ? max_label = P[kk] : 0 ;
            if(P[kk] == L[kk])
            {
                neq++;
            }
        }
        printf("%u / %u predicted correctly (%.2f %%)\n", neq, n_sample,
               100.0 * (double) neq / (double) n_sample);

        /* Construct a confusion matrix and display it */

        if(max_label < 10)
        {
            max_label++;
            u32 * CM = calloc((size_t )max_label* (size_t) max_label, sizeof(u32));
            assert(CM != NULL);
            for(size_t kk = 0; kk < n_sample ; kk++)
            {
                CM[L[kk] + max_label*P[kk]]++;
            }

            printf("\n");
            for(int kk = 0; kk < max_label; kk++)
            {
                for(int ll = 0; ll < max_label; ll++)
                {
                    printf("%4d ", CM[kk+max_label*ll]);
                }
                printf("\n");
            }
            printf("\n");

            free(CM);
        }
    }


    free(P);
    free(L);
    free(F);

    trafo_free(T);

    return;
}


void
subselect_features(u32 * S, u32 id,
                   const double * F,
                   u32 n_sample, u32 n_feature,
                   double * Ftrain, u32 * nFtrain,
                   double * Feval, u32 *nFeval)
{
    u32 n_train = 0;
    u32 n_pred = 0;
    for(size_t kk = 0; kk < n_sample; kk++)
    {
        if(S[kk] == id)
        {
            memcpy(Feval+n_pred*n_feature,
                   F+kk*n_feature,
                   n_feature*sizeof(double));
            n_pred++;
        } else {
            memcpy(Ftrain+n_train*n_feature,
                   F+kk*n_feature,
                   n_feature*sizeof(double));
            n_train++;
        }
    }
    *nFtrain = n_train;
    *nFeval = n_pred;
}

void
subselect_labels(u32 * S, u32 id,
                 const u32 * L,
                 u32 n_sample,
                 u32 * Ltrain, u32 * nLtrain,
                 u32 * Leval, u32 *nLeval)
{
    u32 n_train = 0;
    u32 n_pred = 0;
    for(size_t kk = 0; kk < n_sample; kk++)
    {
        if(S[kk] == id)
        {
            Leval[n_pred++] = L[kk];
        } else {

            Ltrain[n_train++] = L[kk];
        }
    }
    *nLtrain = n_train;
    *nLeval = n_pred;
    return;
}


void xfold(trafo_cli_settings * conf,
           double * F, u32 * L,
           u32 n_sample, u32 n_feature)
{
    char * header = calloc(80, 1);
    assert(header != NULL);
    sprintf(header, "%u-fold cross validation ntree: %u",
            conf->xfold, conf->n_tree);
    print_section(header);
    free(header);

    double * acc = calloc(conf->xfold, sizeof(double));
    assert(acc != NULL);
    u32 * S = calloc(n_sample, sizeof(u32));
    assert(S != NULL);

    for(size_t kk = 0; kk < n_sample; kk++)
    {
        S[kk] = kk % conf->xfold;
    }

    double * Ftrain = calloc(n_sample*n_feature, sizeof(double));
    assert(Ftrain != NULL);
    double * Feval = calloc(n_sample*n_feature, sizeof(double));
    assert(Feval != NULL);
    u32 * Ltrain = calloc(n_sample, sizeof(u32));
    assert(Ltrain != NULL);
    u32 * Leval = calloc(n_sample, sizeof(u32));
    assert(Leval != NULL);

    //double * Ft = transpose_f64(F, n_feature, n_sample);
    for(u32 ff = 0; ff < conf->xfold; ff++)
    {
        if(conf->verbose > 1)
        {
            printf("Selection %u/%u\n", ff+1, conf->xfold);
        }
        u32 nFtrain;
        u32 nFeval;
        subselect_features(S, ff,
                           F,
                           n_sample, n_feature,
                           Ftrain, &nFtrain,
                           Feval, &nFeval);

        u32 nLtrain;
        u32 nLeval;
        subselect_labels(S, ff,
                         L,
                         n_sample,
                         Ltrain, &nLtrain,
                         Leval, &nLeval);
        if(conf->verbose > 1)
        {
            printf("Using %u for training and %u for testing\n",
                   nLtrain, nLeval);
        }

        trafo_settings Tconf = {0};
        Tconf.F_row_major = Ftrain;
        Tconf.label = Ltrain;
        Tconf.n_tree = conf->n_tree;
        Tconf.n_sample = nLtrain;
        Tconf.n_feature = n_feature;
        Tconf.verbose = 0;
        Tconf.tree_n_feature = conf->tree_features;
        Tconf.tree_f_sample = conf->tree_f_sample;
        Tconf.min_samples_leaf = conf->min_leaf_size;
        Tconf.entropy = conf->entropy;

        /* Train */
        trf * T = trafo_fit(&Tconf);
        /* Predict the left out */
        u32 * P = trafo_predict(T, NULL, Feval, nFeval);
        for(size_t kk = 0; kk < nLeval; kk++)
        {
            //printf("%u -- %u\n", P[kk], Leval[kk]);
            acc[ff] += (P[kk] == Leval[kk]);
        }
        acc[ff] /= (double) nLeval;
        acc[ff]*=100.0;
        //printf("Accuracy: %f\n", acc[ff]);
        free(P);
        trafo_free(T);
    }

    //free(Ft);
    double macc = 0;
    for(size_t kk = 0; kk < conf->xfold; kk++)
    {
        macc += acc[kk];
    }

    macc = macc / conf->xfold;
    printf("Mean pred acc: %.2f %%\n", macc);
    free(Feval);
    free(Ftrain);
    free(Leval);
    free(Ltrain);
    free(S);
    free(acc);
    return;
}

void train_file(trafo_cli_settings * conf)
{
    printf("Reading from %s\n", conf->file_train);
    ftab_t * tab = ftab_from_dlm(conf->file_train);

    if(tab == NULL)
    {
        printf("Error reading input file\n");
        return;
    }
    printf("columns: %zu, data rows: %zu\n", tab->ncol, tab->nrow);
    if(conf->verbose > 1)
    {
        printf("Columns:\n");
        for(size_t kk = 0; kk < tab->ncol; kk++)
        {
            printf("%zu: '%s'\n", kk, tab->colnames[kk]);
        }
    }

    u32 n_feature = tab->ncol;

    int c_col = conf->classcol_id;
    if(c_col < 0)
    {
        c_col = ftab_get_col(tab, conf->classcol);
    }
    if(c_col < 0)
    {
        printf("Unable to find any column with the name \"%s\"\n"
               "Check the input data or use another column for class"
               "with the --clascol argument\n",
               conf->classcol);
        ftab_free(tab);
        return;
    } else {
        printf("Found feature column \"%s\"\n",
               conf->classcol);
        n_feature--;
    }

    u32 n_sample = tab->nrow;

    double * F = calloc(tab->ncol*tab->nrow, sizeof(double));
    assert(F != NULL);

    u32 * L = calloc(tab->nrow, sizeof(u32));
    assert(L != NULL);

    for(int row = 0; row < (int) tab->nrow; row ++)
    {
        u32 wcol = 0;
        for(int col = 0; col < (int) tab->ncol; col++)
        {
            if(col == c_col)
            {
                continue;
            }
            F[row*n_feature + wcol] = tab->T[row*tab->ncol + col];
            wcol++;
        }
    }

    if(c_col >= 0)
    {
        for(size_t row = 0; row < tab->nrow; row++)
        {
            L[row] = tab->T[row*tab->ncol + c_col];
        }
    }
    ftab_free(tab);

    if(conf->xfold > 0)
    {
        xfold(conf, F, L, n_sample, n_feature);
        free(L);
        free(F);
        return;
    }

    trafo_settings Tconf = {0};
    Tconf.F_row_major = F;
    Tconf.label = L;
    Tconf.n_tree = conf->n_tree;
    Tconf.n_sample = n_sample;
    Tconf.n_feature = n_feature;
    Tconf.verbose = 1;
    Tconf.tree_n_feature = conf->tree_features;
    Tconf.tree_f_sample = conf->tree_f_sample;
    Tconf.min_samples_leaf = conf->min_leaf_size;
    Tconf.entropy = conf->entropy;

    if(conf->n_tree == 1)
    {
        if(conf->tree_f_sample == 0)
        {
            Tconf.tree_f_sample = 1;
        }
        if(conf->tree_features == 0)
        {
            Tconf.tree_n_feature = n_feature;
        }
        if(conf->min_leaf_size == 0)
        {
            Tconf.min_samples_leaf = 1;
        }
    }

    struct timespec t0, t1;
    size_t mem0, mem1, mem_tmp;
    get_peakMemoryKB(&mem_tmp, &mem0);
    dw_gettime(&t0);
    trf * T = trafo_fit(&Tconf);
    dw_gettime(&t1);
    get_peakMemoryKB(&mem_tmp, &mem1);
    if(conf->verbose > 1)
    {
        printf("deltaRSS: %zu\n", mem1-mem0);
    }
    if(T == NULL)
    {
        fprintf(stderr, "Fatal error: trf_fit returned NULL\n");
        exit(EXIT_FAILURE);
    }

    printf("trafo: Forest training took %.6f s\n",
           timespec_diff(&t1, &t0));

    double * imp = trafo_importance(T);
    if(imp != NULL)
    {

        printf("Feature importance*:\n");

        for(u32 ff = 0; ff < Tconf.n_feature; ff++)
        {
            printf("#%3d : %3.1f %%\n", ff, imp[ff]*100.0);
        }
        free(imp);
    }

    if(conf->classifier_out)
    {
        printf("Writing to %s\n", conf->classifier_out);
        trafo_save(T, conf->classifier_out);
    }

    free(L);
    free(F);

    trafo_free(T);

    return;
}

int main(int argc, char ** argv)
{
#ifndef NDEBUG
    printf("Compiled with DEBUG asserts. Speed is compromised.\n\n");
    printf("sizeof(size_t) == %zu\n", sizeof(size_t));
#endif


    trafo_cli_settings * conf = parse_cli(argc, argv);
    if(conf == NULL)
    {
        return EXIT_FAILURE;
    }
    if(conf->verbose > 1)
    {
        trafo_cli_settings_show(conf);
    }

    if(conf->builtin_tests)
    {
        test_iris_tree(0);
        test_iris_tree(1);
        print_peak_memory();
        test_iris_forest(0);
        test_iris_forest(1);
        print_peak_memory();
        trafo_ut();
        trafo_cli_settings_free(conf);
        return EXIT_SUCCESS;
    }

    if(conf->file_train)
    {
        train_file(conf);
    }

    if(conf->file_predict)
    {
        predict_file(conf);
    }

    trafo_cli_settings_free(conf);

    print_peak_memory();
    return EXIT_SUCCESS;
}
