#ifndef __prf_tree_h__
#define __prf_tree_h__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <malloc.h>
#include "qsort.h"


#define PRF_AUTOMATIC -1
#define PRF_NODE_NO_CLASS -1

typedef struct PrfNode{
    struct PrfNode * parent;
    struct PrfNode * left;
    struct PrfNode * right;

    int class; // -1 unless a final node.

    int var; // Variable, set to -1 if not decided yet
    float th; // Threshold

    size_t N; // Number of rows to train on
    float gini; // Gini (for the split)
    // unique id for this node
    // Set by prf_tree_enumerate
    size_t id;
} PrfNode;



typedef struct {
    int min_size;
    int verbose;
    int maxClass; // max class number
    // Number of features to consider when splitting a node.
    int max_features;

    size_t n_features; // Number of feature columns
    size_t n_samples;

    PrfNode * nodes_root; // The root node
    PrfNode * node_pointer; // Where to insert the next node

    // Buffers
    int * varlist;
    uint32_t * best_idx;
    float * VC; // 2*N, Variable and Class interlaced for Gini calculations
    int * randperm_buff; // 2*N
    uint32_t * feature_map;
    size_t n_nodes; // Number of used nodes
    size_t n_nodes_alloc; // Number of available nodes
} PrfTree;


PrfTree * prf_tree_new();

// Called by prf_tree_new
void prf_tree_set_defaults(PrfTree *T);

void prf_tree_show(FILE *, PrfTree *);
void prf_tree_init_varlist(PrfTree *);

void prf_tree_free(PrfTree * );
/* Create a tree, based on a table with
   N rows and M+1 columns
   The last column are the classes */
void prf_tree_train(PrfTree *,
                         const float * X);

PrfNode * prf_tree_train_with_parent(PrfTree *,
                         PrfNode * parent,
                                     const size_t N, // Number of features
                                     const float * X,
    float gini);


/* Classify one vector, X, by the tree T */
float prf_tree_classify(const PrfTree * T, const float * X);

/* Count all nodes/ set id */
size_t prf_tree_enumerate(PrfTree * T);
/* Free a node and all children recursively */
void prf_node_free(PrfNode * N);


/* A more compact form 196-bit */
// TODO

typedef struct{
    uint32_t left;
    uint32_t right;
    uint32_t var;
    uint32_t class;
    float th;
} PrfNodeRow;

typedef struct{
    size_t Nrows;
    PrfNodeRow * TT;
} PrfTreeTable;

PrfTreeTable * prf_tree_to_tree_table(PrfTree);
void prf_tree_table_free(PrfTreeTable *);
int prf_tree_table_classify(PrfTreeTable *, float * X);
// Read a tree-table from disk
PrfTreeTable * prf_tree_table_from_file(FILE *);
void prf_tree_table_to_file(PrfTreeTable *, FILE *);

#endif
