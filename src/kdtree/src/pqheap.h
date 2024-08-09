#pragma once

/*
 * A max-heap with constant size that can be used as a priority queue
 * keeping the k smallest elements.
*/


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>

// Items in the heap
typedef struct{
    double value;
    uint64_t idx;
} pqheap_item_t;

// Actual heap
struct pqheap {
    pqheap_item_t * G; // graph
    size_t n; // n items in graph/heap
    size_t k; // capacity, max items in G
};

typedef struct pqheap pqheap_t;


// Allocate a binary heap with
// a fixed max capacity, k
pqheap_t * pqheap_new(int k);
void pqheap_free(pqheap_t ** Bp);
// Insert into the heap
void pqheap_insert(pqheap_t * restrict B, double val, uint64_t id);
// Get pop out the smallest item
int pqheap_pop(pqheap_t * restrict B, double * val, uint64_t * id);
void pqheap_cascade_down(pqheap_t * B);
double pqheap_get_max_value(pqheap_t *);

// For testing/debugging:
void pqheap_validate(pqheap_t * B);
void pqheap_show(pqheap_t * B);
