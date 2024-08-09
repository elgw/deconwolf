#include "pqheap.h"

static void pqheap_heapify(pqheap_t * restrict B);

pqheap_t * pqheap_new(int k)
{
    pqheap_t * B = calloc(1, sizeof(pqheap_t));
    assert(B != NULL);
    B->n = 0;
    B->k = k;
    B->G = calloc(k, sizeof(pqheap_item_t));
    assert(B->G != NULL);
    return B;
}

static void pqheap_heapify(pqheap_t * restrict B)
{
    // Something is wrong here, should only need to
    // compare to parents and swap with them ...


#if 1
    pqheap_item_t * G = B->G;
    size_t n = B->n-1;
    size_t p = (n-1)/2;
    while(p != (size_t) -1)
    {
        if(G[p].value < G[n].value)
        {
            pqheap_item_t t = G[p];
            G[p] = G[n];
            G[n] = t;
            n = p;
            p = (n-1)/2;
        } else {
            p = -1;
        }
    }
    //pqheap_validate(B);
    return;
#else
    // Complicates the procedure: not all parents have two children
    //    printf("Before heapify:\n");
    //    pqheap_show(B);
    pqheap_item_t * G = B->G;

    // make sure that the heap is a heap
    // assumes that only one element was inserted since last call
    size_t n = B->n-1;

    size_t p = (n-1)/2;
    size_t c1 = p*2+1;
    size_t c2 = p*2+2;
    //    printf("n=%zu, p = %zu, c1 = %zu, c2 = %zu\n", n, p, c1, c2);
    if(c2 > n)
    {
        c2=c1;
    }
    if(G[c2].value > G[c1].value)
    {
        c1 = c2;
    }
    double cmax = G[c1].value;
    double pvalue = G[p].value;
    //    printf("pvalue[%zu]: %f, cmax[%zu]: %f\n", p, pvalue, c1, cmax);

    while ( pvalue < cmax )
    {
        //        printf("1 pvalue[%zu]: %f, cmax[%zu]: %f\n", p, pvalue, c1, cmax);
        // swap parent with the largest
        pqheap_item_t t = B->G[p];
        B->G[p] = B->G[c1];
        B->G[c1] = t;

        if(p == 0)
        {
            break;
        }
        p = (p-1)/2;
        c1 = p*2+1;
        c2 = p*2+2;

        if(G[c2].value > G[c1].value)
        {
            c1 = c2;
        }
        cmax = G[c1].value;
        pvalue = G[p].value;
        //        printf("2 pvalue[%zu]: %f, cmax[%zu]: %f\n", p, pvalue, c1, cmax);
 }
    //  printf("After heapify\n");
    //pqheap_show(B);
    pqheap_validate(B);
    #endif
}

void pqheap_insert(pqheap_t * restrict B, double val, uint64_t id)
{

    if( B->n == B->k) // if full
    {
        // Either reject if too large
        if(val > B->G[0].value)
        {
            return; // silently reject the offer
        }
        // Or replace the largest with the new value
        B->G[0].value = val;
        B->G[0].idx = id;
        pqheap_cascade_down(B);
        return;
   }

    // Normal insertion

    size_t n = B->n++;
    //    printf("Inserting (%f, %lu) at %zu\n", val, id, n);
    B->G[n].value = val;
    B->G[n].idx = id;
    if(n>0)
    {
        pqheap_heapify(B);
    }
}

void pqheap_show(pqheap_t * B)
{
    for(size_t kk = 0; kk<B->n; kk++)
    {
        printf("G[%zu] = (%f, %lu)\n", kk, B->G[kk].value, B->G[kk].idx);
    }
}

void pqheap_validate(pqheap_t * B)
{
    for(size_t kk = 1; kk<B->n; kk++)
    {
        size_t n = kk;
        size_t parent = (n-1)/2;
        int fail = 0;
        if( !( B->G[n].value < B->G[parent].value ))
        {
            printf("G[%zu].value = %f !< G[%zu].value = %f\n",
                   n, B->G[n].value,
                   parent, B->G[parent].value);
            fail = 1;
        };
        if(fail)
        {
            printf("pqheap_validate failed\n");
            pqheap_show(B);
            assert(0);
        }
    }
}

void pqheap_cascade_down(pqheap_t * B)
{
    pqheap_item_t * G = B->G;
    // cascade down the non-heaped branch
    size_t n = 0;
    while( 2*n+1 < B->n ) // while there are children
    {
        size_t c = 2*n+1;
        size_t c2 = 2*n+2;
        //printf("%d -> (%d, %d)\n", n, c, c2);
        if(c2 < B->n)
        {
            if(B->G[c2].value > B->G[c].value)
            {
                c = c2;
            }
        }
        //printf("n=%d, c=%d\n", n, c);
        if(G[n].value > G[c].value)
        {
            return;
        }
        // swap element n and c
        pqheap_item_t t = G[c];
        G[c] = G[n];
        G[n] = t;
        n = c;
    }
    return;
}

int pqheap_pop(pqheap_t * restrict B, double * v, uint64_t * i)
{
    pqheap_item_t * restrict G = B->G;

    // Store value to be returned
    v[0] = G[0].value;
    i[0] = G[0].idx;


    if(B->n == 1)
    { // If we popped the last element
        B->n = 0;
        return 0;
    }

    // move last element first
    G[0] = G[(B->n--) -1];

    pqheap_cascade_down(B);

    return 0;
}

void pqheap_free(pqheap_t ** Bp)
{
    pqheap_t * B = Bp[0];
    free(B->G);
    free(B);
    Bp[0] = NULL;
}

double pqheap_get_max_value(pqheap_t * pq)
{
    return pq->G[0].value;
}
