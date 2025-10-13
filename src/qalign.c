#include "qalign.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

typedef int64_t i64;
typedef int32_t i32;


qalign_config * qalign_config_new(void)
{
    qalign_config * qconf = malloc(sizeof(qalign_config));
    qconf->similarity_threshold = 0.5;
    qconf->deviation_x = 1;
    qconf->verbose = 1;
    qconf->npoint = 1000;
    return qconf;
}

void qalign_config_free(qalign_config * q)
{
    if(q == NULL)
    {
        return;
    }
    free(q->H);
    free(q);
    return;
}

static int qalign_cmpx( const void * _A, const void * _B)
{
    qalign_pair * A = (qalign_pair*) _A;
    qalign_pair * B = (qalign_pair*) _B;
    if(A->dx > B->dx)
    {
        return 1;
    }
    if(A->dx < B->dx)
    {
        return -1;
    }
    return 0;
}

static void sort_pairs(qalign_pair * P, i64 nP)
{
    /* Sort by dx, smallest first */
    qsort(P, nP, sizeof(qalign_pair), qalign_cmpx);
}

static void get_pairs(const float * X, i64 nX, qalign_pair * P)
{
    /* Construct all $(nX-1)nX/2$ difference vectors from the points in X */

    i64 idx = 0;
    for(i64 kk = 0; kk < nX; kk++)
    {
        for(i64 ll = kk+1; ll < nX; ll++)
        {
            for(i64 ii = 0; ii < 3; ii++)
            {
                P[idx].delta[ii] = X[3*kk + ii] - X[3*ll + ii];
            }
            P[idx].idx1 = kk;
            P[idx].idx2 = ll;
            idx++;
        }
    }
}

static float float3_norm(float * v)
{
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static float qpair_distance(const qalign_pair * A,
                            const qalign_pair * B)
{
    float delta[3] = {
        A->dx - B->dx,
        A->dy - B->dy,
        A->dz - B->dz};

    return float3_norm(delta);
}

static float qpair_distance_swap(const qalign_pair * A,
                                 const qalign_pair * B)
{
    float delta[3] = {
        A->dx + B->dx,
        A->dy + B->dy,
        A->dz + B->dz};

    return float3_norm(delta);
}


static void find_hypotheses(qalign_config * config,
                            const qalign_pair * pA,
                            i64 npA,
                            const qalign_pair * pB,
                            i64 npB,
                            qalign_pair * H, i64 nH_alloc, i64 * nH)
{
    /* Loop over all pairs in pA and see if we can find correponding
       pairs in pB. If found they are added to array of Hypotheses, H
    */
    i64 iB0 = 0; // range in B: [iB0, iB1]
    float deviation_x = config->deviation_x;
    float similarity_threshold = config->similarity_threshold;

    for(i64 iA = 0; iA < npA; iA++)
    {
        qalign_pair A = pA[iA];
        //printf("iA=%ld, %f, %f, %f\n", iA, A.dx, A.dy, A.dz);
        for(i64 iB = iB0; iB < npB; iB++)
        {
            qalign_pair B = pB[iB];
            //printf("iA=%ld, %f, %f, %f ", iA, A.dx, A.dy, A.dz);
            //printf("iB=%ld, %f, %f, %f\n", iB, B.dx, B.dy, B.dz);

            /* Manage the start position in pB */
            if(B.dx + deviation_x < A.dx )
            {
                iB0 = iB+1;
            }
            /* Don't compare more than necessary */
            if(B.dx > A.dx + deviation_x)
            {
                break;
            }


            // TODO:
            // Model points with measurement error e, i.e., something like
            // X = X' + e, i.e. include
            // where X' is the true location and X the measured

            // Handle Direction ambivalence:
            // -- test with both A vs B and A vs -B
            i64 idxA = A.idx1;
            i64 idxB = B.idx1;
            float diff = qpair_distance(&A, &B);
            float diff2 = qpair_distance_swap(&A, &B);
            if(diff2 < diff)
            {
                diff = diff2;
                idxB = B.idx2;
            }

            float d1 = float3_norm(A.delta);
            float d2 = float3_norm(B.delta);
            float dm = d1;
            dm < d2 ? dm = d2 : 0;

            if(dm > 1.0) // very short distance vectors ignored
            {
                if(diff / dm < 1e-4) // relative error
                {

                    H[*nH].dx = config->A[3*idxA + 0] - config->B[3*idxB+0];
                    H[*nH].dy = config->A[3*idxA + 1] - config->B[3*idxB+1];
                    H[*nH].dz = config->A[3*idxA + 2] - config->B[3*idxB+2];
                    H[*nH].idx1 = iA;
                    H[*nH].idx2 = iB;
                    (*nH)++;
                    if(*nH == nH_alloc)
                    {
                        return;
                    }
                }
            }

        }
    }
    return;
}

int qalign(qalign_config * conf)
{
    assert(conf != NULL);
    assert(conf->A != NULL);
    assert(conf->B != NULL);
    if(conf->nA < 2)
    {
        if(conf->verbose > 0)
        {
            printf("qalign: Too few dots in set A\n");
        }
        return EXIT_FAILURE;
    }
    if(conf->nB < 2)
    {
        if(conf->verbose > 0)
        {
            printf("qalign: Too few dots in set B\n");
        }
        return EXIT_FAILURE;
    }
    if(conf->nA > conf->npoint)
    {
        conf->nA = conf->npoint;
    }
    if(conf->nB > conf->npoint)
    {
        conf->nB = conf->npoint;
    }

    i64 npA = conf->nA * (conf->nA-1) / 2;
    qalign_pair * pA = malloc(npA*sizeof(qalign_pair));
    if(pA == NULL)
    {
        return EXIT_FAILURE;
    }

    i64 npB = conf->nB * (conf->nB-1) / 2;
    qalign_pair * pB = malloc(npB*sizeof(qalign_pair));
    if(pB == NULL)
    {
        return EXIT_FAILURE;
    }

    if(conf->verbose > 1)
    {
        printf("Finding pairs from set A\n");
    }
    get_pairs(conf->A, conf->nA, pA);
    if(conf->verbose > 1)
    {
        printf("Finding pairs from set B\n");
    }
    get_pairs(conf->B, conf->nB, pB);

    if(conf->verbose > 1)
    {
        printf("sorting pairs\n");
    }
    sort_pairs(pA, npA);
    sort_pairs(pB, npB);

    i64 nH_alloc = 1000;
    i64 nH = 0;
    qalign_pair * H = malloc(nH_alloc*sizeof(qalign_pair));
    conf->H = H;
    if(H == NULL)
    {
        goto cleanup;
    }

    if(conf->verbose > 1)
    {
        printf("Building hypotheses\n");
    }
    find_hypotheses(conf,
                    pA, npA,
                    pB, npB,
                    H, nH_alloc, &nH);
    conf->nH = nH;
    if(conf->verbose > 1)
    {
        printf("Found %ld hypotheses\n", nH);


        int nshow = nH;
        nshow > 10 ? nshow = 10 : 0;
        for(int kk = 0; kk < nshow; kk++)
        {
            qalign_pair * h = H + kk;
            printf("[%f, %f, %f]\n", h->dx, h->dy, h->dz);
        }
    }

    free(pA);
    free(pB);
    return EXIT_SUCCESS;

cleanup:
    free(pA);
    free(pB);
    return EXIT_FAILURE;
}
