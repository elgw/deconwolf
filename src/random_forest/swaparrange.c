/* Conclusion building a list of values to swap is slower than
making a copy of the data and then writing back. */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <limits.h>
#include <string.h>

int cmp_uint32(const void * a, const void * b)
{
    uint32_t * A = (uint32_t *) a;
    uint32_t * B = (uint32_t *) b;
    return A[0]>B[0];
}

void randperm_uint32(uint32_t * X, const size_t N)
{
    uint32_t * X2 = malloc(2*N*sizeof(uint32_t));
    for(int kk = 0; kk<N; kk++)
    {
        X2[2*kk] = rand() % INT_MAX; // Not truly random even if rand() would be
        X2[2*kk+1] = X[kk];
    }

    qsort(X2, N, 2*sizeof(uint32_t), cmp_uint32);

    for(int kk = 0; kk<N; kk++)
    {
        X[kk] = X2[2*kk+1];
    }
    free(X2);
}

size_t rearrange(float * X, size_t N, float * L, float * R, float th)
{
    memcpy(L, X, N*sizeof(float));
    size_t lpos = 0;
    size_t rpos = 0;
    for(size_t kk = 0; kk<N; kk++)
    {
        if(X[kk] , th)
        {
            L[lpos++] = X[kk];
        } else {
            R[rpos++] = X[kk];
        }
    }
    return lpos;
}

size_t gen_swaplists(float * X, size_t N, uint32_t * idx, uint8_t * B, size_t nL,
                     size_t * S1, size_t * S2)
{
// X[0] .. X[nL-1] should contain only indexes from idx[0] ... idx[nL-1]

    // 1. Create binary list that says where elements should go
    for(size_t kk = 0; kk<nL; kk++)
    {
        B[idx[kk]] = 0;
    }
    for(size_t kk = nL; kk<N; kk++)
    {
        B[idx[kk]] = 1;
    }

    size_t S1pos = 0;
    for(size_t kk = 0; kk<nL; kk++)
    {
        if(B[kk] == 1)
        {
            S1[S1pos++] = kk;
        }
    }

    size_t S2pos = 0;
    for(size_t kk = nL; kk<N; kk++)
    {
        if(B[kk] == 0)
        {
            S2[S2pos++] = kk;
        }
    }

    assert(S1pos == S2pos);
    return S1pos;
}

size_t use_swaplists(float * X, size_t * S1, size_t * S2, size_t nswap)
{
    for(size_t kk = 0; kk<nswap; kk++)
    {
        float t = X[S1[kk]];
        X[S1[kk]] = X[S2[kk]];
        X[S2[kk]] = t;
    }
    return (size_t) X[1];
}

int main(int argc, char ** argv)
{
    size_t N = 10000;
    size_t niter = 1;
    size_t nvar = 2;

    int method = 0;

    if(argc > 1)
    {
        method = atoi(argv[1]);
    }

    if(argc > 2)
    {
        N = atol(argv[2]);
    }
    if(argc > 3)
    {
        niter = atol(argv[3]);
    }

    if(argc > 4)
    {
        nvar = atol(argv[4]);
    }

    if(method == 0)
    {
        printf("Method: Movie all data according to th\n");
    }

    if(method == 1)
    {
        printf("Build and use swaplists\n");
    }

    printf("N = %zu, niter = %zu nvariables: %zu\n", N, niter, nvar);

    float * X = malloc(N*sizeof(float));
    uint32_t * I = malloc(N*sizeof(uint32_t));
    size_t nleft = N/2;

    for(size_t kk = 0; kk<N; kk++)
    {
        X[kk] = (float) rand() / (float) RAND_MAX;
        I[kk] = kk;
    }
    for(int kk = 0; kk<10; kk++)
    {
        printf(" %u ", I[kk]);
    }
    printf("\n");

    randperm_uint32(I, N);

    for(int kk = 0; kk<10; kk++)
    {
        printf(" %u ", I[kk]);
    }
    printf("\n");

    float * L = malloc(N*sizeof(float));
    float * R = malloc(N*sizeof(float));
    uint8_t * B = malloc(N*sizeof(float));
    size_t * S1 = malloc(N*sizeof(size_t));
    size_t * S2 = malloc(N*sizeof(size_t));

    float sum = 0;
    for(size_t kk = 0; kk<N; kk++)
    {
        sum += X[kk];
    }

    float th = sum/ (float) N;

    printf("th: %f\n", th);

    size_t nswap = 0;

    size_t nlpos = 0;
    for(size_t ii = 0; ii< niter; ii++)
    {
        if(method == 0)
        {
            for(int vv = 0; vv<nvar; vv++)
            {
                nlpos += rearrange(X, N, L, R, th);
            }
        }
        if(method == 1)
        {
            nswap = gen_swaplists(X, N, I, B, N/2, S1, S2);
            for(int vv = 0; vv<nvar; vv++)
            {
                nlpos += use_swaplists(X, S1, S2, nswap);
            }
        }
    }
    printf("nlpos: %zu, nswap: %zu\n", nlpos, nswap);

    free(X);
    free(I);
    free(L);
    free(B);
    free(S1);
    free(S2);

    return 0;
}
