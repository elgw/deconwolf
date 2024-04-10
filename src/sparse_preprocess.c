#include "sparse_preprocess.h"

#define LL /* For log-likelihood (not L2 ) */
#define ff

typedef struct {
    const float * restrict image;
    float lambda;
    float lambda_s;
    u64 M;
    u64 N;
    u64 P;
    float * restrict buff;
    int periodic;
} sparse_conf_t;

/** Note: this function does not care about the smoothness prior at
 * the boundary. */
static float
get_error(const sparse_conf_t * s, const float * restrict u)
{
    const size_t n = s->M*s->N*s->P;
    const float * restrict I = s->image;
    const u64 M = s->M;
    const u64 N = s->N;
    const u64 P = s->P;

    double E_data = 0;
#pragma omp parallel for reduction(+:E_data)
    for(size_t kk = 0; kk<n; kk++)
    {
        if(u[kk] > 1 && I[kk] > 1)
        {
#ifdef LL
#ifdef ff
            E_data += I[kk]*logf((I[kk]+1)/(u[kk]+1)) + u[kk] - I[kk];
#else
            E_data += u[kk]*logf(u[kk]/I[kk]) + I[kk] - u[kk];
#endif

#else
            E_data += powf(I[kk]-u[kk], 2)/2;
#endif
        }
    }

    double E_sparse = 0;
#pragma omp parallel for reduction(+:E_sparse)
    for(size_t kk = 0; kk<n; kk++)
    {
        assert(u[kk] >= 0);
        E_sparse += u[kk];
    }

    double E_smooth = 0;

    /* x */
#pragma omp parallel for reduction(+:E_smooth)
    for(size_t pp = 0; pp < P; pp++)
    {
        for(size_t nn = 0; nn < N; nn++)
        {
            const float * restrict lu = u + pp*M*N + nn*M;
            for(size_t mm = 1; mm+1 < M ; mm++)
            {
                E_smooth += fabs(lu[mm-1] + lu[mm+1] - 2*lu[mm]);
            }
        }
    }

    /* y */
#pragma omp parallel for reduction(+:E_smooth)
    for(size_t pp = 0; pp < P; pp++)
    {
        for(size_t nn = 1; nn+1 < N; nn++)
        {
            const float * restrict lu = u + pp*M*N + nn*M;
            for(size_t mm = 0; mm < M ; mm++)
            {
                E_smooth += fabs(lu[mm-M] + lu[mm+M] - 2*lu[mm]);
            }
        }
    }

    /* z  */
    {
        const size_t MN = M*N;
#pragma omp parallel for reduction(+:E_smooth)
        for(size_t pp = 1; pp < P-1; pp++)
        {
            for(size_t nn = 0; nn < N; nn++)
            {
                const float * restrict lu = u + pp*M*N + nn*M;
                for(size_t mm = 0; mm < M ; mm++)
                {
                    E_smooth += fabs(lu[mm-MN] + lu[mm+MN] - 2*lu[mm]);
                }
            }
        }
    }

    //printf("E_data=%e, E_smooth=%e, E_sparse=%e\n", E_data, E_smooth, E_sparse);

    return (s->lambda*E_data + E_smooth + s->lambda_s*E_sparse) / (double) (M*N*P);
}

static float max(const float a, const float b)
{
    if(a > b)
        return a;
    return  b;

}

static float asign(const float x)
{
    return (x >= 0) - (x < 0);
}

static void get_gradient(const sparse_conf_t * restrict s,
                         const float * restrict u,
                         float * restrict dE)
{

    const float * restrict I = s->image;
    const u64 M = s->M;
    const u64 N = s->N;
    const u64 P = s->P;
    const size_t n = M*N*P;

    /* Data fidelity and sparsity */
#pragma omp parallel for
    for(size_t kk = 0 ; kk < n; kk++)
    {
#ifdef LL
#ifdef ff
        /* Prevent the gradient to explode with the max ...  */

        dE[kk] = s->lambda * (1.0 - (I[kk]+1.0)/(u[kk] + 1.0)) + s->lambda_s;

        //dE[kk] = s->lambda * (1.0-I[kk]/max(u[kk], 0.5)) + s->lambda_s;

#else
        dE[kk] = s->lambda * logf(max(u[kk], 1e-4)/max(I[kk], 1e-4)) + s->lambda_s;
#endif
#else
        dE[kk] = s->lambda*(u[kk]-I[kk]) + s->lambda_s;
#endif
    }


    /* x */
#pragma omp parallel for
    for(size_t pp = 0; pp < P; pp++)
    {
        for(size_t nn = 0; nn < N; nn++)
        {
            const float * restrict lu = u + pp*M*N + nn*M;
            float * restrict lbuff = s->buff + pp*M*N + nn*M;
            if(s->periodic)
            {
                lbuff[0] = asign(lu[M-1] + lu[1] - 2*lu[0]);
            } else {
                lbuff[0] = asign(lu[1] - lu[0]);
            }
            for(size_t mm = 1; mm+1 < M ; mm++)
            {
                //lbuff[mm] = copysignf(1, lI[mm-1] + lI[mm+1] - 2*lI[mm]);
                lbuff[mm] = asign(lu[mm-1] + lu[mm+1] - 2*lu[mm]);
            }
            if(s->periodic)
            {
                lbuff[M-1] = asign(lu[M-2] + lu[0] - 2*lu[M-1]);
            } else {
                lbuff[M-1] = asign(lu[M-2]-lu[M-1]);
            }

        }
    }

#pragma omp parallel for
    for(size_t pp = 0; pp < P; pp++)
    {
        for(size_t nn = 0; nn < N; nn++)
        {
            float * restrict lbuff = s->buff + pp*M*N + nn*M;
            float * restrict ldE = dE + pp*M*N + nn*M;;

            if(s->periodic)
            {
                ldE[0] += -2*lbuff[0] + lbuff[1] + lbuff[M-1];
            } else {
                ldE[0] = -lbuff[0] + lbuff[1];
            }

            for(size_t mm = 1; mm+1 < M ; mm++)
            {
                ldE[mm] += (-2*lbuff[mm] + lbuff[mm-1] + lbuff[mm+1]);
            }
            if(s->periodic)
            {
                ldE[M-1] = lbuff[0] - 2*lbuff[M-1] + lbuff[M-2];
            } else {
                ldE[M-1] = -lbuff[M-1]+lbuff[M-2];
            }
        }
    }


    /* y */
#pragma omp parallel for
    for(size_t pp = 0; pp < P; pp++)
    {
        for(size_t nn = 0; nn == 0; nn++)
        {
            const float * restrict lu = u + pp*M*N +nn*M;
            float * restrict lbuff = s->buff + pp*M*N+nn*M;

            if(s->periodic)
            {
                for(size_t mm = 0; mm < M ; mm++)
                {
                    lbuff[mm] = asign(lu[mm+M] - 2*lu[mm] + u[ pp*M*N + M*(N-1) + mm ]);
                }
            } else {
                for(size_t mm = 0; mm < M ; mm++)
                {
                    lbuff[mm] = asign(lu[mm+M] - lu[mm]);
                }
            }
        }

        for(size_t nn = 1; nn+1 < N; nn++)
        {
            const float * restrict lu = u + pp*M*N+nn*M;
            float * restrict lbuff = s->buff + pp*M*N+nn*M;

            for(size_t mm = 0; mm < M ; mm++)
            {
                lbuff[mm] = asign(lu[mm-M] + lu[mm+M] - 2*lu[mm]);
            }
        }

        for(size_t nn = N-1; nn == N-1; nn++)
        {
            const float * restrict lu = u + pp*M*N+nn*M;
            float * restrict lbuff = s->buff + pp*M*N+nn*M;

            if(s->periodic)
            {
                for(size_t mm = 0; mm < M ; mm++)
                {
                    lbuff[mm] = asign(lu[mm-M] - 2*lu[mm] + u[ pp*M*N + mm ]);
                }
            } else {
                for(size_t mm = 0; mm < M ; mm++)
                {
                    lbuff[mm] = asign(lu[mm-M] - lu[mm]);
                }
            }
        }
    }

#pragma omp parallel for
    for(size_t pp = 0; pp < P; pp++)
    {
        for(size_t nn = 0; nn== 0 ; nn++)
        {
            float * restrict lbuff = s->buff + pp*M*N+nn*M;
            float * restrict ldE = dE + pp*M*N+nn*M;
            if(s->periodic){
                for(size_t mm = 0; mm < M ; mm++)
                {
                    ldE[mm] += (-2*lbuff[mm] + lbuff[mm+M] + s->buff[pp*M*N+mm]);
                }
            } else {
                for(size_t mm = 0; mm < M ; mm++)
                {
                    ldE[mm] += (-lbuff[mm] + lbuff[mm+M]);
                }
            }
        }

        for(size_t nn = 1; nn+1 < N; nn++)
        {
            float * restrict lbuff = s->buff + pp*M*N+nn*M;

            float * restrict ldE = dE + pp*M*N+nn*M;;

            for(size_t mm = 0; mm < M ; mm++)
            {
                ldE[mm] += (-2*lbuff[mm] + lbuff[mm-M] + lbuff[mm+M]);
            }
        }

        for(size_t nn = N-1; nn== N-1 ; nn++)
        {
            float * restrict lbuff = s->buff + pp*M*N+nn*M;
            float * restrict ldE = dE + pp*M*N+nn*M;;
            if(s->periodic){
                for(size_t mm = 0; mm < M ; mm++)
                {
                    ldE[mm] += (-2*lbuff[mm] + lbuff[mm-M] + s->buff[pp*M*N+M*(N-1) + mm]);
                }
            } else {

                for(size_t mm = 0; mm < M ; mm++)
                {
                    ldE[mm] += (-lbuff[mm] + lbuff[mm-M]);
                }
            }
        }

    }

    /* z */
    const size_t MN = s->M * s->N;


#pragma omp parallel for
    for(size_t pp = 0; pp < P; pp++)
    {
        if(pp == 0)
        {
            if(s->periodic) {
                for(size_t nn = 0; nn < N; nn++)
                {
                    const float * restrict lu = u + pp*M*N + nn*M;
                    float * restrict lbuff = s->buff + pp*M*N + nn*M;

                    for(size_t mm = 0; mm < M ; mm++)
                    {
                        lbuff[mm] = asign(lu[mm+MN] - 2*lu[mm]+ u[nn*M + mm + M*N*(P-1)]);
                    }
                }
            } else {
                for(size_t nn = 0; nn < N; nn++)
                {
                    const float * restrict lu = u + pp*M*N + nn*M;
                    float * restrict lbuff = s->buff + pp*M*N + nn*M;

                    for(size_t mm = 0; mm < M ; mm++)
                    {
                        lbuff[mm] = asign(lu[mm+MN] - lu[mm]);
                    }
                }
            }
        } else if(pp == P-1)
        {
            if(s->periodic) {
                for(size_t nn = 0; nn < N; nn++)
                {
                    const float * restrict lu = u + pp*M*N + nn*M;
                    float * restrict lbuff = s->buff + pp*M*N + nn*M;

                    for(size_t mm = 0; mm < M ; mm++)
                    {
                        lbuff[mm] = asign(lu[mm-MN] - 2*lu[mm] + u[nn*M + mm]);
                    }
                }
            } else {
                for(size_t nn = 0; nn < N; nn++)
                {
                    const float * restrict lu = u + pp*M*N + nn*M;
                    float * restrict lbuff = s->buff + pp*M*N + nn*M;

                    for(size_t mm = 0; mm < M ; mm++)
                    {
                        lbuff[mm] = asign(lu[mm-MN] - lu[mm]);
                    }
                }
            }
        } else {

            for(size_t nn = 0; nn < N; nn++)
            {
                const float * restrict lu = u + pp*M*N + nn*M;
                float * restrict lbuff = s->buff + pp*M*N + nn*M;

                for(size_t mm = 0; mm < M ; mm++)
                {
                    lbuff[mm] = asign(lu[mm-MN] + lu[mm+MN] - 2*lu[mm]);
                }
            }
        }
    }

#pragma omp parallel for
    for(size_t pp = 0; pp < P; pp++)
    {
        if(pp == 0)
        {
            for(size_t nn = 0; nn < N; nn++)
            {
                float * restrict lbuff = s->buff + pp*M*N + nn*M;
                float * restrict ldE = dE + pp*M*N + nn*M;
                if(s->periodic == 1)
                {
                    for(size_t mm = 0; mm < M ; mm++)
                    {
                        ldE[mm] += (-2*lbuff[mm] + lbuff[mm+MN] + s->buff[(P-1)*M*N + nn*M + mm]);
                    }
                } else {
                    for(size_t mm = 0; mm < M ; mm++)
                    {
                        ldE[mm] += (-lbuff[mm] + lbuff[mm+MN]);
                    }
                }
            }
        } else if(pp == P-1)
        {
            for(size_t nn = 0; nn < N; nn++)
            {
                float * restrict lbuff = s->buff + pp*M*N + nn*M;
                float * restrict ldE = dE + pp*M*N + nn*M;
                if(s->periodic == 1)
                {
                    for(size_t mm = 0; mm < M ; mm++)
                    {
                        ldE[mm] += (-2*lbuff[mm] + lbuff[mm-MN] + s->buff[nn*N+mm]);
                    }
                } else {
                    for(size_t mm = 0; mm < M ; mm++)
                    {
                        ldE[mm] += (-lbuff[mm] + lbuff[mm-MN]);
                    }
                }
            }
        } else {
            for(size_t nn = 0; nn < N; nn++)
            {
                float * restrict lbuff = s->buff + pp*M*N + nn*M;
                float * restrict ldE = dE + pp*M*N + nn*M;

                for(size_t mm = 0; mm < M ; mm++)
                {
                    ldE[mm] += (-2*lbuff[mm] + lbuff[mm-MN] + lbuff[mm+MN]);
                }
            }
        }
    }

    return;
}

/* dtu = u + step*dE */
void
add(float * restrict dtu,
    const float * restrict u,
    const float dt,
    const float * restrict dE,
    const size_t n)
{
#pragma omp parallel for
    for(size_t kk = 0; kk < n ; kk++)
    {
        dtu[kk] = max(u[kk] + dt*dE[kk], 0);
    }
    return;
}


float *
sparse_preprocess(const float * image, const u64 M, const u64 N, const u64 P,
                  const double lambda, const double lambda_s,
                  const int periodic, const u64 iter)
{

#ifdef LL
    printf("Using Log-Likelihood (I-divergence) for the data fidelity\n");
#else
    printf("Using L2-norm for data fidelity\n");
#endif

    printf("lambda = %e, lambda_s = %e\n", lambda, lambda_s);

    sparse_conf_t * s = calloc(1, sizeof(sparse_conf_t));
    assert(s != NULL);
    s->image = image;
    s->lambda = lambda;
    s->lambda_s = lambda_s;
    s->M = M;
    s->N = N;
    s->P = P;
    s->buff = calloc(M*N*P, sizeof(float));
    s->periodic = periodic;


    const size_t n = M*N*P;
    const float min_step = 1e-6;
    float dt = 1;

    float * u = calloc(M*N*P, sizeof(float));
    assert(u != NULL);
    memcpy(u, image, M*N*P*sizeof(float));

    float * dtu = calloc(M*N*P, sizeof(float));
    assert(dtu != NULL);

    float * dE = calloc(M*N*P, sizeof(float));
    assert(dE != NULL);

    float E0 = get_error(s, u);
    float E = 0;

    printf("processing ... \n");
    for(u64 ii = 0; ii<iter; ii++)
    {
        printf("\r%03lu/%03lu: E=%e, step=%f", ii+1, iter, E0, dt);
        fflush(stdout);
        get_gradient(s, u, dE);


        /* dtu = u - dt*dE */

        add(dtu, u, -dt, dE, n);
        E = get_error(s, dtu);
        //        printf(" ---  E=%e\n", E);

        /* If the step was too long and didn't decrease the error,
         * decrease the step size and try again */
        // dt = 1; /* Restart. Does not help */
        while((E >= E0) && (dt > min_step))
        {
            //            printf("step: %f -> %f\n", step, 0.6*step);
            dt *= 0.8;
            /* dtu = u + step*dE */
            add(dtu, u, -dt, dE, n);
            E = get_error(s, dtu);
        }

        { /* Swap to avoid a copy */
            float * t = u;
            u = dtu;
            dtu = t;
        }

        if(dt < min_step)
        {
            printf("\n");
            printf("Stopping because dt < min_step (%e < %e)\n", dt, min_step);
            break;
        }
        E0 = E;
        //dt*=2;
        //dt > 1 ? dt = 1 : 0 ;
    }
    printf("\n");


    free(dE);

    free(dtu);
    free(s->buff);
    free(s);

    return u;
}
