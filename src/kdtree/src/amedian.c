#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <gsl/gsl_statistics_double.h>

static double timespec_diff(struct timespec* end, struct timespec * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}

double amedian(double * A, size_t N, double Amin, double Amax)
{
    // Histogram-based quick approximate median when the bounds are known
    const int Hn = 256;
    size_t * H = calloc(Hn, sizeof(size_t));
    for(size_t kk = 0; kk<N; kk++)
    {
        size_t idx = (A[kk]-Amin)/Amax*Hn;
        H[idx]++;
    }
    // Integrate to find the correct bin
    for(int kk = 1; kk<Hn; kk++)
    {
        H[kk] += H[kk-1];
    }

    for(int kk = 0; kk<Hn; kk++)
    {
        //printf("H[%d] = %zu\n", kk, H[kk]);
        if(H[kk] >= N/2.0)
        {
            return (double) kk/(double) Hn*(Amax-Amin)+Amin;
        }
    }
    assert(0);
    return -1;
}


int main(int argc, char ** argv)
{
    size_t N = 1500000;
    if(argc > 1)
    {
        N = atol(argv[1]);
    }
    double * X = malloc(N*sizeof(double));

    double Xmin = 0;
    double Xmax = 1;
    for(size_t kk = 0; kk<N; kk++)
    {
        X[kk] = (double) rand() / (double) RAND_MAX;
        X[kk] *= (Xmax-Xmin);
        X[kk] += Xmin;
        //X[kk] = kk;
    }
    //Xmin = 0;
    //Xmax = N-1;

    struct timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);
    double median0 = amedian(X, N, Xmin, Xmax);
    clock_gettime(CLOCK_REALTIME, &tend);
    double t0 = timespec_diff(&tend, &tstart);

    clock_gettime(CLOCK_REALTIME, &tstart);
    // note: changes the input array
    double median1 =  gsl_stats_median(X, 1, N);
    clock_gettime(CLOCK_REALTIME, &tend);
    double t1 = timespec_diff(&tend, &tstart);

    printf("gsl_stats_median = %f, t = %f s\n", median1, t1);
    printf("         amedian = %f, t = %f s\n", median0, t0);
    free(X);
    return EXIT_SUCCESS;
}
