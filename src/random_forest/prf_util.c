#include "prf_util.h"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

void copy_row(float * A, const size_t sA,
              const float * B, const size_t sB,
              const size_t M)
{
    for(size_t mm = 0; mm < M; mm++)
    {
        A[mm*sA] = B[mm*sB];
    }
}

void get_subset_k(const float * X, const size_t N, const size_t M,
                  const int K, const int k,
                  // outputs:
                  float ** _Tr, size_t * _nTr,
                  float ** _Va, size_t * _nVa)
{
    assert(k >= 0);
    assert(k < K);

    ldiv_t d = ldiv(N-k, K);

    size_t nVa = d.quot + (d.rem > 0);
    size_t nTr = N - nVa;
    assert(nVa + nTr == N);
    //printf("K = %d, k = %d, N = %zu, nVa = %zu, nTr = %zu\n", K, k, N, nVa, nTr);
    float * Tr0 = malloc(nTr*M*sizeof(float));
    float * Tr = Tr0;
    float * Va0 = malloc(nVa*M*sizeof(float));
    float * Va = Va0;
    size_t t = 0; size_t v = 0;
    for(size_t nn = 0; nn < N; nn++)
    {
        if( nn % (size_t) K == (size_t )k )
        {
            // Write to validation set
            v++; assert(v <= nVa);
            copy_row(Va++, nVa, X+nn, N, M);

        } else {
            // Write to training set
            t++; assert(t <= nTr);
            copy_row(Tr++, nTr, X+nn, N, M);
        }
    }


    _nTr[0] = nTr;
    _nVa[0] = nVa;
    _Tr[0] = Tr0;
    _Va[0] = Va0;
}


void show_row_n(const float * T, size_t N, size_t M, size_t n)
{
    printf("%06zu ", n);
    for(size_t cc = 0; cc<M; cc++)
    {
        printf(" \t%.2f",T[n+cc*N]);
    }
    printf("\n");
}

size_t vector_size_t_argmax(const size_t * V, const size_t N)
{
    size_t max = 0;
    size_t argmax = 0;
    for(size_t kk = 0; kk<N; kk++)
    {
        if(V[kk] > max)
        {
            max = V[kk];
            argmax = kk;
        }
    }
    return argmax;
}


uint8_t * get_random_selection(size_t N, size_t n)
{
    uint8_t * use = malloc((N+1)*sizeof(uint8_t));
    use[N] = 1;

    if(n < N/2)
    {
        for(size_t kk = 0; kk<N; kk++)
        {
            use[kk] = 0;
        }

        size_t sel = 0;
        while(sel < n)
        {
            size_t pos = rand() % N;
            if(use[pos] == 0)
            {
                use[pos] = 1;
                sel++;
            }
        }
    }
    else
    {
        for(size_t kk = 0; kk<N; kk++)
        {
            use[kk] = 1;
        }

        size_t sel = N;
        while(sel != n)
        {
            size_t pos = rand() % N;
            if(use[pos] == 1)
            {
                use[pos] = 0;
                sel--;
            }
        }
    }

    return use;
}

void scramble_feature(float * X, const size_t N, const size_t M, const int f)
{
    assert(f>=0);
    assert(f<(int) M);
    float * XF = X + f*N;
    for(size_t kk = 0; kk<N; kk++)
    {
        XF[kk] = (float) rand() / (float) RAND_MAX;
    }
}


static void fill(int n, char c)
{
    for(int kk = 0; kk < n ; kk++)
        printf("%c", c);
}

static void line(int n)
{
    printf(" ");
    fill(n-1, '-');
    printf("\n");
}

void print_section(char * msg)
{

    int w = 80;
    line(w);

    int w1 = (w - 6 - strlen(msg))/2;
    int w2 = w - w1 - strlen(msg)-6;

    assert(3 + w1 + strlen(msg) + w2 + 3 == 80);
    printf(" --");
    fill(w1, ' ');
    printf("%s", msg);
    fill(w2, ' ');
    printf(" --\n");
    line(w);
}

#ifdef __APPLE__
size_t get_peakMemoryKB(void)
{
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);
  return (size_t) round((float) r_usage.ru_maxrss/1024.0);
}
#endif

#ifndef __APPLE__
size_t get_peakMemoryKB(void)
{
  char * statfile = malloc(100*sizeof(char));
  sprintf(statfile, "/proc/%d/status", getpid());
  FILE * sf = fopen(statfile, "r");
  if(sf == NULL)
  {
    fprintf(stderr, "Failed to open %s\n", statfile);
    free(statfile);
    return 0;
  }

  char * peakline = NULL;

  char * line = NULL;
  size_t len = 0;

  while( getline(&line, &len, sf) > 0)
  {
    if(strlen(line) > 6)
    {
      if(strncmp(line, "VmPeak", 6) == 0)
      {
        peakline = strdup(line);
      }
    }
  }
  free(line);
  fclose(sf);
  free(statfile);

  // Parse the line starting with "VmPeak"
  // Seems like it is always in kB
  // (reference: fs/proc/task_mmu.c)
  // actually in kiB i.e., 1024 bytes
  // since the last three characters are ' kb' we can skip them and parse in between
  size_t peakMemoryKB = 0;
  //  printf("peakline: '%s'\n", peakline);
  if(strlen(peakline) > 11)
  {
    peakline[strlen(peakline) -4] = '\0';

    //    printf("peakline: '%s'\n", peakline+7);
    peakMemoryKB = (size_t) atol(peakline+7);
  }

  free(peakline);
  return peakMemoryKB;
}
#endif


float timespec_diff(struct timespec* end, struct timespec * start)
{
    float elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}

void vector_show(float * v, int n)
{
    for(int kk = 0; kk<n; kk++)
        printf("%f ", v[kk]);
    printf("\n");
}

void print_todo(char * msg)
{
    printf("\e[4;33m" "TODO" ANSI_COLOR_RESET ": %s\n", msg);
}

void print_ok()
{
    printf(ANSI_COLOR_GREEN "ok\n" ANSI_COLOR_RESET);
}

void print_error(char * msg)
{
    printf(ANSI_COLOR_RED "ERROR" ANSI_COLOR_RESET ": %s\n", msg);
}

void print_warning(char * msg)
{
    printf(ANSI_COLOR_RED "WARNING" ANSI_COLOR_RESET ": %s\n", msg);
}
