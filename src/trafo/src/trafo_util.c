#include "trafo_util.h"

double timespec_diff(struct timespec* end, struct timespec * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}

/* Assumes that the unit will be kb in /proc/self/status
 *
 * VmHWM = Peak Resident Set Size (peak RSS) i.e. how much physical memory
 *
 * VmPeak = peak Virtual Memory Size (peak VMZ) i.e. how much that is
 * allocated, although possibly never used.
 *
 */

int get_peakMemoryKB(size_t * _VmPeak, size_t * _VmHWM)
{
    FILE * sf = fopen("/proc/self/status", "r");
    if(sf == NULL)
    {
        fprintf(stderr, "Failed to open /proc/self/status\n");
        return 0;
    }
    *_VmPeak = 0;
    *_VmHWM = 0;

    char * VmPeak = NULL;
    char * VmHWM = NULL;

    char * line = NULL;
    size_t len = 0;

    while( getline(&line, &len, sf) > 0)
    {
        if(strlen(line) > 6)
        {
            if(strncmp(line, "VmPeak", 6) == 0)
            {
                free(VmPeak);
                VmPeak = strdup(line);
            }
            if(strncmp(line, "VmHWM", 5) == 0)
            {
                free(VmHWM);
                VmHWM = strdup(line);
            }
        }
    }
    free(line);
    fclose(sf);

    if((VmPeak != NULL) && (strlen(VmPeak) > 11))
    {
        VmPeak[strlen(VmPeak) - 4] = '\0';

        //    printf("peakline: '%s'\n", peakline+7);
        *_VmPeak = (size_t) atol(VmPeak+7);
    }
    free(VmPeak);

    if((VmHWM != NULL) && (strlen(VmHWM) > 10))
    {
        VmHWM[strlen(VmHWM) - 4] = '\0';

        //    printf("peakline: '%s'\n", peakline+7);
        *_VmHWM = (size_t) atol(VmHWM+7);
    }
    free(VmHWM);

    return 1;
}

void print_peak_memory(void)
{
    size_t VmPeak;
    size_t VmHWM;
    if(get_peakMemoryKB(&VmPeak, &VmHWM))
    {
        printf("VmPeak: %zu (kb) VmHWM: %zu (kb)\n", VmPeak, VmHWM);
    }
    return;
}

static void fill(int n, char c)
{
    for(int kk = 0; kk < n ; kk++)
        printf("%c", c);
}

static void line3(int w, const char * start, const char * fill, const char *end)
{
    printf("%s", start);
    for(int kk = 1; kk+1 < w; kk++)
    {
        printf("%s", fill);
    }
    printf("%s\n", end);
}


void print_section(const char * msg)
{

#if 1
    // Single line, curved corner
    const char nw[] = "\u256D"; // north west ...
    const char ne[] = "\u256e";
    const char sw[] = "\u2570";
    const char se[] = "\u256F";
    const char ho[] = "\u2500"; // horizontal
    const char ve[] = "\u2502"; // vertical
#else
    // Double lines, sharp corners
    const char nw[] = "\u2554"; // north west ...
    const char ne[] = "\u2557";
    const char sw[] = "\u255A";
    const char se[] = "\u255D";
    const char ho[] = "\u2550"; // horizontal
    const char ve[] = "\u2551"; // vertical
#endif

    // Text output width
    int w = 80;

    line3(w, nw, ho, ne);

    // Central line: ho w1 x ' ' msg w2 x ' ' ho
    int w1 = (w - 4 - strlen(msg))/2;
    int w2 = w - w1 - strlen(msg)-2;

    assert(1 + w1 + strlen(msg) + w2 + 1 == 80);
    printf("%s", ve);
    fill(w1, ' ');
    printf("%s", msg);
    fill(w2, ' ');
    printf("%s\n", ve);
    line3(w, sw, ho, se);
    printf("\n");
}
