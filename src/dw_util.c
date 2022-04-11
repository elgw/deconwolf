/*    Copyright (C) 2020 Erik L. G. Wernersson
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "dw_util.h"

int dw_get_threads(void)
{
    int nThreads = 4;
#ifndef WINDOWS
/* Reports #threads, typically 2x#cores */
    nThreads = sysconf(_SC_NPROCESSORS_ONLN)/2;
#endif
#ifdef OMP
/* Reports number of cores */
    nThreads = omp_get_num_procs();
#endif
    return nThreads;
}

int dw_file_exist(char * fname)
{
    if( access( fname, F_OK ) != -1 ) {
        return 1; // File exist
    } else {
        return 0;
    }
}

void dw_nullfree(void * p)
{
    if(p == NULL)
        return;

    free(p);
    return;
}

float timespec_diff(struct timespec* end, struct timespec * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}

#ifdef WINDOWS
size_t get_peakMemoryKB(void)
{
    return 0;
}
#else

#ifdef __APPLE__
size_t get_peakMemoryKB(void)
{
    struct rusage r_usage;
    getrusage(RUSAGE_SELF, &r_usage);
    return (size_t) round((double) r_usage.ru_maxrss/1024.0);
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
#endif
