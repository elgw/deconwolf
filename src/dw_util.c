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

int isdir(char * dir)
{

    DIR* odir = opendir(dir);
    if (odir) {
        /* Directory exists. */
        closedir(odir);
        return 1;
    } else if (ENOENT == errno) {
        /* Directory does not exist. */
        return 0;
    } else {
        /* opendir() failed for some other reason. */
        return 0;
    }
}



int ensuredir(char * dir)

{
    if(isdir(dir) == 1)
    {
        return 0;
    }

#ifdef WINDOWS
    if(_mkdir(dir) == ENOENT)
    {
        return 0;
    }
#else
    if(mkdir(dir, 0700) == 0)
    {
        return 0;
    }
#endif

    return 1;
}



int ptr_alignment_B(const void * p)
{
    int align = 1;
    size_t address = (size_t) p;
    if(address == 0)
    {
        return 0;
    }
    while(address % (align*2) == 0)
    {
        align*=2;
    }
    return align;
}

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
    assert(statfile != NULL);
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

    if(peakline == NULL)
    {
        return 0;
    }
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


float dw_read_scaling(char * file)
{
    float scaling = 1.0;
    char * logfile = malloc(strlen(file)+32);
    assert(logfile != NULL);
    sprintf(logfile, "%s.log.txt", file);
    if( ! dw_file_exist(logfile))
    {
        goto leave;
    }

    FILE * fid = fopen(logfile, "r");
    if(fid == NULL)
    {
        goto leave;
    }

    char * line = NULL;
    size_t len = 0;

    while( getline(&line, &len, fid) > 0)
    {
        if(strlen(line) > 8)
        {
            if(strncmp(line, "scaling:", 8) == 0)
            {
                scaling = atof(line + 8);
            }
        }
    }

    free(line);
    fclose(fid);
leave:
    free(logfile);
    return scaling;
}
