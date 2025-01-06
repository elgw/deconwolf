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

#ifdef WINDOWS

#endif

int isdir(const char * dir)
{
#ifdef WINDOWS
    if( _taccess_s( dir, 0 ) == 0 )
    {
        struct _stat status;
        _tstat( dir, &status );
        return (status.st_mode & _S_IFDIR) != 0;
    }

    return 0;
#else

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
#endif
}



int ensuredir(const char * dir)

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

void dw_gettime(struct timespec * t)
{
#ifdef WINDOWS
    timespec_get(t, TIME_UTC); // since C11
#else
    clock_gettime(CLOCK_REALTIME, t);
#endif
    return;
}

char * dw_dirname(const char * path)
{
#ifdef WINDOWS
    size_t maxlen = strlen(path)+1;
    char * drive = malloc(maxlen);
    char * dir = malloc(maxlen);

    _splitpath(path, drive, dir, NULL, NULL);
    char * outpath = malloc(maxlen);
    _makepath(outpath, drive, dir, NULL, NULL);
    // Adds a trailing '\'
    free(drive);
    free(dir);
    return outpath;
#else
    char * t = strdup(path);
    char * _dir = dirname(t); // should not be freed
    // Does not add a trailing '/'
    char * dir = strdup(_dir);
    free(t);
    return dir;
#endif
}

char * dw_basename(const char * path)
{
#ifdef WINDOWS
    size_t maxlen = strlen(path)+1;
    char * fname = malloc(maxlen);
    char * ext = malloc(maxlen);
    char * outname = malloc(maxlen);
    _splitpath(path, NULL, NULL, fname, ext);
    _makepath(outname, NULL, NULL, fname, ext);
    free(fname);
    free(ext);
    return outname;
#else
    char * t = strdup(path);
    char * _fname = basename(t); // should not be freed
    char * fname = strdup(_fname);
    free(t);
    return fname;
#endif
}

char * dw_getcwd(char * buf, size_t size)
{
#ifdef WINDOWS
    return _getcwd(buf, size);
#else
    return getcwd(buf, size);
#endif
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
#else
    nThreads = atoi(getenv("NUMBER_OF_PROCESSORS"))/2;
#endif
#ifdef OMP
    /* Reports number of cores */
    nThreads = omp_get_num_procs();
#endif
    nThreads < 1 ? nThreads = 1 : 0;
    return nThreads;
}


int dw_isfile(const char * fname)
{
#ifndef WINDOWS
    struct stat buf;
    if(stat(fname, &buf) < 0) {
        return 0;
    }

    if(S_ISDIR(buf.st_mode) != 0) {
        return 0;
    }
    return 1;
#else
    FILE * fid = fopen(fname, "r");
    if(fid == NULL)
    {
        return 0;
    }
    fclose(fid);
    return 1;
#endif
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
                free(peakline);
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


float dw_read_scaling(const char * file)
{
    float scaling = 1.0;
    char * logfile = malloc(strlen(file)+32);
    assert(logfile != NULL);
    sprintf(logfile, "%s.log.txt", file);
    if( ! dw_isfile(logfile))
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


#ifdef WINDOWS
int getline(char **lineptr, size_t *n, FILE *stream)
{
    static char line[256];
    char *ptr;
    unsigned int len;

    if (lineptr == NULL || n == NULL)
    {
        errno = EINVAL;
        return -1;
    }

    if (ferror (stream))
        return -1;

    if (feof(stream))
        return -1;

    fgets(line,256,stream);

    ptr = strchr(line,'\n');
    if (ptr)
        *ptr = '\0';

    len = strlen(line);

    if ((len+1) < 256)
    {
        ptr = realloc(*lineptr, 256);
        if (ptr == NULL)
            return(-1);
        *lineptr = ptr;
        *n = 256;
    }

    strcpy(*lineptr,line);
    return(len);
}
#endif



char *
dw_prefix_file(const char * inFile, const char * prefix)
{
#ifdef WINDOWS

    char* drive = calloc(strlen(inFile) + 16, 1);
    char* dir = calloc(strlen(inFile) + 16, 1);
    char* fname = calloc(strlen(inFile) + 16, 1);
    char* ext = calloc(strlen(inFile) + 16, 1);

    _splitpath(
        inFile,
        drive,
        dir,
        fname,
        ext
    );

    char* pre_fname = calloc(strlen(fname) + strlen(prefix) + 16, 1);
    sprintf(pre_fname, "%s_%s", prefix, fname);
    char* outFile = calloc(strlen(inFile) + strlen(prefix) + 128, 1);

    _makepath(
        outFile,
        drive,
        dir,
        pre_fname,
        ext
    );

    free(drive);
    free(dir);
    free(fname);
    free(pre_fname);
    free(ext);
    return outFile;
#else
    char * dname = dw_dirname(inFile);
    assert(dname != NULL);
    char * fname = dw_basename(inFile);
    assert(fname != NULL);
    char * outFile = calloc(strlen(inFile) + strlen(prefix)+16, 1);
    assert(outFile != NULL);

    if(strlen(dname) > 0)
    {
        sprintf(outFile, "%s%c%s_%s", dname, FILESEP, prefix, fname);
    } else {
        sprintf(outFile, "%s_%s", prefix, fname);
    }

    free(dname);
    free(fname);

    return outFile;
#endif
}

float abbe_res_xy(float lambda, float NA)
{
    return lambda/(2.0*NA);
}


float abbe_res_z(float lambda, float NA)
{
    return 2.0*lambda / pow(NA, 2.0);
}
