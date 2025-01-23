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

#ifndef WINDOWS
#include <unistd.h>
#endif
#include "dw.h"

/* GLOBALS */
/* used by tiffErrHandler */
FILE * logfile = NULL;
/* END GLOBALS */


dw_iterator_t * dw_iterator_new(const dw_opts * s)
{
    dw_iterator_t * it = calloc(1, sizeof(dw_iterator_t));
    assert(it != NULL);
    it->error = 1;
    it->lasterror = 1;
    it->itertype = s->iter_type;
    it->abserror = s->err_abs;
    it->relerror = s->err_rel;
    it->iter = -1;
    it->niter = s->nIter;
    if(it->itertype != DW_ITER_FIXED)
    {
        it->niter = s->maxiter;
    }

    return it;
}

int dw_iterator_next(dw_iterator_t * it)
{
    it->iter++;
    switch(it->itertype)
    {
    case DW_ITER_FIXED:
        //printf("DW_ITER_FIXED\n");
        break;
    case DW_ITER_REL:
        //printf("DW_ITER_REL, %e, %e\n", fabs(it->error - it->lasterror)/it->error, it->relerror);
        if(it->iter < 2)
        {
            it->lasterror = 2*it->error*it->relerror;
        }
        if(fabs(it->error - it->lasterror)/it->error < it->relerror)
        {
            return -1;
        }
        break;
    case DW_ITER_ABS:
        //printf("DW_ITER_ABS %f < %f ?\n", it->error, it->abserror);
        if(it->iter > 0 && it->error < it->abserror)
        {
            return -1;
        }
        break;
    }

    if(it->iter >= it->niter)
    {
        return -1;
    }

    return it->iter;
}

void dw_iterator_set_error(dw_iterator_t * it, float err)
{
    it->lasterror = it->error;
    it->error = err;
}

void dw_iterator_show(dw_iterator_t * it, const dw_opts *s)
{
    if(s->verbosity > 0){
        printf("\r                                             ");
        if(s->metric == DW_METRIC_MSE)
        {
            printf("\rIteration %3d/%3d, fMSE=%.3e ",
                   it->iter+1, it->niter, it->error);
        }
        if(s->metric == DW_METRIC_IDIV)
        {
            printf("\rIteration %3d/%3d, Idiv=%.3e ",
                   it->iter+1, it->niter, it->error);
        }
        if(it->itertype == DW_ITER_REL && it->iter > 1)
        {
            double rel = fabs(it->error - it->lasterror)/it->error;
            printf("(%.3e", rel);
            if(rel > it->relerror)
            {
                printf(" > ");
            } else {
                printf(" < ");
            }
            printf("%.3e) ", it->relerror);
        }
        if(it->itertype == DW_ITER_ABS && it->iter > 0)
        {
            printf("(");
            if(it->error > it->abserror)
            {
                printf(" > ");
            } else {
                printf(" < ");
            }
            printf("%.3e) ", it->abserror);
        }
        fflush(stdout);
    }

    if(s->log != NULL && s->log != stdout)
    {
        if(s->metric == DW_METRIC_MSE)
        {
            fprintf(s->log, "Iteration %3d/%3d, fMSE=%e\n",
                    it->iter+1, it->niter, it->error);
        }
        if(s->metric == DW_METRIC_IDIV)
        {
            fprintf(s->log, "Iteration %3d/%3d, Idiv=%e\n",
                    it->iter+1, it->niter, it->error);
        }
        fflush(s->log);
    }
}
void dw_iterator_free(dw_iterator_t * it)
{
    free(it);
}



dw_opts * dw_opts_new(void)
{
    dw_opts * s = calloc(1, sizeof(dw_opts));
    assert(s != NULL);
    s->nThreads_FFT = dw_get_threads();
    s->nThreads_OMP = s->nThreads_FFT;

    s->nThreads_FFT < 1 ? s->nThreads_FFT = 1 : 0;
    s->nThreads_OMP < 1 ? s->nThreads_OMP = 1 : 0;

    s->fft_inplace = 1;

    s->nIter = 1; /* Always overwritten if used */
    s->maxiter = 250;
    s->err_rel = 0.02;
    s->err_abs = 1; /* Always overwritten if used */
    s->nIter_auto = 1;
    s->imFile = NULL;
    s->psfFile = NULL;
    s->outFile = NULL;
    s->logFile = NULL;
    s->refFile = NULL;
    s->tsvFile = NULL;
    s->iter_type = DW_ITER_REL;
    s->tsv = NULL;
    s->ref = NULL;
    s->prefix = malloc(10*sizeof(char));
    assert(s->prefix != NULL);
    sprintf(s->prefix, "dw");
    s->log = NULL;
    s->color = 1;
    s->verbosity = 1;
    s->showTime = 0;
    s->overwrite = 0;
    s->tiling_maxSize = -1;
    s->tiling_padding = 20;
    s->method = DW_METHOD_SHB;
    s->fun = deconvolve_shb;
    s->iterdump = 0;
    s->xycropfactor = 0.001;
    s->commandline = NULL;
    s->onetile = 0;
    s->borderQuality = 2;
    s->outFormat = 16; // write 16 bit int
    s->scaling = -1.0;
    s->experimental1 = 0;
    s->fulldump = 0;
    s->positivity = 1;
    s->bg = 1e-2; /* Should be strictly positive or pixels will be freezed */
    s->bg_auto = 1;
    s->offset = 5;
    s->flatfieldFile = NULL;
    s->lookahead = 0;
    s->psigma = 0;
    s->biggs = 1;
    s->metric = DW_METRIC_IDIV;
    dw_gettime(&s->tstart);
    s->fftw3_planning = FFTW_MEASURE;
    s->alphamax = 1;

    /* https://no-color.org/ */
    char *no_color = getenv("NO_COLOR");
    if (no_color != NULL && no_color[0] != '\0')
    {
        s->color = 0;
    }

    s->start_condition = DW_START_FLAT;

    return s;
}

void dw_show_iter(dw_opts * s, int it, int nIter, float err)
{
    if(s->verbosity > 0){
        printf("\r                                             ");
        if(s->metric == DW_METRIC_MSE)
        {
            printf("\rIteration %3d/%3d, fMSE=%e ", it+1, nIter, err);
        }
        if(s->metric == DW_METRIC_IDIV)
        {
            printf("\rIteration %3d/%3d, Idiv=%e ", it+1, nIter, err);
        }
        fflush(stdout);
    }

    if(s->log != NULL && s->log != stdout)
    {
        if(s->metric == DW_METRIC_MSE)
        {
            fprintf(s->log, "Iteration %3d/%3d, fMSE=%e\n", it+1, nIter, err);
        }
        if(s->metric == DW_METRIC_IDIV)
        {
            fprintf(s->log, "Iteration %3d/%3d, Idiv=%e\n", it+1, nIter, err);
        }
        fflush(s->log);
    }
}


char * gen_iterdump_name(
    __attribute__((unused)) const dw_opts * s,
    int it)
{
    // Generate a name for the an iterdump file
    // at iteration it
    char * name = malloc(strlen(s->outFolder) + 100*sizeof(char));
    assert(name != NULL);
    sprintf(name, "%sitd%05d.tif", s->outFolder, it);
    return name;
}

int64_t int64_t_max(int64_t a, int64_t b)
{
    if( a > b)
        return a;
    return b;
}

void dw_opts_free(dw_opts ** sp)
{
    dw_opts * s = sp[0];
    free(s->imFile);
    free(s->psfFile);
    free(s->outFile);
    free(s->outFolder);
    free(s->logFile);
    free(s->flatfieldFile);
    free(s->prefix);
    free(s->commandline);
    free(s->ref);
    free(s->refFile);
    free(s->tsvFile);
    if(s->tsv != NULL)
    {
        fclose(s->tsv);
    }
    free(s);
}

void dw_opts_fprint(FILE *f, dw_opts * s)
{
    f == NULL ? f = stdout : 0;

    fprintf(f, "> Settings:\n");
    fprintf(f, "image:  %s\n", s->imFile);
    if(s->flatfieldFile != NULL)
    {
        fprintf(f, "flat field: %s\n", s->flatfieldFile);
    }
    fprintf(f, "psf:    %s\n", s->psfFile);
    fprintf(f, "output: %s\n", s->outFile);
    fprintf(f, "log file: %s\n", s->logFile);
    fprintf(f, "nIter:  %d\n", s->nIter);
    fprintf(f, "nThreads for FFT: %d\n", s->nThreads_FFT);
    fprintf(f, "nThreads for OMP: %d\n", s->nThreads_OMP);
    fprintf(f, "verbosity: %d\n", s->verbosity);
    if(s->bg_auto == 1)
    {
        fprintf(f, "background level: auto\n");
    } else {
        fprintf(f, "background level: %f\n", s->bg);
    }

    switch(s->method)
    {
    case DW_METHOD_RL:
        fprintf(f, "method: Richardson-Lucy (RL)\n");
        break;
    case DW_METHOD_ID:
        fprintf(f, "method: Identity, doing nothing (ID)\n");
        break;
    case DW_METHOD_SHB:
        fprintf(f, "method: Scaled Heavy Ball (SHB)\n");
        break;
#ifdef OPENCL
    case DW_METHOD_SHBCL:
        fprintf(f, "method: Scaled Heavy Ball + OpenCL for FFT (SHBCL)\n");
        break;
    case DW_METHOD_SHBCL2:
        fprintf(f, "method: Scaled Heavy Ball + OpenCL (SHBCL2)\n");
        break;
#endif
    }
    switch(s->metric)
    {
    case DW_METRIC_MSE:
        fprintf(f, "metric: MSE\n");
        break;
    case DW_METRIC_IDIV:
        fprintf(f, "metric: Idiv\n");
        break;
    }
    switch(s->iter_type)
    {
    case DW_ITER_ABS:
        fprintf(f, "Stopping on absolute error: %e or %d iterations\n",
                s->err_abs, s->maxiter);
        break;
    case DW_ITER_REL:
        fprintf(f, "Stopping on relative error: %e or %d iterations\n",
                s->err_rel, s->maxiter);
        break;
    case DW_ITER_FIXED:
        fprintf(f, "Stopping after %d iterations\n",
                s->nIter);
        break;
    }
    if(s->psigma > 0)
    {
        fprintf(f, "pre-filtering enabled, sigma = %f\n", s->psigma);
    }
    if(s->overwrite == 1)
    { fprintf(f, "overwrite: YES\n"); } else
    { fprintf(f, "overwrite: NO\n"); }

    if(s->tiling_maxSize > 0)
    {
        fprintf(f, "tiling, maxSize: %d\n", s->tiling_maxSize);
        fprintf(f, "tiling, padding: %d\n", s->tiling_padding);
    } else {
        fprintf(f, "tiling: OFF\n");
    }
    fprintf(f, "XY crop factor: %f\n", s->xycropfactor);
    fprintf(f, "Offset: %f\n", s->offset);
    fprintf(f, "Output Format: ");
    switch(s->outFormat){
    case 16:
        fprintf(f, "16 bit integer\n");
        break;
    case 32:
        fprintf(f, "32 bit float\n");
        break;
    default:
        fprintf(f, "ERROR: Unknown\n");
        break;
    }
    if(s->outFormat == 16)
    {
        if(s->scaling <= 0)
        {
            fprintf(f, "Scaling: Automatic\n");
        } else {
            fprintf(f, "Scaling value: %f\n", s->scaling);
        }

    }

    fprintf(f, "Border Quality: ");
    switch(s->borderQuality){
    case 0:
        fprintf(f, "0 No boundary handling\n");
        break;
    case 1:
        fprintf(f, "1 Somewhere between 0 and 2\n");
        break;
    case 2:
        fprintf(f, "2 Minimal boundary artifacts\n");
        break;
    default:
        ;
    }
    fprintf(f, "FFT lookahead: %d", s->lookahead);

    if(s->onetile == 1)
    {
        fprintf(f, "DEBUG OPTION: ONETILE = TRUE (only first tile will be deconvolved)\n");
    }
    fprintf(f, "\n");

    fprintf(f, "FFTW3 plan: ");
    switch(s->fftw3_planning)
    {
    case FFTW_ESTIMATE:
        fprintf(f, "FFTW_ESTIMATE (--noplan)\n");
        break;
    case FFTW_MEASURE:
        fprintf(f, "FFTW_MEASURE\n");
        break;
    case FFTW_PATIENT:
        fprintf(f, "FFTW_PATIENT\n");
        break;
    case FFTW_EXHAUSTIVE:
        fprintf(f, "FFTW_EXHAUSTIVE\n");
        break;
    default:
        fprintf(f, "UNKNOWN/Wrongly Set\n");
        break;
    }

    fprintf(f, "Initial guess: ");
    switch(s->start_condition)
    {
    case DW_START_FLAT:
        fprintf(f, "Flat average\n");
        break;
    case DW_START_IDENTITY:
        fprintf(f, "Identity\n");
        break;
    case DW_START_LP:
        fprintf(f, "Low pass filtered\n");
    }
    return;
}

void warning(FILE * fid)
{
    //fprintf(fid, ANSI_UNDERSCORE " ! " ANSI_COLOR_RESET );
    fprintf(fid, " ! ");
    return;
}


void fulldump(dw_opts * s, float * A, size_t M, size_t N, size_t P, char * name)
{
    /* Write A to disk as fulldump_<name>.tif if s->fulldump = 1 */
    if(s->fulldump != 1)
    {
        return;
    }
    assert(name != NULL);

    if(A != NULL)
    {
        printf("Dumping to %s\n", name);
        fim_tiff_write_float(name, A, NULL, M, N, P);
    }
    return;
}

void dw_fprint_info(FILE * f, dw_opts * s)
{
    f == NULL ? f = stdout : 0;

    fprintf(f, "deconwolf: '%s'\n", deconwolf_version);

    if(f != stdout)
    {
#ifndef WINDOWS
        fprintf(f, "PID: %d\n",  (int) getpid());
#endif
    }

    if(f != stdout)
    {
        char cwd[1024];
        if (dw_getcwd(cwd, sizeof(cwd)) != NULL) {
            fprintf(f, "PWD: %s\n", cwd);
        }
    }

    if(f != stdout)
    {
        if(s->commandline != NULL)
        {
            fprintf(f, "CMD: %s\n", s->commandline);
        }
    }

#ifdef GIT_VERSION
    fprintf(f, "GIT_VERSION: '%s'\n", GIT_VERSION);
#endif

#ifdef CC_VERSION
    fprintf(f, "COMPILER: '%s'\n", CC_VERSION);
#endif

    fprintf(f, "BUILD_DATE: '%s'\n", __DATE__);
#ifdef CUDA
    fprintf(f, "FFT Backend: 'cuFFT\n");
#else
#ifndef WINDOWS
    fprintf(f, "FFT Backend: '%s'\n", fftwf_version);
#endif
#endif
    fprintf(f, "TIFF Backend: '%s'\n", TIFFGetVersion());

    if(f != stdout)
    {
#ifndef WINDOWS
        char * user = getenv("USER");
        if(user != NULL)
        { fprintf(f, "USER: '%s'\n", user); }
#endif

#ifndef WINDOWS
        char * hname = malloc(1024*sizeof(char));
        assert(hname != NULL);
        if(gethostname(hname, 1023) == 0)
        {
            fprintf(f, "HOSTNAME: '%s'\n", hname);
            free(hname);
        }
#endif
    }

#ifdef _OPENMP
    fprintf(f, "OpenMP: YES\n");
#endif

#ifdef OPENCL
    fprintf(f, "OpenCL: YES\n");
#else
    fprintf(f, "OpenCL: No. (Rebuild with OPENCL=1 to enable)\n");
#endif

#ifdef VKFFT
    fprintf(f, "VkFFT: YES\n");
#endif

#ifndef NDEBUG
    fprintf(f, "Warning: Built in debug mode (NDEBUG not defined)\n");
#endif

    if(s->verbosity > 1)
    {
        fprintf(f, "sizeof(int) = %zu\n", sizeof(int));
        fprintf(f, "sizeof(float) = %zu\n", sizeof(float));
        fprintf(f, "sizeof(double) = %zu\n", sizeof(double));
        fprintf(f, "sizeof(size_t) = %zu\n", sizeof(size_t));
    }

    fprintf(f, "\n");
    fflush(f);
    return;
}

static void
getCmdLine(int argc, char ** argv, dw_opts * s)
{
    // Copy the command line to s->commandline
    size_t lcmd = 2;
    for(int kk = 0; kk<argc; kk++)
    {
        lcmd += strlen(argv[kk]) + 4;
    }

    s->commandline = calloc(lcmd, 1);
    assert(s->commandline != NULL);
    snprintf(s->commandline, lcmd, "%s", "");

    for(int kk = 0; kk<argc; kk++)
    {
        size_t pos = strlen(s->commandline);
        snprintf(s->commandline+pos, lcmd, "'%s' ", argv[kk]);
    }

    return;
}


void dw_argparsing(int argc, char ** argv, dw_opts * s)
{

    getCmdLine(argc, argv, s);


    struct option longopts[] = {
        { "no-inplace",no_argument,       NULL, '1' },
        { "ompthreads",required_argument, NULL, '2' },
        { "psf-pass",  required_argument, NULL, '3' },
	{ "cldevice",  required_argument, NULL, '4' },
        { "start_flat",no_argument,       NULL, '5' },
        { "start_id",  no_argument,       NULL, '6' },
        { "start_lp",  no_argument,       NULL, '7' },
        { "temp",      required_argument, NULL, '9' },
        { "noplan",    no_argument,       NULL, 'a' },
        { "bg",        required_argument, NULL, 'b' },
        { "threads",   required_argument, NULL, 'c' },
        { "tsv",       required_argument, NULL, 'd' },
        { "abserror",  required_argument, NULL, 'e' },
        { "prefix",    required_argument, NULL, 'f' },
        { "times",     no_argument,       NULL, 'g' },
        { "help",      no_argument,       NULL, 'h' },
        { "iterdump",  no_argument,       NULL, 'i' },
        { "relerror",  required_argument, NULL, 'j' },
        { "verbose",   required_argument, NULL, 'l' },
        { "method",    required_argument, NULL, 'm' },
        { "iter",      required_argument, NULL, 'n' },
        { "out",       required_argument, NULL, 'o' },
        { "tilepad",   required_argument, NULL, 'p' },
        { "offset",    required_argument, NULL, 'q' },
        { "tilesize",  required_argument, NULL, 's' },
        { "test",      no_argument,       NULL, 't' },
        { "version",   no_argument,       NULL, 'v' },
        { "overwrite", no_argument,       NULL, 'w' },
        { "xyfactor",  required_argument, NULL, 'x' },
        { "az",        required_argument, NULL, 'A' },
        { "bq",        required_argument, NULL,  'B' },
        { "flatfield", required_argument, NULL,  'C' },
        { "fulldump",  no_argument,       NULL,  'D' },
        { "float",     no_argument,       NULL,  'F' },
        { "gpu",       no_argument,       NULL,  'G' },
        { "niterdump", required_argument, NULL,  'I' },
        { "lookahead", required_argument, NULL,  'L' },
        { "mse",       no_argument,       NULL,  'M' },
        { "maxiter",   required_argument, NULL,  'N' },
        { "periodic",  no_argument,       NULL,  'O' },
        { "ref",       required_argument, NULL,  'R' },
        { "scaling",   required_argument, NULL,  'S' },
        { "onetile",   no_argument,       NULL,  'T' },
        { "nopos",     no_argument,       NULL,  'P' },
        { "psigma",    required_argument, NULL,  'Q' },
        { "expe1",     no_argument,       NULL,  'X' },
        { "cz",        required_argument, NULL,  'Z' },
        { NULL,           0,                 NULL,   0   }
    };


    int known_method = 1;
    int ch;
    int prefix_set = 0;
    int use_gpu = 0;
    while((ch = getopt_long(argc, argv,
                            "12345679ab:c:f:Gghil:m:n:o:p:q:s:tvwx:A:B:C:DFI:L:MOR:S:TPQ:X:Z:",
                            longopts, NULL)) != -1)
    {
        switch(ch) {
        case '1':
            s->fft_inplace = 0;
            break;
        case '2':
            s->nThreads_OMP = atoi(optarg);
            break;
        case '3':
            s->psf_pass = atof(optarg);
            break;
	case '4':
            s->cl_device = atoi(optarg);
            break;
        case '5':
            s->start_condition = DW_START_FLAT;
            break;
        case '6':
            s->start_condition = DW_START_IDENTITY;
            break;
        case '7':
            s->start_condition = DW_START_LP;
            break;
        case '9':
            s->alphamax = atof(optarg);
            break;
        case 'C':
            s->flatfieldFile = strdup(optarg);
            break;
        case 'a':
            s->fftw3_planning = FFTW_ESTIMATE;
            break;
        case 'A':
            s->auto_zcrop = atoi(optarg);
            break;
        case 'b':
            s->bg = atof(optarg);
            s->bg_auto = 0;
            break;
        case 'd': /* diagnostics */
            s->tsvFile = strdup(optarg);
            break;
        case 'e': /* --abserror */
            s->err_abs = atof(optarg);
            s->iter_type = DW_ITER_ABS;
            break;
        case 'g':
            s->showTime = 1;
            break;
        case 'G':
            use_gpu = 1;
            break;
        case 'j': /* --relerror */
            s->err_rel = atof(optarg);
            s->iter_type = DW_ITER_REL;
            break;
        case 'L':
            s->lookahead = atoi(optarg);
            break;
        case 'M':
            s->metric = DW_METRIC_MSE;
            break;
        case 'N':
            s->maxiter = atoi(optarg);
            break;
        case 'O':
            s->borderQuality = 0;
            break;
        case 'P':
            s->positivity = 0;
            printf("Turning off positivity constraint!\n");
            break;
        case 'q':
            s->offset = atof(optarg);
            break;
        case 'Q':
            s->psigma = atof(optarg);
            break;
        case 'D':
            s->fulldump = 1;
            break;
        case 'F':
            s->outFormat = 32;
            break;
        case 'B':
            s->borderQuality = atoi(optarg);
            break;
        case 'v':
            dw_fprint_info(NULL, s);
            exit(0);
            break;
        case 'h':
            dw_usage(argc, argv, s);
            exit(0);
            break;
        case 'o':
            free(s->outFile);
            s->outFile = strdup(optarg);
            assert(s->outFile != NULL);
            break;
        case 'n':
            s->nIter = atoi(optarg);
            s->iter_type = DW_ITER_FIXED;
            break;
        case 'c':
            s->nThreads_FFT = atoi(optarg);
            s->nThreads_OMP = atoi(optarg);
            break;
        case 'i':
            s->iterdump = 1;
            break;
        case 'I':
            s->iterdump = atoi(optarg);
            break;
        case 'l':
            s->verbosity = atoi(optarg);
            break;
        case 't':
            dw_unittests();
            exit(0);
            break;
        case 's':
            s->tiling_maxSize = atoi(optarg);
            break;
        case 'S':
            s->scaling = atof(optarg);
            break;
        case 'p':
            s->tiling_padding = atoi(optarg);
            break;
        case 'w':
            s->overwrite = 1;
            break;
        case 'f':
            free(s->prefix);
            s->prefix = malloc(strlen(optarg) + 1);
            assert(s->prefix != NULL);
            strcpy(s->prefix, optarg);
            prefix_set = 1;
            break;
        case 'm':
            known_method = 0;
            if(strcmp(optarg, "rl") == 0)
            {
                s->method = DW_METHOD_RL;
                if(prefix_set == 0)
                {
                    sprintf(s->prefix, "rl");
                }
                s->fun = &deconvolve_rl;
                known_method = 1;
            }
            if(strcmp(optarg, "id") == 0)
            {
                s->method = DW_METHOD_ID;
                if(prefix_set == 0)
                {
                    sprintf(s->prefix, "id");
                }
                s->fun = &deconvolve_identity;
                known_method = 1;
            }
            if(strcmp(optarg, "shb") == 0)
            {
                s->method = DW_METHOD_SHB;
                if(prefix_set == 0)
                {
                    sprintf(s->prefix, "shb");
                }
                s->fun = &deconvolve_shb;
                known_method = 1;
            }
#ifdef OPENCL
            if(strcmp(optarg, "shbcl") == 0)
            {
                s->method = DW_METHOD_SHBCL;
                if(prefix_set == 0)
                {
                    sprintf(s->prefix, "shbcl");
                }
                s->fun = &deconvolve_shb_cl;
                known_method = 1;
            }
            if(strcmp(optarg, "shbcl2") == 0)
            {
                s->method = DW_METHOD_SHBCL2;
                if(prefix_set == 0)
                {
                    sprintf(s->prefix, "shbcl2");
                }
                s->fun = &deconvolve_shb_cl2;
                known_method = 1;
            }
#endif
            if(known_method == 0)
            {
                fprintf(stderr, "--method %s is unknown. Please specify ",  optarg);
#ifdef OPENCL
                fprintf(stderr, "shbcl, shbcl2, ");
#endif
                fprintf(stderr, "shb (default), rl or id\n");
                exit(EXIT_FAILURE);
            }
            break;
        case 'x':
            s->xycropfactor = atof(optarg);
            if(s->xycropfactor > 1 || s->xycropfactor < 0)
            {
                fprintf(stderr, "The crop factor in x and y has to be => 0 and < 1\n");
                exit(1);
            }
            break;
        case 'R':
            s->refFile = strdup(optarg);
            break;
        case 'T':
            s->onetile = 1;
            break;
        case 'X':
            s->experimental1 = 1;
            break;
        case 'Z':
            s->zcrop = atoi(optarg);
            break;
        default:
            fprintf(stderr, "dw got an unknown command line argument. Exiting!\n");
            exit(EXIT_FAILURE);
        }
    }


    if(s->offset < 0)
    {
        printf("WARNING: A negative offset not allowed, setting it to 0\n");
        s->offset = 0;
    }

    if(s->zcrop < 0)
    {
        fprintf(stderr,
                "the zcrop value (--cz) can not be negative\n");
        exit(EXIT_FAILURE);
    }

    if(s->zcrop > 0)
    {
        if(s->auto_zcrop > 0)
        {
            fprintf(stderr, "zcrop and auto_zcrop can not be combined\n");
            exit(EXIT_FAILURE);
        }
        if(s->tiling_maxSize > 0)
        {
            fprintf(stderr, "zcrop can not be combined with tiling\n");
        }
    }

    if( (s->auto_zcrop > 0) && (s->tiling_maxSize > 0))
    {
        fprintf(stderr, "auto_zcrop can not be combined with tiling\n");
        exit(EXIT_FAILURE);
    }

    if(use_gpu)
    {
#ifdef OPENCL
        if(s->method == DW_METHOD_SHB)
        {
            s->method = DW_METHOD_SHBCL2;
            s->fun = &deconvolve_shb_cl2;
        }
#else
        printf("WARNING: dw was not compiled with GPU support\n");
#endif
    }

    /* Take care of the positional arguments */
    if(optind + 2 != argc)
    {
        printf("Deconwolf: To few input arguments.\n");
        printf("See `%s --help` or `man dw`.\n", argv[0]);
        exit(1);
    }

#ifdef WINDOWS
    /* TODO, see GetFullPathNameA in fileapi.h */
    s->imFile = strdup(argv[optind]);
#else
    s->imFile = realpath(argv[optind], 0);
#endif

    if(s->imFile == NULL)
    {
        fprintf(stderr, "ERROR: Can't read %s\n", argv[optind]);
        exit(1);
    }

#ifdef WINDOWS
    s->psfFile = strdup(argv[++optind]);
#else
    s->psfFile = realpath(argv[++optind], 0);
#endif

    if(s->psfFile == NULL)
    {
        fprintf(stderr, "ERROR: Can't read %s\n", argv[optind]);
        exit(1);
    }


    /* Set s->outFile and s->outFolder based on s->imFile */
    if(s->outFile == NULL)
    {
        s->outFile = dw_prefix_file(s->imFile, s->prefix);

        char * dname = dw_dirname(s->imFile);
        s->outFolder = malloc(strlen(dname) + 16);
        assert(s->outFolder != NULL);
        sprintf(s->outFolder, "%s%c", dname, FILESEP);
        free(dname);
        if(s->verbosity > 1)
        {
            printf("outFile: %s, outFolder: %s\n", s->outFile, s->outFolder);
        }
    } else {
        if( isdir(s->outFile) )
        {
            if(s->verbosity > 0)
            {
                printf("--out describes a folder");
            }
            free(s->outFolder);
            s->outFolder = malloc(strlen(s->outFile) + 8);
            sprintf(s->outFolder, "%s%c", s->outFile, FILESEP);
            free(s->outFile);
            char * basename = dw_basename(s->imFile);
            char * outfile0 = malloc(strlen(basename) + strlen(s->outFolder) + 8);
            assert(outfile0 != NULL);
            sprintf(outfile0, "%s%c%s", s->outFolder, FILESEP, basename);
            s->outFile = dw_prefix_file(outfile0, s->prefix);
            free(outfile0);
        } else {
            char * dname = dw_dirname(s->outFile);
            free(s->outFolder);
            s->outFolder = malloc(strlen(dname) + 16);
            assert(s->outFolder != NULL);
            sprintf(s->outFolder, "%s%c", dname, FILESEP);
            free(dname);
        }
    }

    if(! s->iterdump)
    {
        if( s->overwrite == 0 && dw_isfile(s->outFile))
        {
            printf("%s already exist. Use --overwrite to overwrite existing files.\n",
                   s->outFile);
            exit(0);
        }
    }

    if(s->nThreads_FFT < 1 || s->nThreads_OMP < 1)
    {
        printf("Invalid number of threads (%d), "
               "please verify your command line\n", s->nThreads_FFT);
        exit(EXIT_FAILURE);
    }

    s->logFile = malloc(strlen(s->outFile) + 10);
    assert(s->logFile != NULL);
    sprintf(s->logFile, "%s.log.txt", s->outFile);



    if(s->tsvFile != NULL)
    {
        s->tsv = fopen(s->tsvFile, "w");
        if(s->tsv == NULL)
        {
            fprintf(stderr, "Failed to open %s for writing\n", s->tsvFile);
            exit(EXIT_FAILURE);
        }
        fprintf(s->tsv, "iteration\ttime\tKL\n");
    }

    /* Set the plan to be used with fftw3 */
    fft_set_plan(s->fftw3_planning);
    fft_set_inplace(s->fft_inplace);
    if(s->verbosity > 2)
    {
        printf("Command line arguments accepted\n");
    }
}


void fsetzeros(const char * fname, size_t N)
/* Create a new file and fill it with N bytes of zeros
 */
{
    size_t bsize = 1024*1024;
    char * buffer = malloc(bsize);
    assert(buffer != NULL);
    memset(buffer, 0, bsize);
    FILE * fid = fopen(fname, "wb");
    if(fid == NULL)
    {
        fprintf(stderr, "%s (%d): Unable to open %s\n",
                __FILE__, __LINE__, fname);
        exit(EXIT_FAILURE);
    }
    size_t written = 0;
    while(written + bsize < N)
    {
        fwrite(buffer, bsize, 1, fid);
        written += bsize;
    }
    //  printf("\r %zu / %zu \n", written, N); fflush(stdout);

    fwrite(buffer, N-written, 1, fid);
    fclose(fid);
    free(buffer);
}

void benchmark_write(dw_opts * s, int iter, double fMSE,
                     const float * x0, // current guess of work size
                     const int64_t M, const int64_t N, const int64_t P,
                     const int64_t wM, const int64_t wN, const int64_t wP)
{
    if(s->tsv == NULL)
    {
        return;
    }

    float * x = fim_subregion(x0, wM, wN, wP, M, N, P);

    size_t MNP = M*N*P;
    double KL = 0;
    if(s->ref != NULL)
    {
        for(size_t kk = 0; kk<MNP; kk++)
        {
            if(s->ref[kk] > 0)
            {
                KL += log( x[kk]/s->ref[kk]) * x[kk];
            }
        }
    }
    free(x);
    struct timespec tnow;
    dw_gettime(&tnow);
    double time = clockdiff(&tnow, &s->tstart);
    fprintf(s->tsv, "%d\t%f\t%f\t%f\n", iter, time, fMSE, KL);

    return;
}

float getErrorX(const float * restrict y, const float * restrict g, const int64_t M, const int64_t N, const int64_t P, const int64_t wM, const int64_t wN, const int64_t wP)
{
    /* Same as getError with the difference that G is expanded to MxNxP */
    if(M > wM || N > wN || P > wP)
    {
        fprintf(stderr,"Something is funky with the dimensions of the images.\n");
        exit(-1);
    }

    double e = 0;
    for(size_t c = 0; c < (size_t) P; c++)
    {
        for(size_t b = 0; b < (size_t) N; b++)
        {
            for(size_t a = 0; a < (size_t) M; a++)
            {
                double yval = y[a + b*wM + c*wM*wN];
                double gval = g[a + b*wM + c*wM*wN];
                e+=pow(yval-gval, 2);
            }
        }
    }
    e/=(M*N*P);
    return (float) e;
}


float getError_idiv(const float * restrict y, const float * restrict g,
                    const int64_t M, const int64_t N, const int64_t P,
                    const int64_t wM, const int64_t wN, const int64_t wP)
{
    /* Csiszar’s I-Divergence between the input, y, and the forward propagated
     * current guess */

    if(M > wM || N > wN || P > wP)
    {
        fprintf(stderr, "Something is funky with the dimensions of the images.\n");
        exit(-1);
    }
    float idiv = 0;
#pragma omp parallel for reduction(+: idiv)
    for(int64_t c = 0; c<P; c++)
    {
        for(int64_t b = 0; b<N; b++)
        {
            for(int64_t a = 0; a<M; a++)
            {
                float obs =  y[a + b*wM + c*wM*wN];
                float est = g[a + b*M + c * M*N];
                if(est > 0)
                {
                    idiv += obs*logf(obs/est) - obs + est;
                }
            }
        }
    }

    float idiv_mean = idiv / (float) (M*N*P);
    return idiv_mean;
}

/* Return the "error" or distance between the input image and the
   current guess convolved with the PSF. Also known as the forward
   error */
float getError(const float * restrict y, const float * restrict g,
               const int64_t M, const int64_t N, const int64_t P,
               const int64_t wM, const int64_t wN, const int64_t wP,
               dw_metric metric)
{
    /* Idiv is a better distance measurement but it takes longer to
       calculate due to the log function.
       These calculations needs to be done with double precision or the results
       will look different depending on the number of threads used (cancellation)
    */

    float error = 0;
    switch(metric)
    {
    case DW_METRIC_MSE:
        error = get_fMSE(y, g, M, N, P, wM, wN, wP);
        break;
    case DW_METRIC_IDIV:
        error = get_fIdiv(y, g, M, N, P, wM, wN, wP);
        break;
    }

    return error;
}

float get_fMSE(const float * restrict y, const float * restrict g,
               const int64_t M, const int64_t N, const int64_t P,
               const int64_t wM, const int64_t wN, const int64_t wP)
{
    /* Mean squared error between the input, y, and the forward propagated
     * current guess */

    if(M > wM || N > wN || P > wP)
    {
        fprintf(stderr, "Something is funky with the dimensions of the images.\n");
        exit(-1);
    }
    double e = 0;
    //#pragma omp parallel for reduction(+: e) shared(y, g)
    for(int64_t c = 0; c<P; c++)
    {
        for(int64_t b = 0; b<N; b++)
        {
            for(int64_t a = 0; a<M; a++)
            {
                double yval = y[a + b*wM + c*wM*wN];
                double gval = g[a + b*M + c * M*N];
                e += pow(yval-gval, 2);
            }
        }
    }

    double mse = e / (double) (M*N*P);
    return (double) mse;
}

float get_fIdiv(const float * restrict y, const float * restrict g,
                const int64_t M, const int64_t N, const int64_t P,
                const int64_t wM, const int64_t wN, const int64_t wP)
{
    /* Mean squared error between the input, y, and the forward propagated
     * current guess */

    if(M > wM || N > wN || P > wP)
    {
        fprintf(stderr, "Something is funky with the dimensions of the images.\n");
        exit(-1);
    }
    double I = 0;
#pragma omp parallel for reduction(+: I) shared(y, g)
    for(int64_t c = 0; c<P; c++)
    {
        for(int64_t b = 0; b<N; b++)
        {
            for(int64_t a = 0; a<M; a++)
            {
                double yval = y[a + b*wM + c*wM*wN];
                double gval = g[a + b*M + c * M*N];
                if(yval > 0 && gval > 0)
                {
                    I += gval*logf(gval/yval) - (gval-yval);
                }
            }
        }
    }

    double mI = I / (double) (M*N*P);
    return (float) mI;
}


void putdot(const dw_opts *s)
{
    if(s->verbosity > 0)
    {
        if(s->color)
        {
            printf(ANSI_COLOR_GREEN "." ANSI_COLOR_RESET);
        } else {
            printf(".");
        }
        fflush(stdout);
    }
    return;
}

float getError_ref(const float * restrict y,
                   const float * restrict g,
                   int64_t M, int64_t N, int64_t P,
                   int64_t wM, int64_t wN, int64_t wP)
{
    /* Note: this is just for comparison with getError, not used
     * in production
     **/
    if(M > wM || N > wN || P > wP)
    {
        fprintf(stderr, "Something is funky with the dimensions of the images.\n");
        exit(-1);
    }


    double e = 0;
    for(int64_t c = 0; c<P; c++)
    {
        for(int64_t b = 0; b<N; b++)
        {
            for(int64_t a = 0; a<M; a++)
            {
                double yval = y[a + b*wM + c*wM*wN];
                double gval = g[a + b*M + c * M*N];
                e+=pow(yval-gval, 2);
            }
        }
    }

    e/=(M*N*P);
    return (float) e;
}



void dw_usage(__attribute__((unused)) const int argc, char ** argv, const dw_opts * s)
{
    printf("deconwolf: %s\n", deconwolf_version);
    printf("usage: %s [<options>] image.tif psf.tif\n", argv[0]);

    printf("\n");
    printf(" Options:\n");
    printf(" --version\n\t Show version info\n");
    printf(" --help\n\t Show this message\n");
    printf(" --out file\n\t Specify output image name. If not set the input image "
           "will be prefixed with dw_\n.");
    printf(" --iter N\n\t "
           "Specify the number of iterations to use (default: %d)\n", s->nIter);
    printf(" --gpu\n\t Use GPU processing\n");
    printf(" --cldevice n\n\t Use OpenCL device #n\n");
    printf(" --threads N\n\t Specify the number of CPU threads to use\n");
    printf(" --verbose N\n\t Set verbosity level (default: %d)\n", s->verbosity);
    printf(" --test\n\t Run unit tests\n");
    printf(" --tilesize N\n\t Enables tiling mode and sets the largest tile size\n\t"
           "to N voxels in x and y.\n");
    printf(" --tilepad N\n\t Sets the tiles to overlap by N voxels in tile mode \n\t"
           "(default: %d)\n", s->tiling_padding);
    printf(" --prefix str\n\t Set the prefix of the output files (default: '%s')\n",
           s->prefix);
    printf(" --overwrite\n\t Allows deconwolf to overwrite already existing output files\n");

    printf(" --psigma s\n\t"
           "Pre filter the image with a Gaussian filter of sigma=s while Anscombe "
           "transformed\n");

    printf(" --cz n\n\t"
           "remove n zplanes from the top and bottom of the image.\n\t"
           "has no effect when tiling is used\n");

    printf(" --az n\n\t"
           "Automatically crop the image to n z-planes. The plane with the largest\n\t"
           "will be placed in the centre.\n\t"
           "has no effect when tiling is used\n");

    printf(" --bq Q\n\t Set border handling to \n\t"
           "0 'none' i.e. periodic\n\t"
           "1 'compromise', or\n\t"
           "2 'normal' which is default\n");
    printf("--periodic\n\t Equivalent to --bq 0\n");

    printf("--scale s\n\t"
           "Set the scaling factor for the output image manually to s.\n\t"
           "Warning: Might cause clipping or discretization artifacts\n\t"
           "This option is only used when the output format is 16-bit (i.e. \n\t"
           "not with --float)\n");
    printf(" --float\n\t Set output format to 32-bit float (default is 16-bit \n\t"
           "int) and disable scaling\n");
    printf(" --bg l\n\t Set background level, l\n");
    printf(" --offset l\n\t"
           "Set a positive offset that will be added to the image during\n\t"
           "processing and remove before saving to disk. Can help to mitigate\n\t"
           "some of the detector noise (non-Poissonian)\n");
    printf(" --flatfield image.tif\n\t"
           " Use a flat field correction image. Deconwolf will divide each plane of the\n\t"
           " input image, pixel by pixel, by this correction image.\n");
    printf(" --lookahead N"
           "\n\t Try to do a speed-for-memory trade off by using a N pixels larger"
           "\n\t job size that is better suited for FFT.\n");
    printf("--method name\n\t"
           "Select what method to use. Valid options: rl, shb, shbcl2\n");
    printf("--start_id\n\t"
           "Set the input image as the initial guess\n");
    printf("--start_flat\n\t"
           "Use the average of the input image as the initial guess. Default\n");
    printf("--start_lp\n\t"
           "Use a low passed version of the input image as the initial guess.n");
    printf(" --noplan\n\t Don't use any planning optimization for fftw3\n");
    printf(" --no-inplace\n\t Disable in-place FFTs (for fftw3), uses more\n\t"
           "memory but could potentially be faster for some problem sizes.\n");
    printf("\n");

    printf("Additional commands with separate help sections:\n");
    printf("   maxproj      maximum Z-projections\n");
    printf("   merge        merge individual slices to volume\n");
#ifdef dw_module_dots
    printf("   dots         detect dots with sub pixel precision\n");
#endif
#ifdef dw_module_psf
    printf("   psf          generate PSFs for Widefield and Confocal\n");
#endif
#ifdef dw_module_psf_sted
    printf("   psf-STED     PSFs for 3D STED\n");
#endif
#ifdef dw_module_nuclei
    printf("   nuclei       pixel classifier\n");
#endif
#ifdef dw_module_background
    printf("   background   vignetting/background estimation\n");
#endif
    printf("\n");
    printf("see: %s [command] --help\n", argv[0]);
    printf("\n");

    printf("Web page: https://www.github.com/elgw/deconwolf/\n");
}




double clockdiff(struct timespec* end, struct timespec * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}


void testfinite(float * x, size_t N)
{
    for(size_t kk = 0; kk<N; kk++)
    {
        if(!isfinite(x[kk]))
        {
            printf("Not finite\n");
            exit(1);
        }
    }
}





float * psf_autocrop_centerZ(float * psf, int64_t * pM, int64_t * pN, int64_t * pP,  // psf and size
                             dw_opts * s)
{

    const int64_t m = pM[0];
    const int64_t n = pN[0];
    const int64_t p = pP[0];

    const int64_t midm = (m-1)/2;
    const int64_t midn = (n-1)/2;
    const int64_t midp = (p-1)/2;

    //  printf("m: %d, n:%d, p:%d\n", m, n, p);
    //  printf("midm: %d, midn: %d, midp: %d\n", midm, midn, midp);

    float maxvalue = -INFINITY;
    int64_t maxp = -1;

    for(int64_t pp = 0; pp<p; pp++)
    {
        size_t idx = midm + midn*m + pp*m*n;
        if(psf[idx] > maxvalue)
        {
            maxp = pp;
            maxvalue = psf[idx];
        }
    }

    if(maxp == midp)
    {
        if(s->verbosity > 2)
        {
            printf("PSF is Z-centered :)\n");
        }
        return psf;
    }


    int64_t m0 = 0, m1 = m-1;
    int64_t n0 = 0, n1 = n-1;
    int64_t p0 = maxp, p1 = maxp;

    while(p0 > 1 && p1+2 < p)
    {
        p0--; p1++;
    }
    if(s->verbosity > 0)
    {
        printf("PSF has %" PRId64 " slices\n", p);
        printf("brightest at plane %" PRId64 "\n", maxp);
        printf("Selecting Z-planes: %" PRId64 " -- %" PRId64 "\n", p0, p1);
    }

    fprintf(s->log, "Selecting Z-planes %" PRId64 " -- %" PRId64 "\n", p0, p1);

    float * psf_cropped = fim_get_cuboid(psf, m, n, p,
                                         m0, m1, n0, n1, p0, p1);
    free(psf);
    pP[0] = p1-p0+1;
    return psf_cropped;

}

float * psf_autocrop_byImage(float * psf,/* psf and size */
                             int64_t * pM, int64_t * pN, int64_t * pP,
                             int64_t M, int64_t N, int64_t P, /* image size */
                             dw_opts * s)
{

    /* Purpose:
     * For bq=0
     *  Crop the PSF so it is not larger than the image.
     * For bq > 0
     *   Crop the PSF if it is more
     *   than 2x the size of the image in any dimension.
     */

    const int64_t m = pM[0];
    const int64_t n = pN[0];
    const int64_t p = pP[0];


    /* Optimal size */
    int64_t mopt = (M-1)*2 + 1;
    int64_t nopt = (N-1)*2 + 1;
    int64_t popt = (P-1)*2 + 1;

    if(s->borderQuality == 0)
    {
        mopt = M;
        nopt = N;
        popt = P;
    }

    if(p < popt)
    {
        warning(stdout);
        fprintf(stdout, "The PSF has only %" PRId64
                " slices, %" PRId64 " would be better.\n", p, popt);
        fprintf(s->log, "WARNING: The PSF has only %" PRId64 " slices, %" PRId64 " would be better.\n", p, popt);
        return psf;
    }

    if((p % 2) == 0)
    {
        if(s->verbosity > 0)
        {
            fprintf(stderr, " ! The PSF should have odd number of slices\n");
            fprintf(stderr, "   Possibly it will be auto-cropped wrong\n");
        }
        fprintf(s->log, " ! The PSF should have odd number of slices\n");
        fprintf(s->log, "   Possibly it will be auto-cropped wrong\n");
    }


    if(m > mopt || n > nopt || p > popt)
    {
        int64_t m0 = 0, m1 = m-1;
        int64_t n0 = 0, n1 = n-1;
        int64_t p0 = 0, p1 = p-1;
        if(m > mopt)
        {
            m0 = (m-mopt)/2;
            m1 = m1-(m-mopt)/2;
        }
        if(n > nopt)
        {
            n0 = (n-nopt)/2;
            n1 = n1-(n-nopt)/2;
        }
        if(p > popt)
        {
            p0 = (p-popt)/2;
            p1 = p1-(p-popt)/2;
        }
        if(s->verbosity > 2)
        {
            printf("! %" PRId64 " %" PRId64 " : %" PRId64 " %" PRId64 " : %" PRId64 " %" PRId64 "\n", m0, m1, n0, n1, p0, p1);
        }
        float * psf_cropped = fim_get_cuboid(psf, m, n, p,
                                             m0, m1, n0, n1, p0, p1);
        fim_free(psf);

        pM[0] = m1-m0+1;
        pN[0] = n1-n0+1;
        pP[0] = p1-p0+1;

        if(s->verbosity > 0)
        {
            fprintf(stdout, "PSF Z-crop [%" PRId64 " x %" PRId64 " x %" PRId64 "] -> [%" PRId64 " x %" PRId64 " x %" PRId64 "]\n",
                    m, n, p, pM[0], pN[0], pP[0]);
        }
        fprintf(s->log, "PSF Z-crop [%" PRId64 " x %" PRId64 " x %" PRId64 "] -> [%" PRId64 " x %" PRId64 " x %" PRId64 "]\n",
                m, n, p, pM[0], pN[0], pP[0]);

        return psf_cropped;
    } else {
        return psf;
    }
}

float * psf_autocrop_XY(float * psf, int64_t * pM, int64_t * pN, int64_t * pP,  // psf and size
                        __attribute__((unused))    int64_t M,
                        __attribute__((unused)) int64_t N,
                        __attribute__((unused)) int64_t P, // image size
                        dw_opts * s)
{
    int64_t m = pM[0];
    int64_t n = pN[0];
    int64_t p = pP[0];

    // Find the y-z plane with the largest sum
    float maxsum = 0;
    for(int64_t xx = 0; xx<pM[0]; xx++)
    {
        float sum = 0;
        for(int64_t yy = 0; yy<pN[0]; yy++)
        {
            for(int64_t zz = 0; zz<pP[0]; zz++)
            {
                sum += psf[xx + yy*pM[0] + zz*pM[0]*pN[0]];
            }
        }
        sum > maxsum ? maxsum = sum : 0;
    }

    //  printf("X maxsum %f\n", maxsum);

    int64_t first=-1;
    float sum = 0;

    while(sum < s->xycropfactor * maxsum)
    {
        first++;
        sum = 0;
        int64_t xx = first;
        for(int64_t yy = 0; yy<pN[0]; yy++)
        {
            for(int64_t zz = 0; zz<pP[0]; zz++)
            {
                sum += psf[xx + yy*pM[0] + zz*pM[0]*pN[0]];
            }
        }
    }

    if(first < 1)
    {
        if(s->verbosity > 1)
        {
            printf("PSF X-crop: Not cropping\n");
        }
        return psf;
    }

    // Benchmark FFTW to figure out a good compromise between cropping
    // and a size that will be fast for fftw.
    if(s->lookahead > 0)
    {
        fprintf(s->log, "> Testing lookahead up to %d\n", s->lookahead);
        printf("Suggested PSF size: %" PRId64 "\n", pM[0] - 2*first);
        int64_t imsize = M + (pM[0]-1) - 2*first;
        int64_t imsize_max = M + (pM[0]-1);
        printf("Gives images size: %" PRId64 " -- %" PRId64 "\n", imsize, imsize_max);
        if(imsize_max - imsize > s->lookahead)
        {
            imsize_max = imsize+s->lookahead;
        }
        printf("lookahead gives: images size: %" PRId64 " -- %" PRId64 "\n", imsize, imsize_max);

        double * times = fft_bench_1d(imsize, imsize_max, 100);
        int bestAdd = 0;
        double bestTime = INFINITY;
        for(int kk = imsize; kk<= imsize_max; kk+=2)
        {
            if(times[kk-imsize] < bestTime)
            {
                bestTime = times[kk-imsize];
                bestAdd = (kk-imsize)/2;
            }
            printf("job x-size: %d: %f\n", kk, times[kk-imsize]);
        }
        fflush(stdout);
        double gain = times[0]/bestTime;
        fprintf(s->log, "Lookahead gain: %f\n", gain);
        if(s->verbosity > 1)
        {
            printf("Lookahead gain: %f\n", gain);
        }
        free(times);
        first = first-bestAdd;
    }

    /* Allocated with fftwf */
    float * crop = fim_get_cuboid(psf, pM[0], pN[0], pP[0],
                                  first, pM[0] - first -1,
                                  first, pN[0] - first -1,
                                  0, pP[0]-1);
    pM[0] -= 2*first;
    pN[0] -= 2*first;

    if(s->verbosity > 0)
    {
        fprintf(stdout, "PSF XY-crop [%" PRId64 " x %" PRId64 " x %" PRId64 "] -> [%" PRId64 " x %" PRId64 " x %" PRId64 "]\n",
                m, n, p, pM[0], pN[0], pP[0]);
    }
    fprintf(s->log, "PSF XY-crop [%" PRId64 " x %" PRId64 " x %" PRId64 "] -> [%" PRId64 " x %" PRId64 " x %" PRId64 "]\n",
            m, n, p, pM[0], pN[0], pP[0]);

    fim_free(psf);
    return crop;
}

float * psf_autocrop(float * psf, int64_t * pM, int64_t * pN, int64_t * pP,  // psf and size
                     int64_t M, int64_t N, int64_t P, // image size
                     dw_opts * s)
{
    float * p = psf;
    assert(pM[0] > 0);

    /* Crop the PSF if it is larger than necessary */
    p = psf_autocrop_byImage(p, pM, pN, pP,
                             M, N, P,
                             s);

    /* Crop the PSF by removing outer planes that has very little information.
     * Only if the PSF is larger than the image in some dimension. */
    if((s->borderQuality > 0))
    {
        if(s->xycropfactor > 0)
        {
            p = psf_autocrop_XY(p, pM, pN, pP, M, N, P, s);
        }
    }
    assert(pM[0] > 0);
    assert(p != NULL);
    return p;
}


static int
deconvolve_tiles(const int64_t M, const int64_t N, const int64_t P,
                 const float * restrict psf, const int64_t pM, const int64_t pN, const int64_t pP,
                 dw_opts * s)
{

    tiling * T = tiling_create(M,N,P, s->tiling_maxSize, s->tiling_padding);
    if( T == NULL)
    {
        fprintf(stderr, "Tiling failed, please check your settings\n");
        exit(1);
    }

    if( T->nTiles == 1)
    {
        fprintf(stderr, "\n"
                "ERROR: Only one tile! Please omit the `--tilesize` parameter if "
                "that is what you intended to to or decrease the value if you "
                "want to process the image in tiles."
                "\n\n");
        exit(1);
    }

    if(s->verbosity > 0)
    {
        printf("-> Divided the [%" PRId64 " x %" PRId64 " x %" PRId64 "] image into %d tiles\n", M, N, P, T->nTiles);
    }

    /* Output image initialize as zeros
     * will be updated block by block
     */
    char * tfile = malloc(strlen(s->outFile)+10);
    assert(tfile != NULL);
    sprintf(tfile, "%s.raw", s->outFile);

    if(s->verbosity > 0)
    {
        printf("Initializing %s to 0\n", tfile); fflush(stdout);
    }
    fsetzeros(tfile, (size_t) M* (size_t) N* (size_t) P*sizeof(float));

    char * imFileRaw = malloc(strlen(s->imFile) + 10);
    assert(imFileRaw != NULL);
    sprintf(imFileRaw, "%s.raw", s->imFile);

    if(s->verbosity > 0)
    {
        printf("Dumping %s to %s (for quicker io)\n", s->imFile, imFileRaw);
    }

    fim_tiff_to_raw(s->imFile, imFileRaw);
    if(s->verbosity > 10){
        printf("Writing to imdump.tif\n");
        fim_tiff_from_raw("imdump.tif", M, N, P, imFileRaw, NULL);
    }

    //fim_tiff_write_zeros(s->outFile, M, N, P);
    if(s->verbosity > 0)
    {
        printf("\n"); fflush(stdout);
    }

    int nTiles = T->nTiles;
    if(s->onetile == 1)
    {
        nTiles = 1;
        fprintf(s->log, "DEBUG: only the first tile to be deconvolved\n");
        fprintf(stdout, "DEBUG: only the first tile to be deconvolved\n");
    }

    for(int tt = 0; tt < nTiles; tt++)
    {
        // Temporal copy of the PSF that might be cropped to fit the tile
        float * tpsf = fim_copy(psf, pM*pN*pP);
        int64_t tpM = pM, tpN = pN, tpP = pP;

        if(s->verbosity > 0)
        {
            printf("-> Processing tile %d / %d\n", tt+1, T->nTiles);
            fprintf(s->log, "-> Processing tile %d / %d\n", tt+1, T->nTiles);
        }

        //    tictoc
        //   tic
        //float * im_tile = tiling_get_tile_tiff(T, tt, s->imFile);
        float * im_tile = tiling_get_tile_raw(T, tt, imFileRaw);
        //    toc(tiling_get_tile_tiff)

        int64_t tileM = T->tiles[tt]->xsize[0];
        int64_t tileN = T->tiles[tt]->xsize[1];
        int64_t tileP = T->tiles[tt]->xsize[2];

        if(0)
        {
            printf("writing to tiledump.tif\n");
            fim_tiff_write("tiledump.tif", im_tile, NULL, tileM, tileN, tileP);
            getchar();
        }

        fim_normalize_sum1(tpsf, tpM, tpN, tpP);

        tpsf = psf_autocrop(tpsf, &tpM, &tpN, &tpP,
                            tileM, tileN, tileP, s);

        fim_normalize_sum1(tpsf, tpM, tpN, tpP);

        if(s->offset > 0)
        {
            fim_add_scalar(im_tile, tileM*tileN*tileP, s->offset);
        }

        float * dw_im_tile = s->fun(im_tile, tileM, tileN, tileP, // input image and size
                                    tpsf, tpM, tpN, tpP, // psf and size
                                    s);
        fim_free(im_tile);
        if(s->offset > 0)
        {
            fim_add_scalar(dw_im_tile, tileM*tileN*tileP, -s->offset);
            fim_project_positive(dw_im_tile, tileM*tileN*tileP);
        }
        if(s->verbosity > 1)
        {
            printf("Saving tile to disk\n");
        }
        tiling_put_tile_raw(T, tt, tfile, dw_im_tile);
        fim_free(dw_im_tile);
        // free(tpsf);
    }
    tiling_free(T);
    free(T);

    if(s->verbosity > 2)
    {
        printf("converting %s to %s\n", tfile, s->outFile);
    }

    fim_tiff_from_raw(s->outFile, M, N, P, tfile, s->imFile);

    if(s->verbosity > 1)
    {
        printf("conversion done\n");
    }

    if(s->verbosity < 5)
    {
        remove(tfile);
    } else {
        printf("Keeping %s for inspection, remove manually\n", tfile);
    }

    if(s->verbosity > 2)
    {
        printf("freeing up\n");
    }
    free(tfile);
    remove(imFileRaw);
    free(imFileRaw);
    if(s->verbosity > 2)
    {
        printf("Done with tiling\n");
    }
    return 0;
}



void timings()
{
    printf("-> Timings\n");
    tictoc
        int64_t M = 1024, N = 1024, P = 50;

#ifndef WINDOWS
    tic
        usleep(1000);
    toc(usleep_1000)
#endif
        tic
        float * V = fim_malloc(M*N*P*sizeof(float));
    toc(malloc)

        float * A = fim_malloc(M*N*P*sizeof(float));

    tic
        memset(V, 0, M*N*P*sizeof(float));
    toc(memset)
        memset(A, 0, M*N*P*sizeof(float));


    tic
        for(size_t kk = 0; kk < (size_t) M*N*P; kk++)
        {
            A[kk] = (float) rand()/(float) RAND_MAX;
        }
    toc(rand)

        // ---
        tic
        fftwf_plan p = fftwf_plan_dft_r2c_3d(P, N, M,
                                             V, NULL,
                                             FFTW_WISDOM_ONLY | FFTW_MEASURE);
    fftwf_destroy_plan(p);
    toc(fftwf_plan_create_and_destroy)


        // ---
        tic
        fim_flipall(V, A, M, N, P);
    toc(fim_flipall)

        // ---
        tic
        float e1 = getError_ref(V, A, M, N, P, M, N, P);
    toc(getError_ref)
        V[0]+= e1;

    tic
        float e2 = getError(V, A, M, N, P, M, N, P, DW_METRIC_MSE);
    toc(getError)
        V[0]+=e2;

    printf("e1: %f, e2: %f, fabs(e1-e2): %f\n", e1, e1, fabs(e1-e2));

    // ---
    tic
        float * S1 = fim_subregion(V, M, N, P, M-1, N-1, P-1);
    toc(fim_subregion)

        tic
        float * S2 = fim_subregion_ref(V, M, N, P, M-1, N-1, P-1);
    toc(fim_subregion_ref)
        printf("S1 - S2 = %f\n", getError(S1, S1, M-1, N-1, P-1, M-1, N-1, P-1, DW_METRIC_MSE));
    free(S1);
    free(S2);

    // ---
    tic
        fim_insert(V, M, N, P, A, M-1, N-1, P-1);
    toc(fim_subregion)

        tic
        fim_insert_ref(V, M, N, P, A, M-1, N-1, P-1);
    toc(fim_subregion_ref)

        // ---


        ((float volatile *)V)[0] = V[0];
    printf("V[0] = %f\n", V[0]);
    free(A);
    free(V);
}

void dw_unittests()
{
    fprint_peak_memory(stdout);
    timings();

    //fim_ut();
    fim_tiff_ut();
    fft_ut();
    printf("done\n");
}

void show_time(FILE * f)
{
    f == NULL ? f = stdout : 0;
    int buf_len = 256;
    char buf[256];

    time_t now = time(NULL);
    struct tm * ptm = localtime(&now);
    strftime(buf, buf_len, "%FT%T", ptm);

    fprintf(f, "%s\n", buf);
    return;
}

float * psf_makeOdd(float * psf, int64_t * pM, int64_t * pN, int64_t *pP)
{
    // Expand the psf so that it had odd dimensions it if doesn't already have that
    int64_t m = pM[0];
    int64_t n = pN[0];
    int64_t p = pP[0];
    int reshape = 0;
    if(m % 2 == 0)
    { m++; reshape = 1;}
    if(n % 2 == 0)
    { n++; reshape = 1;}
    if(p % 2 == 0)
    { p++; reshape = 1;}

    if(reshape == 0)
    {  return psf; }
    // printf("%d %d %d -> %d %d %d\n", pM[0], pN[0], pP[0], m, n, p);
    float * psf2 = fim_zeros(m*n*p);
    fim_insert(psf2, m, n, p, psf, pM[0], pN[0], pP[0]);
    free(psf);
    pM[0] = m;
    pN[0] = n;
    pP[0] = p;
    return psf2;
}


void dcw_init_log(dw_opts * s)
{
    s->log = fopen(s->logFile, "w");
    if (s->log == NULL)
    {
        fprintf(stderr, "Unable to open %s for writing\n", s->logFile);
        fprintf(stderr,
                "Please check that you have permissions to write to the folder\n"
                "and that the drive is not full\n");
        exit(EXIT_FAILURE);
    }
    assert(s->log != NULL);
    show_time(s->log);
    dw_opts_fprint(s->log, s);
    dw_fprint_info(s->log, s);
}

void dcw_close_log(dw_opts * s)
{
    fprint_peak_memory(s->log);
    show_time(s->log);
    fclose(s->log);
}

double get_nbg(float * I, size_t N, float bg)
{ // count the number of pixels in the image
    // that are set to the background level
    double nbg = 0;
    for(size_t kk = 0; kk<N; kk++)
    {
        if(I[kk] == bg)
        {
            nbg++;
        }
    }
    return nbg;
}


void flatfieldCorrection(dw_opts * s, float * im, int64_t M, int64_t N, int64_t P)
{
    printf("Experimental: applying flat field correction using %s\n",
           s->flatfieldFile);
    ttags * T = ttags_new();
    int64_t m = 0, n = 0, p = 0;
    float * C = fim_tiff_read(s->flatfieldFile, T, &m, &n, &p, s->verbosity);
    ttags_free(&T);

    assert(m == M);
    assert(n == N);
    assert(p == 1);

    // TODO:
    // check that it is positive

    for(int64_t zz = 0; zz<P; zz++)
    {
        for(size_t pos = 0; pos < (size_t) M*N; pos++)
        {
            im[M*N*zz + pos] /= C[pos];
        }
    }
    free(C);
}


static void prefilter(dw_opts * s,
                      float * im, int64_t M, int64_t N, int64_t P,
                      __attribute__((unused)) float * psf,
                      __attribute__((unused)) int64_t pM,
                      __attribute__((unused)) int64_t pN,
                      __attribute__((unused)) int64_t pP)
{

    if(s->psigma<= 0)
    {
        return;
    }

    fim_anscombe(im, M*N*P);
    fim_gsmooth(im, M, N, P, s->psigma);
    fim_ianscombe(im, M*N*P);

    return;
}

/** Create initial guess: the fft of an image that is 1 in MNP and 0
 * outside M, N, P is the dimension of the microscopic image
 *
 * Possibly more stable to use the mean of the input image rather than 1
 */
fftwf_complex * initial_guess(const int64_t M, const int64_t N, const int64_t P,
                              const int64_t wM, const int64_t wN, const int64_t wP)
{
    assert(wM >= M); assert(wN >= N); assert(wP >= P);

    float * one = fim_zeros(wM*wN*wP);

#pragma omp parallel for shared(one)
    for(int64_t cc = 0; cc < P; cc++) {
        for(int64_t bb = 0; bb < N; bb++) {
            for(int64_t aa = 0; aa < M; aa++) {
                one[aa + wM*bb + wM*wN*cc] = 1;
            }
        }
    }
    //printf("writing to one.tif");
    //fim_tiff_write_float("one.tif", one, NULL, wM, wN, wP);
    //  writetif("one.tif", one, wM, wN, wP);

    fftwf_complex * Fone = fft(one, wM, wN, wP);

    fim_free(one);
    return Fone;
}

static void dw_set_omp_threads(const dw_opts *s)
{
#ifdef _OPENMP
#ifdef MKL
    mkl_set_num_threads(s->nThreads_FFT);
    fprintf(s->log, "Set the number of MKL threads to %d\n", s->nThreads_FFT);
#else
    omp_set_num_threads(s->nThreads_OMP);
    fprintf(s->log, "Set the number of OMP threads to %d\n", s->nThreads_OMP);
    /* Fastest of static, dynamic and guided in limited tests */
    omp_set_dynamic(false);
    omp_set_schedule(omp_sched_static, 0);
    fprintf(s->log, "Using static scheduling for OMP\n");
    /* Better on heavily loaded machine? */
    //omp_set_schedule(omp_sched_guided, 0);
#endif
#endif
    return;
}

int dw_run(dw_opts * s)
{
    struct timespec tstart, tend;
    dw_gettime(&tstart);
    dcw_init_log(s);

    if(s->verbosity > 1)
    {
        dw_opts_fprint(NULL, s);
        printf("\n");
    }

    s->verbosity > 1 ? dw_fprint_info(NULL, s) : 0;
    dw_set_omp_threads(s);

    logfile = stdout;

    fim_tiff_init();
    fim_tiff_set_log(s->log);
    if(s->verbosity > 2)
    {
        fim_tiff_set_log(s->log);
    }

    int64_t M = 0, N = 0, P = 0;
    if(fim_tiff_get_size(s->imFile, &M, &N, &P))
    {
        printf("Failed to open %s\n", s->imFile);
        return -1;
    } else {
        if(s->verbosity > 3)
        {
            printf("Got image info from %s\n", s->imFile);
        }
    }

    if(s->verbosity > 1)
    {
        printf("Image dimensions: %" PRId64 " x %" PRId64 " x %" PRId64 "\n",
               M, N, P);
    }


    int tiling = 0;
    if(s->tiling_maxSize > 0 && (M > s->tiling_maxSize || N > s->tiling_maxSize))
    {
        tiling = 1;
    }

    float * im = NULL;
    ttags * T = ttags_new();

    if(tiling == 0)
    {
        if(s->verbosity > 0 )
        {
            printf("Reading %s\n", s->imFile);
        }

        im = fim_tiff_read(s->imFile, T, &M, &N, &P, s->verbosity);
        if(im == NULL)
        {
            fprintf(stderr, "Failed to open %s\n", s->imFile);
            exit(EXIT_FAILURE);
        }

        if(s->verbosity > 4)
        {
            if(M > 9)
            {
                printf("image data: ");
                for(size_t kk = 0; kk<10; kk++)
                {
                    printf("%f ", im[kk]);
                }
                printf("\n");
            }
            printf("Done reading\n"); fflush(stdout);
        }

        if(s->auto_zcrop > 0)
        {
            if(s->verbosity > 0)
            {
                printf("Cropping the image to %" PRId64 " x %" PRId64 " x %d\n",
                       M, N, s->auto_zcrop);
            }
            float * zim = fim_auto_zcrop(im, M, N, P, s->auto_zcrop);
            if(zim == NULL)
            {
                fprintf(stderr,
                        "Automatic cropping failed\n");
                exit(EXIT_FAILURE);
            }
            fim_free(im);
            im = zim;
            P = s->auto_zcrop;

        }

        if(s->zcrop > 0)
        {
            if(2*s->zcrop > P)
            {
                fprintf(stderr, "Impossible to remove 2x%d planes from an image with %zu planes\n",
                        s->zcrop, P);
            }
            if(s->verbosity > 0)
            {
                printf("Removing %d planes from the top and bottom of the image\n",
                       s->zcrop);
            }
            float * zim = fim_zcrop(im, M, N, P, (size_t )s->zcrop);
            if(zim == NULL)
            {
                fprintf(stderr,
                        "Automatic cropping failed\n");
                exit(EXIT_FAILURE);
            }
            fim_free(im);
            im = zim;
            P = P - 2*s->zcrop;
            if(s->verbosity > 0)
            {
                printf("New image size: [%" PRId64 "x %" PRId64 "x %" PRId64 "]\n",
                       M, N, P);
            }
        }


        if(fim_min(im, M*N*P) < 0)
        {
            fprintf(stderr,
                    "ERROR: The image contains negative values, can not continue!\n");
            fprintf(s->log,
                    "ERROR: The image contains negative values, can not continue!\n");
            exit(EXIT_FAILURE);
        }
        float maxval = fim_max(im, M*N*P);
        if(maxval < 1.0)
        {
            fprintf(stderr,
                    "ERROR: The image has too low intensity, can not continue!\n");
            fprintf(s->log,
                    "ERROR: The image has too low intensity, can not continue!\n");
            fprintf(s->log, "The largest value of the input image is %f\n",
                    maxval);
            exit(EXIT_FAILURE);
        }

        if(maxval < 100)
        {
            if(s->verbosity > 0)
            {
                printf("WARNING: The largest value of the input image is %f\n",
                       maxval);
            }
            fprintf(s->log,
                    "WARNING: The largest value of the input image is %f\n",
                    maxval);
        }


        if(s->refFile != NULL)
        {
            int64_t rM = 0, rN = 0, rP = 0;
            s->ref = fim_tiff_read(s->refFile, NULL, &rM, &rN, &rP, s->verbosity);
            if(s->ref == NULL)
            {
                fprintf(stderr, "Failed to open %s\n", s->imFile);
                exit(1);
            }
            if( (rM != M) || (rN != N) || (rP != P) )
            {
                fprintf(stderr, "Image and reference image does not have matching size\n");
                exit(1);
            }
        }

    }

    // Set up the string for the TIFFTAG_SOFTWARE
    char * swstring = malloc(1024);
    assert(swstring != NULL);
    sprintf(swstring, "deconwolf %s", deconwolf_version);
    ttags_set_software(T, swstring);
    free(swstring);

    // fim_tiff_write("identity.tif", im, M, N, P);

    int64_t pM = 0, pN = 0, pP = 0;
    float * psf = NULL;

    if(s->verbosity > 0)
    {
        printf("Reading %s\n", s->psfFile);
    }
    psf = fim_tiff_read(s->psfFile, NULL, &pM, &pN, &pP, s->verbosity);
    if(psf == NULL)
    {
        fprintf(stderr, "Failed to open %s\n", s->psfFile);
        exit(1);
    }
    if(s->verbosity > 4)
    {
        if(M > 9)
        {
            printf("image data: ");
            for(size_t kk = 0; kk<10; kk++)
            {
                printf("%f ", psf[kk]);
            }
            printf("\n");
        }
    }

    if(fim_maxAtOrigo(psf, pM, pN, pP) == 0)
    {
        /* It might still be centered between pixels */
        if(s->verbosity > 0)
        {
            warning(stdout);
            printf("The PSF is not centered!\n");
        }
        fprintf(s->log, " ! The PSF is not centered\n");
    }

    // Possibly the PSF will be cropped even more per tile later on

    fim_normalize_sum1(psf, pM, pN, pP);
    if(1)
    {
        psf = psf_autocrop(psf, &pM, &pN, &pP,
                           M, N, P, s);
    }


    if(s->verbosity > 0)
    {
        printf("Output: %s(.log.txt)\n", s->outFile);
    }

    myfftw_start(s->nThreads_FFT, s->verbosity, s->log);

    float * out = NULL;

    if(tiling)
    {
        if(s->flatfieldFile != NULL)
        {
            warning(stdout);
            printf("Flat-field correction can't be used in tiled mode\n");
        }
        deconvolve_tiles(M, N, P,
                         psf, pM, pN, pP, // psf and size
                         s);// settings
        fim_free(psf);
    } else {
        fim_normalize_sum1(psf, pM, pN, pP);
        if(s->flatfieldFile != NULL)
        {
            flatfieldCorrection(s, im, M, N, P);
        }

        /* Pre filter by psigma */
        prefilter(s, im, M, N, P, psf, pM, pN, pP);

        if(s->offset > 0)
        {
            fim_add_scalar(im, M*N*P, s->offset);
        }

        /* Note: psf is freed bu the deconvolve_* functions*/
        out = s->fun(im, M, N, P, // input image and size
                     psf, pM, pN, pP, // psf and size
                     s);// settings

        if(s->offset > 0)
        {
            fim_add_scalar(im, M*N*P, -s->offset);
            fim_project_positive(im, M*N*P);
        }

        psf = NULL;
    }

    if(tiling == 0)
    {
        fim_free(im);


        if(out == NULL)
        {
            if(s->verbosity > 0)
            {
                printf("Nothing to write to disk :(\n");
            }
        } else
        {
            double nZeros = get_nbg(out, M*N*P, s->bg);
            fprintf(s->log, "%f%% pixels at bg level in the output image.\n", 100*nZeros/(M*N*P));
            if(s->verbosity > 1)
            {
                printf("%f%% pixels at bg level in the output image.\n", 100*nZeros/(M*N*P));
                printf("Writing to %s\n", s->outFile); fflush(stdout);
            }

            //    floatimage_normalize(out, M*N*P);
            if(s->outFormat == 32)
            {
                if(s->iterdump)
                {
                    char * outFile = gen_iterdump_name(s, s->nIter);
                    fim_tiff_write_float(outFile, out, T, M, N, P);
                    free(outFile);
                } else {
                    fim_tiff_write_float(s->outFile, out, T, M, N, P);
                }
            } else {
                if(s->iterdump)
                {
                    char * outFile = gen_iterdump_name(s, s->nIter);
                    fim_tiff_write(outFile, out, T, M, N, P);
                    free(outFile);
                } else {
                    fim_tiff_write_opt(s->outFile, out, T, M, N, P, s->scaling);
                }
            }
        }
    }

    ttags_free(&T);

    if(s->verbosity > 1)
    {
        printf("Finalizing "); fflush(stdout);
    }

    fim_free(out);
    myfftw_stop();

    dw_gettime(&tend);
    fprintf(s->log, "Took: %f s\n", timespec_diff(&tend, &tstart));
    dcw_close_log(s);

    if(s->verbosity > 1)
    { fprint_peak_memory(stdout); }

    if(s->verbosity > 0)
    { printf("Done!\n"); }

    dw_opts_free(&s);

    return 0;
}
