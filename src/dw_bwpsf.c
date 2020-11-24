#include "dw_bwpsf.h"
#include "li.c"

pthread_mutex_t stdout_mutex;

void bw_conf_printf(FILE * out, bw_conf * conf)
{
    fprintf(out, "cmd: %s\n", conf->cmd);
    fprintf(out, "lambda = %.2f nm\n", conf->lambda*1e9);
    fprintf(out, "NA = %f\n", conf->NA);
    fprintf(out, "ni = %f\n", conf->ni);
    fprintf(out, "TOL = %f\n", conf->TOL);
    fprintf(out, "resLateral = %.2f nm\n", conf->resLateral*1e9);
    fprintf(out, "resAxial = %.2f nm\n", conf->resAxial*1e9);
    fprintf(out, "out size: [%d x %d x %d] pixels\n", conf->M, conf->N, conf->P);
    fprintf(out, "nThreads: %d\n", conf->nThreads);
    fprintf(out, "Verbosity: %d\n", conf->verbose);
    fprintf(out, "Overwrite: %d\n", conf->overwrite);
    fprintf(out, "File: %s\n", conf->outFile);
    fprintf(out, "Log: %s\n", conf->logFile);
    if(conf->Simpson)
    {
        fprintf(out, "Pixel integration: %dx%dx%d samples per pixel\n", conf->Simpson_N, conf->Simpson_N, conf->Simpson_N);
        fprintf(out, "Oversampling of radial profile: %d X\n", conf->oversampling_R);
    }
}


bw_conf * bw_conf_new()
{
    bw_conf * conf = malloc(sizeof(bw_conf));
    conf->lambda = 600*1e-9;
    conf->NA = 1.45;
    conf->ni = 1.515;
    conf->TOL = 1e-1;
    conf->K = 9; // corresponding to "best" in PSFGenerator
    conf->M = 181;
    conf->N = 181;
    conf->P = 181;
    conf->resAxial = 300*1e-9;
    conf->resLateral = 130*1e-9;
    conf->nThreads = 4;
    conf->V = NULL;
    conf->verbose = 1;
    conf->V = NULL;
    conf->outFile = NULL;
    conf->logFile = NULL;
    conf->overwrite = 0;
    conf->Simpson = 2;
    conf->Simpson_N = 5;
    conf->oversampling_R = 10;
    return conf;
}

static double timespec_diff(struct timespec* end, struct timespec * start)
{
    double elapsed = (end->tv_sec - start->tv_sec);
    elapsed += (end->tv_nsec - start->tv_nsec) / 1000000000.0;
    return elapsed;
}

void getCmdLine(int argc, char ** argv, bw_conf * s)
{
    // Copy the command line to s->cmd
    int lcmd=0;
    for(int kk = 0; kk<argc; kk++)
    {
        lcmd += strlen(argv[kk]);
    }
    lcmd += argc+2;
    s->cmd = malloc(lcmd);
    int pos = 0;
    for(int kk = 0; kk<argc; kk++)
    {
        sprintf(s->cmd+pos, "%s ", argv[kk]);
        pos += strlen(argv[kk])+1;
    }
    s->cmd[pos-1] = '\0';
}


int file_exist(char * fname)
{
    if( access( fname, F_OK ) != -1 ) {
        return 1; // File exist
    } else {
        return 0;
    }
}

void fprint_time(FILE * f)
{
    f == NULL ? f = stdout : 0;
    time_t now = time(NULL);
    char * tstring = ctime(&now);
    fprintf(f, "%s\n", tstring);
}

void bw_version(FILE* f)
{
    fprintf(f, "deconwolf: '%s' PID: %d\n", deconwolf_version, (int) getpid());
}

void usage(__attribute__((unused)) int argc, char ** argv)
{
    printf("Usage:\n");
    printf("\t$ %s <options> psf.tif \n", argv[0]);
    printf("\n");
    printf(" Options:\n");
    printf(" --version\n\t Show version info\n");
    printf(" --help\n\t Show this message\n");
    printf(" --verbose L\n\t Set verbosity level to L\n");
    printf(" --NA\n\t Set numerical aperture\n");
    printf(" --lamda\n\t Set emission wavelength [nm]\n");
    printf(" --ni\n\t Set refractive index\n");
    printf(" --threads\n\t Set number of threads\n");
    printf(" --resxy\n\t Set pixel size in x-y [nm]\n");
    printf(" --resz\n\t Set pixel size in z [nm]\n");
    printf(" --size N \n\t Set output size to N x N x N [pixels]\n");
    printf(" --overwrite \n\t Toggles overwriting of existing files to YES\n");
    printf(" --quality N\n\t Sets the integration quality to NxNxN samples per pixel\n");
    printf("\b");
    return;
}

void bw_argparsing(int argc, char ** argv, bw_conf * s)
{

    if(argc < 2)
    {
        // Well it could run with the defaults but that does not make sense
        usage(argc, argv);
        exit(-1);
    }

    getCmdLine(argc, argv, s);

    struct option longopts[] = {
                                { "version",     no_argument,       NULL,   'v' },
                                { "help",         no_argument,       NULL,   'h' },
                                // Settings
                                { "lambda", required_argument, NULL, 'l'},
                                { "threads",      required_argument, NULL,   't' },
                                { "verbose",      required_argument, NULL,   'p' },
                                { "overwrite",   no_argument,        NULL,   'w' },
                                { "NA", required_argument, NULL, 'n' },
                                { "ni", required_argument, NULL, 'i' },
                                { "resxy", required_argument, NULL, 'x'},
                                { "resz", required_argument, NULL, 'z'},
                                { "size", required_argument, NULL, 'N'},
                                { "quality", required_argument, NULL, 'q'},
                                { NULL,           0,                 NULL,   0   }
    };
    int ch;
    while((ch = getopt_long(argc, argv, "q:vho:t:p:w:n:i:x:z:N:l:", longopts, NULL)) != -1)
    {
        switch(ch) {
        case 'q':
            s->Simpson_N = atoi(optarg);
            if(s->Simpson_N == 0)
            {
                s->Simpson = 0;
            }
            break;
        case 'v':
            bw_version(stdout);
            exit(0);
            break;
        case 'h':
            usage(argc, argv);
            exit(0);
        case 'o':
            s->outFile = malloc(strlen(optarg) + 1);
            sprintf(s->outFile, "%s", optarg);
            break;
        case 't':
            s->nThreads = atoi(optarg);
            break;
        case 'p':
            s->verbose = atoi(optarg);
            break;
        case 'w':
            s->overwrite = 1;
            break;
        case 'n':
            s->NA = atof(optarg);
            break;
        case 'i':
            s->ni = atof(optarg);
            break;
        case 'x':
            s->resLateral = atof(optarg)*1e-9;
            break;
        case 'z':
            s->resAxial = atof(optarg)*1e-9;
            break;
        case 'l':
            s->lambda = atof(optarg)*1e-9;
            break;
        case 'N':
            s->M = atoi(optarg);
            s->N = s->M;
            s->P = s->M;
            break;
        }
    }

    if(s->lambda < 50*1e-9)
    {
        printf("Error: lambda has be be at least 50 nm\n");
        exit(-1);
    }

    if(s->M % 2 == 0)
    {
        printf("Error: The size has to be odd, 1, 3, ...\n");
        exit(-1);
    }

    if(optind + 1 == argc)
    {
        s->outFile = malloc(strlen(argv[argc-1]) + 1);
        sprintf(s->outFile, "%s", argv[argc-1]);
    }

    if(s->outFile == NULL)
    {
        s->outFile = malloc(100*sizeof(char));
        sprintf(s->outFile, "PSFBW_%.2f_%.2f_%.1f_%.1f_%.1f.tif",
                s->NA, s->ni, s->lambda*1e9, s->resLateral*1e9, s->resAxial*1e9);
    }

    s->logFile = malloc(strlen(s->outFile) + 10);
    sprintf(s->logFile, "%s.log.txt", s->outFile);
}

void * BW_thread(void * data)
{
    // Entry point for pthread_create
    bw_conf * conf = (bw_conf *) data;

    if(conf->verbose > 3)
    {
        pthread_mutex_lock(&stdout_mutex);
        printf("-> From thread %d\n", conf->thread);
        bw_conf_printf(stdout, conf);
        pthread_mutex_unlock(&stdout_mutex);
    }

    float * V = conf->V;
    int M = conf->M;
    int N = conf->N;
    int P = conf->P;
    size_t MN = M*N;

    for (int z = conf->thread; z <= (P-1)/2; z+=conf->nThreads) {
        float defocus = conf->resAxial * (z - (P - 1.0) / 2.0);
        BW_slice(V + z*MN, defocus, conf);
    }
    return NULL;
}

void BW(bw_conf * conf)
{
    assert(conf->V == NULL);
    assert(conf->nThreads > 0);
    conf->V = malloc(conf->M*conf->N*conf->P*sizeof(float));

    float * V = conf->V;
    int M = conf->M;
    int N = conf->N;
    int P = conf->P;
    size_t MN = M*N;

    int nThreads = conf->nThreads;

    if(nThreads == 1)
    {
        for (int z = 0; z <= (P-1)/2; z++) {
            float defocus = conf->resAxial * (z - (P - 1.0) / 2.0);
            BW_slice(conf->V + z*M*N, defocus, conf);
        }
    } else {
        pthread_t * threads = malloc(nThreads*sizeof(pthread_t));
        bw_conf ** confs = malloc(nThreads*sizeof(bw_conf*));

        for(int kk = 0; kk<nThreads; kk++)
        {
            confs[kk] = (bw_conf*) malloc(sizeof(bw_conf));
            memcpy(confs[kk], conf, sizeof(bw_conf));
            confs[kk]->thread = kk;
            pthread_create(&threads[kk], NULL, BW_thread, (void *) confs[kk]);
        }

        for(int kk = 0; kk<nThreads; kk++)
        {
            pthread_join(threads[kk], NULL);
            free(confs[kk]);
        }
        free(confs);
        free(threads);
    }

    // symmetry in Z
    for(int z = 0; z<(P-1)/2; z++)
    {
        memcpy(V+(P-z-1)*MN, V+z*MN, MN*sizeof(float));
    }
}


float complex integrand(float rho, float r, float defocus, bw_conf * conf) {

    // 'rho' is the integration parameter.
    // 'r' is the radial distance of the detector relative to the optical
    // axis.
    // NA is assumed to be less than 1.0, i.e. it assumed to be already
    // normalized by the refractive index of the immersion layer, ni.
    // The return value is a complex number.

    assert(rho<=1); assert(rho>=0);

    float k0 = 2.0 * M_PI / conf->lambda;
    // j0 has smaller errors but I don't think that we need that precision.
    // j0f is much faster

    float BesselValue = j0f(k0 * conf->NA * r * rho);

    // Optical path difference
    float OPD = pow(conf->NA,2) * defocus * pow(rho,2) / (2.0 * conf->ni);
    // Phase aberrations.
    float W = k0 * OPD;

    float complex y = BesselValue*cos(W)*rho - I*BesselValue*sin(W)*rho;
    return y;
}



// Simpson approximation for the Kirchhoff diffraction integral
// 'r' is the radial distance of the detector relative to the optical axis.
float calculate(float r, float defocus, bw_conf * conf) {
    // is the stopping criterion really good?
    // why del^2 -- isn't it a 1D integral?
    // Better set a tolerance relative to the wanted precision (never more than 32-bit)
    // doubling the number of points 9 times!!

    float curDifference = conf->TOL; // Stopping criterion

    complex float value = 0;

    float curI = 0.0, prevI = 0.0;

    // Initialization of the Simpson sum (first iteration)
    int64_t N = 2; // number of sub-intervals
    double del = 1.0/N; // sub-interval length


    // Initial 3-point estimation
    complex float sumOddIndex = integrand(0.5, r, defocus, conf);
    complex float sumEvenIndex = 0;
    complex float valueX0 = integrand(0.0, r, defocus, conf);
    complex float valueXn = integrand(1.0, r, defocus, conf);
    float complex sum = valueX0 + 2*sumEvenIndex + 4*sumOddIndex + valueXn;
    curI = (pow(creal(sum),2) + pow(cimag(sum), 2)) * pow(del,2);

    prevI = curI;

    // Increase the number of points until desired precision

    size_t k = 0; // number of consecutive successful approximations
    size_t iteration = 1;
    size_t max_iterations = 16;
    while (k < conf->K && iteration <= max_iterations) {
        //    printf("%zu\n", iteration); fflush(stdout);
        iteration++;
        N *= 2;
        del /= 2;

        sumEvenIndex = sumEvenIndex + sumOddIndex;
        sumOddIndex = 0 + I*0;

        for (int64_t n = 1; n < N; n = n + 2) {
            double rho = n * del;
            value = integrand(rho, r, defocus, conf);
            sumOddIndex += value;
        }

        complex float sum = valueX0 + 2*sumEvenIndex + 4*sumOddIndex + valueXn;
        curI = (pow(creal(sum),2) + pow(cimag(sum), 2)) * pow(del, 2);

        // Relative error between consecutive approximations
        if (prevI == 0.0)
        {
            curDifference = fabs((prevI - curI) / 1E-5);
        } else {
            curDifference = fabs((prevI - curI) / curI);
        }
        if (curDifference <= conf->TOL)
        {
            k++;
        } else {
            k = 0;
        }

        if(fabs(curI-prevI) < 1e-12)
        {
            // This gives a quick finish for close-to-zero pixels
            //   printf("Absolute error break\n");
            //   printf("Iteration: %zu, curDifference: %e, curI: %e\n", iteration, curDifference, curI);
            break;
        }
        prevI = curI;
    }

    return curI/9;
}

double simp1_weight(size_t k, size_t N)
{
    // Weights for 1D-Simpson
    // [1, 4, 1] N=3
    // [1, 4, 2, 4, 1] N=5
    // [1, 4, 2, 4, 2, 5, 1] N=7, ...
    assert(N%2 == 1);
    if(k == 0 || k == N-1)
        return 1;
    if(k%2 == 0)
        return 2;
    return 4;
}

double simp2(int OVER_SAMPLING, double xc, double yc, double * H, double * R, double x0, double x1, double y0, double y1, int N)
{

    double dx = (x1-x0)/(double) N;
    double dy = (y1-y0)/(double) N;
    double v = 0;
    size_t sW = 0;
    for(int kk = 0; kk<N; kk++)
    {
        double wx = simp1_weight(kk, N);
        for(int ll = 0; ll<N; ll++)
        {
            double wy = simp1_weight(ll, N);
            double w = wx*wy;
            double x = x0+kk*dx;
            double y = y0+kk*dy;
            double r = sqrt(pow(x-xc, 2) + pow(y-yc, 2));
            int index = (int) floor(r*OVER_SAMPLING);
            sW+=w;
            v+=w*(H[index] + (H[index+1] - H[index]) * (r-R[index]));
        }
    }
    v = v/sW; // Pretty sure we don't need to count them ...

    return v;
}

void BW_slice(float * V, float z, bw_conf * conf)
{
    // Calculate the PSF for a single plane/slice
    // Todo: incorporate the integration element, dx*dy*dz
    // So that actual values are comparable

    // V is a pointer to the _slice_ not the whole volume
    // The center of the image in units of [pixels]
    double x0 = (conf->M - 1) / 2.0;
    double y0 = (conf->N - 1) / 2.0;


    int maxRadius = (int) round(sqrt(pow(conf->M - x0, 2) + pow(conf->N - y0, 2))) + 1;
    int OVER_SAMPLING = 10;

    size_t nr = maxRadius*conf->oversampling_R;
    double * r = malloc(nr * sizeof(double));
    double * h = malloc(nr * sizeof(double));

    /* Calculate the radial profile
       possibly by integrating over Z */

    for(size_t n = 0; n < nr; n++)
    {
      r[n] = ((double) n) / ((double) OVER_SAMPLING);
    }

    if(conf->Simpson == 0)
    {
        for (size_t n = 0; n < nr; n++) {
            h[n] = calculate(r[n] * conf->resLateral, z, conf);
        }
    }

    if(conf->Simpson == 1)
    {
        int NS = conf->Simpson_N;
    for (size_t n = 0; n < nr; n++) {
            h[n] = 0;
            double last = calculate(r[n] * conf->resLateral, z, conf);
            double now = last;

            last = now;
            now = 0;
            double dz = (conf->resAxial)/(NS-1);
            size_t sW = 0;
            for(int kk = 0; kk < NS; kk++)
            {
                double zs = z + kk*dz -(NS-1)/2*dz;
                double w = simp1_weight(kk, NS);
                sW += w;
                now += w*calculate(r[n] * conf->resLateral, zs, conf);
            }
            now = now/sW;
            h[n] = now;
        }
    }

    if(conf->Simpson == 2) // Using Li
    {
        int NS = conf->Simpson_N;
        memset(h, 0, nr*sizeof(double));
        double dz = (conf->resAxial)/(NS-1);

        for(int kk = 0; kk < NS; kk++) // Over z
        {
            double zs = z + kk*dz - (NS-1)/2*dz;
            double w = simp1_weight(kk, NS);

            li_conf * L = li_new(zs*1e9);
            L->lambda = conf->lambda*1e9;
            L->NA = conf->NA;
            L->ni = conf->ni;
            for(size_t n = 0; n < nr; n++) // over r
            {
                //printf("lambda = %e, r = %e\n", L->lambda, r[n]*conf->resLateral*1e9);
                h[n] += w*li_calc(L, r[n]*conf->resLateral*1e9)/(double) NS;
            }
            li_free(&L);
        }
        double wtot = 0;
        for(int kk = 0; kk < NS; kk++)
            wtot += simp1_weight(kk, NS);
        for(size_t n = 0; n < nr; n++)
            h[n]/=wtot;
    }

    if(0){
    for(size_t kk = 0; kk < nr; kk++)
    {
        printf("%f ", h[kk]);
    }
    printf("\n");
    exit(1);
    }
    assert(r[0] == 0);
    //exit(1);
    for (int x = 0; 2*x <= conf->M; x++) {
        for (int y = 0; 2*y <= conf->N; y++) {
            double v = 0;
            if(conf->Simpson > 0)
            {
                v = simp2(OVER_SAMPLING, x0, y0, h, r, x - 0.5, x + 0.5, y - 0.5, y + 0.5, conf->Simpson_N);
            }
            else
            {
                // radius of the current pixel in units of [pixels]
                double rPixel = sqrt(pow((double) x - x0, 2) + pow((double) y - y0, 2));
                // Index of nearest coordinate from below (replace by pixel integration)
                size_t index = (int) floor(rPixel * OVER_SAMPLING);
                assert(index < nr);
                // Interpolated value.
                v = h[index] + (h[index + 1] - h[index]) * (rPixel - r[index]) * OVER_SAMPLING;
            }


            // Use symmetry to fill the output image
            int xf = conf->M-x-1;
            int yf = conf->N-y-1;

            V[x + conf->M * y] = v;
            V[y + conf->M * x] = v;
            V[xf + conf->M * y] = v;
            //      V[yf + conf->M * x] = v;
            V[x + conf->M * yf] = v;
            //      V[y + conf->M * xf] = v;
            V[xf + conf->M * yf] = v;
            V[yf + conf->M * xf] = v;


        }
    }
    free(r);
    free(h);
    return;
}

int main(int argc, char ** argv)
{
    struct timespec tstart, tend;
    clock_gettime(CLOCK_REALTIME, &tstart);

    // Use defaults
    bw_conf * conf = bw_conf_new();
    bw_argparsing(argc, argv, conf);


    if( conf->overwrite == 0 && file_exist(conf->outFile))
    {
        printf("%s already exist. Doing nothing\n", conf->outFile);
        exit(0);
    }

    if(conf->verbose > 0)
    {
        bw_conf_printf(stdout, conf);
    }
    conf->log = fopen(conf->logFile, "w");
    if(conf->log == NULL)
    {
        printf("ERROR: Failed to open %s for writing\n", conf->logFile);
        exit(-1);
    }
    fprint_time(conf->log);
    bw_conf_printf(conf->log, conf);

    if(0){
        printf("i(rho=.5, r=0, z=0) = %f + %fi\n", creal(integrand(.5, 0, 0, conf)), cimag(integrand(0,0,0,conf)));
        printf("I(r = 0, z = 0) = %f\n", calculate(0,0, conf));
        printf("I(r = 300, z = 0) = %f\n", calculate(300*1e-9,0, conf));
        printf("I(r = 600, z = 0) = %f\n", calculate(600*1e-9,0, conf));
        exit(1);
    }
    // Run
    BW(conf);

    // Write to disk
    if(conf->verbose > 0) {
        printf("Writing as 32-bit floats to %s\n", conf->outFile);
    }

    fim_tiff_write_float(conf->outFile, conf->V, conf->M, conf->N, conf->P, conf->log);

    fprint_time(conf->log);
    clock_gettime(CLOCK_REALTIME, &tend);
    fprintf(conf->log, "Took: %f s\n", timespec_diff(&tend, &tstart));

    fprintf(conf->log, "done!\n");

    // Clean up
    free(conf->V);
    fclose(conf->log);
    free(conf->outFile);
    free(conf->logFile);
    free(conf->cmd);
    free(conf);


    return 0;
}
