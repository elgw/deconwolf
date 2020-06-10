#include "dw_bwpsf.h"

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
}


bw_conf * bw_conf_new()
{
  bw_conf * conf = malloc(sizeof(bw_conf));
  conf->lambda = 600*1e-9;
  conf->NA = 1.45;
  conf->ni = 1.52;
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

void usage(int argc, char ** argv)
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
    { NULL,           0,                 NULL,   0   }
  };
  int ch;
  while((ch = getopt_long(argc, argv, "vho:t:p:w:n:i:x:z:N:l:", longopts, NULL)) != -1)
  {
    switch(ch) {
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
        s->NA, s->ni, s->lambda*1e-9, s->resLateral*1e9, s->resAxial*1e9);
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
  float BesselValue = j0(k0 * conf->NA * r * rho);

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

  float del = 0.5; // integration interval
  float curDifference = conf->TOL; // Stopping criterion

  complex float value = 0;

  float curI = 0.0, prevI = 0.0;

  // Initialization of the Simpson sum (first iteration)
  int N = 2; // number of sub-intervals
  int k = 0; // number of consecutive successful approximations
  float rho = 0.5;

  complex float sumOddIndex = integrand(rho, r, defocus, conf);
  complex float sumEvenIndex = 0;

  complex float valueX0 = integrand(0.0, r, defocus, conf);
  complex float valueXn = integrand(1.0, r, defocus, conf);

  float complex sum = valueX0 + 2*sumEvenIndex + 4*sumOddIndex + valueXn;
  curI = (pow(creal(sum),2) + pow(cimag(sum), 2)) * pow(del,2);

  prevI = curI;

  // Finer sampling grid until we meet the TOL value with the specified
  // number of repetitions, K
  size_t iteration = 1;
  while (k < conf->K && iteration < 10000) {
    iteration++;
    N *= 2;
    del /= 2;

    sumEvenIndex = sumEvenIndex + sumOddIndex;
    sumOddIndex = 0 + I*0;

    for (int n = 1; n < N; n = n + 2) {
      rho = n * del;
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
    //    if(fabs(curI-prevI) < 1e-12)
    //    {
    //      break;
    //    }
    prevI = curI;
    //    printf("Iteration: %zu, curDifference: %e, curI: %e\n", iteration, curDifference, curI);
  }
  //  printf("\n");
  //  printf("@%d: curI = %f\n", k, curI);
  return curI;
}


void BW_slice(float * V, float z, bw_conf * conf)
{
  // V is a pointer to the _slice_ not the whole volume
  // The center of the image in units of [pixels]
  float x0 = (conf->M - 1) / 2.0;
  float y0 = (conf->N - 1) / 2.0;

  // Radial locations.
  // float xpAbs = Math.abs(xp), ypAbs = Math.abs(yp);
  // float maxRadialDistanceInPixels =
  // Math.round(Math.sqrt((xpAbs+nx-x0)*(xpAbs+nx-x0)+(ypAbs+ny-y0)*(ypAbs+ny-y0)))+1;
  int maxRadius = (int) round(sqrt(pow(conf->M - x0, 2) + pow(conf->N - y0, 2))) + 1;
  int OVER_SAMPLING = 1;
  size_t nr = maxRadius*OVER_SAMPLING;
  float * r = malloc(nr * sizeof(float));
  float * h = malloc(nr * sizeof(float));

  for (int n = 0; n < nr; n++) {
    r[n] = ((float) n) / ((float) OVER_SAMPLING);
    h[n] = calculate(r[n] * conf->resLateral, z, conf);
    //    printf("r[%d=%e] = %e, h[%d] = %e\n", n, r[n]*conf->resLateral, r[n], n, h[n]);
  }

  assert(r[0] == 0);
  //exit(1);
  for (int x = 0; 2*x <= conf->M; x++) {
    for (int y = 0; 2*y <= conf->N; y++) {
      // radius of the current pixel in units of [pixels]
      float rPixel = sqrt(pow(x-x0, 2) + pow(y-y0, 2));
      // Index of nearest coordinate from below (replace by pixel integration)
      int index = (int) floor(rPixel * OVER_SAMPLING);
      assert(index < nr);
      // Interpolated value.
      float v = h[index] + (h[index + 1] - h[index]) * (rPixel - r[index]) * OVER_SAMPLING;

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

  // Run
  BW(conf);

  // Write to disk
  if(conf->verbose > 0) {
    printf("Writing as 32-bit floats to %s\n", conf->outFile);
  }

  fim_tiff_write_float(conf->outFile, conf->V, conf->M, conf->N, conf->P);
  
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
