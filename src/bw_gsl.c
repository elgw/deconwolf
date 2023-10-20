#include "bw_gsl.h"

bw_gsl_conf_t * bw_gsl_new(size_t limit)
{
    bw_gsl_conf_t * conf = malloc(sizeof(bw_gsl_conf_t));
    assert(conf != NULL);
    conf->limit = limit;
    conf->ncalls = 0;
    conf->epsabs = 1e-6; /* 1e-9 would not work ... */
    conf->epsrel = 1e-6;

    /* These selections are based on benchmarks */
    conf->key = GSL_INTEG_GAUSS15;
    /* The larger (r,z) the better with many points */
    conf->keybw = GSL_INTEG_GAUSS61;

    conf->wspx = gsl_integration_workspace_alloc(conf->limit);
    conf->wspy = gsl_integration_workspace_alloc(conf->limit);
    conf->wspz = gsl_integration_workspace_alloc(conf->limit);
    conf->wspb = gsl_integration_workspace_alloc(conf->limit);

    return conf;
}

void bw_gsl_free(bw_gsl_conf_t * conf)
{
    gsl_integration_workspace_free(conf->wspx);
    gsl_integration_workspace_free(conf->wspy);
    gsl_integration_workspace_free(conf->wspz);
    gsl_integration_workspace_free(conf->wspb);
    free(conf);
    return;
}

void fprintf_GSL_INTEG_KEY(FILE * fid, int key)
{
switch(key)
{
case 1:
    fprintf(fid, "GSL_INTEG_GAUSS15\n");
    break;
case 2:
    fprintf(fid, "GSL_INTEG_GAUSS21\n");
    break;
case 3:
    fprintf(fid, "GSL_INTEG_GAUSS31\n");
    break;
case 4:
    fprintf(fid, "GSL_INTEG_GAUSS41\n");
    break;
case 5:
    fprintf(fid, "GSL_INTEG_GAUSS51\n");
    break;
case 6:
    fprintf(fid, "GSL_INTEG_GAUSS61\n");
    break;
default:
    fprintf(fid, "Key does not correpond to any GSL_INTEG_* key\n");
    exit(EXIT_FAILURE);
}
return;
}

void bw_gsl_fprint(FILE *fid, bw_gsl_conf_t *conf)
{
    fprintf(fid, "limit: %zu (sub intervals)\n", conf->limit);
    fprintf(fid, "epsabs: %e\n", conf->epsabs);
    fprintf(fid, "epsrel: %e\n", conf->epsrel);
    fprintf(fid, "spatial key: %d ", conf->key);
    fprintf_GSL_INTEG_KEY(fid, conf->key);
    fprintf(fid, "bw key: %d ", conf->keybw);
    fprintf_GSL_INTEG_KEY(fid, conf->keybw);
    return;
}

/* BW Integration */

/* Real part of bw */
static double my_f_real(double rho, void * p)
{
    bw_gsl_conf_t * params = (bw_gsl_conf_t *) p;
    params->ncalls++;
    double k0 = 2.0 * M_PI / params->lambda;

    double NA = params->NA;
    double ni = params->ni;

    double bessel = j0f(k0*NA*params->r*rho);
    // double OPD = pow(q,2)*params->z*pow(rho, 2) / 2.0; // old
    double OPD = pow(NA,2)*params->z*pow(rho, 2) / 2.0/ni; // new
    double W = k0*OPD;

    return bessel*cos(W)*rho;
}

/* imaginary part of bw */
static double my_f_imag(double rho, void * p)
{
    bw_gsl_conf_t * params = (bw_gsl_conf_t *) p;
    params->ncalls++;
    double k0 = 2.0 * M_PI / params->lambda;
    //double q = params->NA / params->ni;
    double NA = params->NA;
    double ni = params->ni;

    double bessel = j0f(k0*NA*params->r*rho);
    // double OPD = pow(q,2)*params->z*pow(rho, 2) / 2.0; // old
    double OPD = pow(NA,2)*params->z*pow(rho, 2) / 2.0/ni; // new
    double W = k0*OPD;
    return -bessel*sin(W)*rho;
}


double bw_gsl_integrate(bw_gsl_conf_t * conf, double r, double z)
{

    conf->r = r;
    conf->z = z;


    if(0) {
        /* Pinhole of 1 Airy Unit */
        double k0 = 2.0 * M_PI / conf->lambda;
        if(k0*conf->NA*r > 2.4048)
            return 0;
        if(z != 0)
            return 0;
    }

    double result_real = 0;
    double abserr_real = 0;
    double result_imag = 0;
    double abserr_imag = 0;

    gsl_function fun;
    fun.params = conf;

    fun.function = &my_f_real;
    gsl_integration_qag(&fun, 0, 1,
                        conf->epsabs, conf->epsrel,
                        conf->limit, conf->keybw, conf->wspb,
                        &result_real, &abserr_real);

    fun.function = &my_f_imag;
    gsl_integration_qag(&fun, 0, 1,
                        conf->epsabs, conf->epsrel,
                        conf->limit, conf->keybw, conf->wspb,
                        &result_imag, &abserr_imag);


    if(0){
    printf("result: %f, %f, abserr: %f, %f\n", result_real, result_imag,
           abserr_real, abserr_imag);
    }
    double result = pow(result_real, 2) + pow(result_imag, 2);
    return result;
}

/* Integrate BW integral over region */

static double my_f_z(double z, void * p)
{
    /* Function to be integrated over z */
    bw_gsl_conf_t * conf = (bw_gsl_conf_t *) p;
    double r = conf->r;
    double v = bw_gsl_integrate(conf, r, z);
    //printf("z=%f, v=%f\n", z, v);
    return v;
}

static double my_f_y(double y, void * p)
{
    /* Function to be integrated over z */
    bw_gsl_conf_t * conf = (bw_gsl_conf_t *) p;
    double r = sqrt( pow(conf->x,2) + pow(y, 2));
    double v = bw_gsl_integrate(conf, r, conf->z);
    return v;
}


static double my_f_z_xy(double z, void * p)
{
    /* Integrate xy over z */
    bw_gsl_conf_t * conf = (bw_gsl_conf_t *) p;

    double v = bw_gsl_integrate_xy(conf,
                                   conf->x0, conf->x1,
                                   conf->y0, conf->y1,
                                   z);
    return v;
}

static double my_f_x_y(double x, void * p)
{
    /* Integrate x  */
    bw_gsl_conf_t * conf = (bw_gsl_conf_t *) p;
    conf->x = x;
    double v = bw_gsl_integrate_y(conf,
                                  conf->x,
                                  conf->y0, conf->y1,
                                  conf->z);
    return v/(conf->x1-conf->x0);
}


double bw_gsl_integrate_z(bw_gsl_conf_t * conf,
                          double r, double z0, double z1)
{
/* Integration of the BW integral over a range of z (constant r)
 * Note: division by dz at the end, i.e. returns the average */

    conf->r = r;

    gsl_function fun;
    fun.function = &my_f_z;
    fun.params = conf;

    double result = 0;
    double abserr = 0;



    /* We ignore the status from gsl_integration_qag */
    gsl_integration_qag(&fun, z0, z1,
                        conf->epsabs, conf->epsrel,
                        conf->limit, conf->key, conf->wspz,
                        &result, &abserr);

    return result/(z1-z0);
}


/* Integrate III bw dx*dy*dz */
double bw_gsl_integrate_xyz(bw_gsl_conf_t * conf,
                            double x0, double x1,
                            double y0, double y1,
                            double z0, double z1)
{
    gsl_function fun;
    fun.function = &my_f_z_xy;
    fun.params = conf;
    conf->x0 = x0;
    conf->x1 = x1;
    conf->y0 = y0;
    conf->y1 = y1;
    conf->z0 = z0;
    conf->z1 = z1;

    double result = 0;
    double abserr = 0;

    gsl_integration_qag(&fun, z0, z1,
                        conf->epsabs, conf->epsrel,
                        conf->limit, conf->key, conf->wspz,
                        &result, &abserr);

    return result/(z1-z0);
}


double bw_gsl_integrate_xy(bw_gsl_conf_t * conf,
                           double x0, double x1,
                           double y0, double y1,
                           double z)
{
    gsl_function fun;
    fun.function = &my_f_x_y;
    fun.params = conf;

    conf->y0 = y0;
    conf->y1 = y1;
    conf->z = z;

    double result = 0;
    double abserr = 0;

    gsl_integration_qag(&fun, x0, x1,
                        conf->epsabs, conf->epsrel,
                        conf->limit, conf->key, conf->wspx,
                        &result, &abserr);

    return result;
}

double bw_gsl_integrate_y(bw_gsl_conf_t * conf,
                          double x, double y0, double y1, double z)
{

    conf->x = x;
    conf->z = z;

    gsl_function fun;
    fun.function = &my_f_y;
    fun.params = conf;

    double result = 0;
    double abserr = 0;

    /* We ignore the status from gsl_integration_qag */
    //printf("I %f, [%f, %f], %f\n", x, y0, y1, z);
    gsl_integration_qag(&fun, y0, y1,
                        conf->epsabs, conf->epsrel,
                        conf->limit, conf->key, conf->wspy,
                        &result, &abserr);

    return result/(y1-y0);
}
