struct my_params {
    double NA;
    double ni;
    double defocus; // z
    double r;
    double lambda;
    size_t ncalls; /* Counting the number of function calls */
    size_t ninterval;      /* Max number of sub-intervals */
};

double my_f_real(double rho, void * p)
{
    struct my_params * params = (struct my_params *) p;
    params->ncalls++;
    double k0 = 2.0 * M_PI / params->lambda;
    double q = params->NA / params->ni;
    double bessel = j0f(k0*q*params->r*rho);
    double OPD = pow(q,2)*params->defocus*pow(rho, 2) / 2.0;
    double W = k0*OPD;

    return bessel*cos(W)*rho;
}

double my_f_imag(double rho, void * p)
{
    struct my_params * params = (struct my_params *) p;
    params->ncalls++;
    double k0 = 2.0 * M_PI / params->lambda;
    double q = params->NA / params->ni;
    double bessel = j0f(k0*q*params->r*rho);
    double OPD = pow(q,2)*params->defocus*pow(rho, 2) / 2.0;
    double W = k0*OPD;

    return -bessel*sin(W)*rho;
}

double integrate_bw_gsl(struct my_params * params, gsl_integration_workspace * workspace)
{

    gsl_function fun;
    fun.function = &my_f_real;
    fun.params = params;


    double result_real = 0;
    double abserr_real = 0;
    double result_imag = 0;
    double abserr_imag = 0;

    double epsabs = 1e-6; /* 1e-9 would not work ... */
    double epsrel = 1e-8;
    size_t limit = params->ninterval;
    int key = GSL_INTEG_GAUSS15; /* Default in MATLAB */
    gsl_integration_qag(&fun, 0, 1,
                        epsabs, epsrel,
                        limit, key, workspace,
                        &result_real, &abserr_real);
    fun.function = &my_f_imag;
    gsl_integration_qag(&fun, 0, 1,
                        epsabs, epsrel,
                        limit, key, workspace,
                        &result_imag, &abserr_imag);


    if(0){
    printf("result: %f, %f, abserr: %f, %f\n", result_real, result_imag,
           abserr_real, abserr_imag);
    }
    double result = pow(result_real, 2) + pow(result_imag, 2);
    return result;
}
