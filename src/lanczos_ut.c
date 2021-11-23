#include "lanczos.h"


/*
 * gcc -Wall -Wextra lanczos_ut.c lanczos.c -lm -o lanczos_ut
*/
int main(int argc, char ** argv)
{
    int nV = 6;
    double * v = malloc(nV*sizeof(double));
    for(int kk = 0; kk<nV; kk++)
    {
        v[kk] = kk;
        printf("v[%d] = %f\n", kk, v[kk]);
    }

    for(double x = 0; x<1; x+=0.1)
    {
        printf("l3(%f) = %f, l5(%f) = %f\n",
               x, lanczos3(v, nV, x),
               x, lanczos5(v, nV, x));
    }
    exit(EXIT_SUCCESS);
}
