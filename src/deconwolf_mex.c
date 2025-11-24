/* This file provides a matlab interface to a few of the internal
 * functions of deconwolf (or libdeconwolf) and is used only for
 * internal testing at this point.
 *
 * Requires that libdeconwolf is installed. And unless the
 * library paths are set up correct you might need to start matlab with
 * something like
 * $ LD_LIBRARY_PATH=/usr/local/lib matlab
 *
 * To build (in MATLAB) use something like
 * >> mex('deconwolf_mex.c', '-L/usr/local/lib/', '-ldeconwolf', ...
 *    '-lgsl', '-lgslcblas', '-lm', '-lfftw3', '-lfftw3f', '-lfftw3f_omp', ...
 *    '-ltiff', '-output', 'deconwolf')
 *
 * check with pkg-config to get it right.
 *
 */

#include "mex.h"
#include <deconwolf/fim.h>

typedef struct{
    float * image;
    int M;
    int N;
    int P;
} matlab_image;

static void return_scalar(double v, mxArray * plhs[])
{
    const mwSize outDims[1] = {1};
    plhs[0] = mxCreateNumericArray(1, outDims, mxDOUBLE_CLASS, mxREAL);
    double * P = mxGetPr(plhs[0]);
    P[0] = v;
    return;
}

matlab_image get_matlab_float_image(const mxArray * A)
{
    if( ! mxIsSingle( A ) )
    {
        mexErrMsgTxt("Image is not of type float");
    }


    int ndim = mxGetNumberOfDimensions( A );
    if(ndim > 3)
    {
        mexErrMsgTxt("Image has more than 3 dimensions\n");
    }

    if(ndim == 0)
    {
        mexErrMsgTxt("Image has 0 dimensions\n");
    }

    matlab_image mi = {};
    mi.N = 1;
    mi.P = 1;

    const mwSize * dims = mxGetDimensions(A);

    mi.M = dims[0];
    if(ndim > 1)
    {
        mi.N = dims[1];
    }
    if(ndim > 2)
    {
        mi.P = dims[2];
    }

    mi.image = (float*) mxGetPr(A);

    return mi;
}

void matlab_return_image(mxArray *plhs[], int lhs_idx, const float * V,
                    const size_t M, const size_t N, const size_t P)
{
    const mwSize outDims[3] = {M, N, P};

    plhs[lhs_idx] = mxCreateNumericArray(3, outDims, mxSINGLE_CLASS, mxREAL);
    float * O = (float *) mxGetPr(plhs[lhs_idx]);

    for(size_t kk = 0; kk<M*N*P; kk++)
    {
        O[kk] = V[kk];
    }
    return;
}

void
mexFunction(int nlhs, mxArray *plhs[],
            int nrhs, const mxArray *prhs[])
{
    if(nrhs < 1)
    {
        mexErrMsgTxt("Usage:\n\t"
                     "deconwolf('function', args ...)\n");
    }

    if(!mxIsChar(prhs[0]))
    {
        mexErrMsgTxt("First argument, i.e. the function name must be a string\n");
    }

    if (mxGetM(prhs[0])!=1)
        mexErrMsgIdAndTxt( "MATLAB:revord:inputNotVector",
                           "Input must be a row vector.");

    char * cmd = (char*) mxArrayToString(prhs[0]);

    if(cmd == NULL)
    {
        mexErrMsgTxt("Error getting the first argument (string)\n");
    }


    if(strcmp(cmd, "dot_lateral_circularity") == 0)
    {
        if(nrhs < 3)
        {
            mexErrMsgTxt("At least two arguments required");
        }

        matlab_image im = get_matlab_float_image(prhs[1]);
        //printf("Got %d x %d x %d\n", im.M, im.N, im.P);
        double circ = fim_dot_lateral_circularity(im.image,
                                                  im.M, im.N, im.P,
                                                  10, 10, 0,
                                                  3);
        return_scalar(circ, plhs);
        return;
    }

    if(strcmp(cmd, "DoH") == 0)
    {
        if(nrhs < 2)
        {
            mexErrMsgTxt("At least one arguments required args: image, sigmaxy, sigmaz");
        }
        matlab_image im = get_matlab_float_image(prhs[1]);
        float * DoH = fim_DoH(im.image,
                              im.M, im.N, im.P,
                              2, 2);
        matlab_return_image(plhs, 0,
                            DoH, im.M, im.N, im.P);
        free(DoH);
        return;
    }

    if(strcmp(cmd, "LoG") == 0)
    {
        if(nrhs < 2)
        {
            mexErrMsgTxt("At least one arguments required args: image, sigmaxy, sigmaz");
        }
        matlab_image im = get_matlab_float_image(prhs[1]);
        float * DoH = fim_LoG_S2(im.image,
                              im.M, im.N, im.P,
                              2, 2);
        matlab_return_image(plhs, 0,
                            DoH, im.M, im.N, im.P);
        free(DoH);
        return;
    }

    printf("%s is an unknown command", cmd);
    mexErrMsgTxt("Incorrect usage\n");
}
