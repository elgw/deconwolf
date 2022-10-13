__kernel void cl_error_idiv(__global float * F, // Forward guess
                            __global float * I, // Reference image
                            __global float * args)
{
    // Assumes F_size[i] >= I_size[i], i=0,1,2
    // I.e. F is a larger image than I

    const int idx_I = get_global_id(0);

    // id is assumed to be a position inside I
    // then we need to figure out what position in F that it corresponds to

    int M = args[0];
    int N = args[1];
    int P = args[2];
    int wM = args[3];
    int wN = args[4];
    int wP = args[5];

    // [m, n, p] = ind2subs(idx_I)
    int idx = idx_I;
    int p = idx / (M*N);
    idx = idx - p*M*N;
    int m = idx / M;
    idx = idx - m*M;
    int n = idx;

    int idx_F = m + n*wM + p*wM*wN;
    float obs = F[idx_F];
    float est = I[idx_I];
    if(est > 0)
    {
        args[6] += 1.0; //obs*log(obs/est) - obs + est;
    }

}
