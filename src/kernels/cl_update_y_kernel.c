kernel void update_y_kernel( global float * y,
                             global const float * image)
{
    int gid = get_global_id(0);
    size_t p = gid / (wM * wN);
    size_t rem = gid - p*wM*wN;
    size_t n = rem / wM;
    size_t m = rem - n*wM;

    const float mindiv = 1e-6;
    size_t imIdx = m + n*M + p*M*N;
    if(m < M0 && n < N && p < P )
    {
        if(y[gid] < mindiv)
        {
            y[gid] = mindiv;
        }
        y[gid] = image[imIdx]/y[gid];
    } else {
        y[gid] = 0;
    }
}
