// Based on https://web.engr.oregonstate.edu/~mjb/cs575/Handouts/opencl.reduction.1pp.pdf
// Modified so that it works when the number of elements isn't a multiple
// of the work-group size.

// see https://dean-shaff.github.io/blog/c++/opencl/2020/03/29/opencl-reduction-sum.html

kernel void idiv_kernel( global const float * forward,
                         global const float * image,
                          local float *wgBuff, /* For this work group only */
                          global float *partialSums )
{
    int gid = get_global_id( 0 ); // 0 .. total_array_size-1
    int numItems = get_local_size( 0 ); // # work-items per work-group
    int tnum = get_local_id( 0 ); // thread (i.e., work-item) number in this work-group
// 0 .. numItems-1
    int wgNum = get_group_id( 0 ); // which work-group number this is in
    wgBuff[ tnum ] = 0.0;

    size_t p = gid / (wM * wN);
    size_t rem = gid - p*wM*wN;
    size_t n = rem / wM;
    size_t m = rem - n*wM;

    size_t imIdx = m + n*M + p*M*N;
    if(m < M0 && n < N && p < P )
    {
        float est = forward[gid];
        float obs = image[imIdx];
        if(est > 0 && obs > 0)
        {
            wgBuff[ tnum ] = obs*log(obs/est) - obs + est;
        } else {
            wgBuff[ tnum ] = 0.0;
        }
    }

// all threads execute this code simultaneously:
    for( int offset = 1; offset < numItems; offset *= 2 )
    {
        int mask = 2*offset - 1;
        barrier( CLK_LOCAL_MEM_FENCE ); // wait for all threads to get here
        if( ( tnum & mask ) == 0 ) // bit-by-bit and’ing tells us which
        { // threads need to do work now
            wgBuff[ tnum ] += wgBuff[ tnum + offset ];
        }
    }
    barrier( CLK_LOCAL_MEM_FENCE );
    if( tnum == 0 )
        partialSums[ wgNum ] = wgBuff[ 0 ];
}
