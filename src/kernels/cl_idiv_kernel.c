// Each work item will read nPerItem elements from the input data
// directly. Hence the number of work items should be N/nPerItem


// Sync free tail reduction.
void warpReduce(volatile local float* wgBuff, int tnum)
{
    wgBuff[tnum] += wgBuff[tnum + 32];
    wgBuff[tnum] += wgBuff[tnum + 16];
    wgBuff[tnum] += wgBuff[tnum + 8];
    wgBuff[tnum] += wgBuff[tnum + 4];
    wgBuff[tnum] += wgBuff[tnum + 2];
    wgBuff[tnum] += wgBuff[tnum + 1];
}


kernel void idiv_kernel( global const float * forward,
                         global const float * image,
                          local float *wgBuff, /* For this work group only */
                          global float *partialSums )
{
    // compile time constant:
    // nPerItem

    int gid = get_global_id( 0 ); // 0 .. total_array_size-1
    int numItems = get_local_size( 0 ); // # work-items per work-group
    int tnum = get_local_id( 0 ); // thread (i.e., work-item) number in this work-group
// 0 .. numItems-1
    int wgNum = get_group_id( 0 ); // which work-group number this is in

    // Initialize the local memory to 0
    wgBuff[ tnum ] = 0.0;

    // Part I : Transfer from global to local memory in an attempt to read
    // as sequential as possible.
    for(int kk = 0; kk<nPerItem; kk++)
    {
        // Advance by the work group size each iter
        size_t idx = nPerItem*wgNum*numItems + tnum + kk*numItems;

        size_t p = idx / (wM * wN);
        size_t rem = idx - p*wM*wN;
        size_t n = rem / wM;
        size_t m = rem - n*wM;

        if(m < M0 && n < N && p < P )
        {
            size_t imIdx = m + n*M + p*M*N;
            float est = forward[idx];
            float obs = image[imIdx];
            if(est > 0 && obs > 0)
            {
                wgBuff[ tnum ] += obs*log(obs/est) - obs + est;
            }
        }
    }

    // Part II : reduce the work group buffer

    // Sequential Addressing
    // read s items, write to s/2 items
    // if numItems = 8 (it is typically 256)
    // WWWW0000 s=8
    // WW000000 s=4
    // W0000000 s=2
    for (unsigned int s=numItems/2; s>32; s>>=1) {
        barrier( CLK_LOCAL_MEM_FENCE ); // wait for all threads to get here
        if (tnum < s) {
            wgBuff[tnum] += wgBuff[tnum + s];
        }
    }

    barrier( CLK_LOCAL_MEM_FENCE );
    if(tnum < 32)
    {
        warpReduce(wgBuff, tnum);
    }

    // Part II : Write the sum for this work group to the output
    if( tnum == 0 )
        partialSums[ wgNum ] = wgBuff[ 0 ];

}
