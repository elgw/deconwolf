__kernel void cl_positivity(__global float * A,
                            const __global size_t * real_size,
                            const __global float * threshold)
{
    const size_t idx = get_global_id(0);

    if(idx < real_size[0])
    {
        if(A[idx] < threshold[0])
        {
            A[idx] = threshold[0];
        }
    }
}
