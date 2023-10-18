__kernel void cl_positivity(__global float * A,
                            const __global float * threshold)
{
    const size_t idx = get_global_id(0);

    if(idx < wMNP)
    {
        if(A[idx] < threshold[0])
        {
            A[idx] = threshold[0];
        }
    }
}
