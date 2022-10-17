__kernel void cl_real_mul_inplace(const __global float * A,
                                  __global float * B,
                                  const __global size_t * real_size)
{
    const size_t idx = get_global_id(0);

    if(idx < real_size[0])
    {
        B[idx] = A[idx]*B[idx];
    }
}
