__kernel void cl_real_mul_inplace(const __global float * A,
                                  __global float * B)
{
    const size_t idx = get_global_id(0);

    if(idx < wMNP) /* Needed when % wgsize != 0 */
    {
        B[idx] = A[idx]*B[idx];
    }
}
