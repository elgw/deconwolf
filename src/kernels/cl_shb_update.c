__kernel void cl_shb_update(__global float * P,
                            const __global float * x,
                            const __global float * xp,
                            const __global size_t * real_size,
                            const __global float alpha)
{
    const size_t idx = get_global_id(0);

    if(idx < real_size[0])
    {
        P[idx] = x[idx] + alpha[0]*(x[idx]-xp[idx]);
    }
}
