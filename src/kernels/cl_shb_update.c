__kernel void cl_shb_update(__global float * P,
                            const __global float * x,
                            const __global float * xp,
                            const __global float * alpha)
{
    const size_t idx = get_global_id(0);

    if(idx < wMNP)
    {
        P[idx] = x[idx] + alpha[0]*(x[idx]-xp[idx]);
    }
}
