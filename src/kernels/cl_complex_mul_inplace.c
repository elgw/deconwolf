// wMNP defined at compile time
__kernel void cl_complex_mul_inplace(const __global float * A,
                                     __global float * B)
{
    // B = A.*B
    const int id = get_global_id(0);
    if(id < cMNP)
    {
        float re = A[2*id]*B[2*id] - A[2*id+1]*B[2*id+1];
        float im = A[2*id]*B[2*id+1] + A[2*id+1]*B[2*id];
        B[2*id]=re;
        B[2*id+1]=im;
    }
}
