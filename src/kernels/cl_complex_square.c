__kernel void cl_complex_square(__global float *A)
{
    int id = get_global_id(0);


    float re = A[2*id]*A[2*id] - A[2*id+1]*A[2*id+1];
    float im = A[2*id]*A[2*id+1] + A[2*id+1]*A[2*id];
    A[2*id]=re;
    A[2*id+1]=im;

}
