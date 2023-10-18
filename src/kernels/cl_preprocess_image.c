__kernel void cl_preprocess_image(__global float *fft_image,
                                 const __global float * fft_PSF,
                                 const __global float * value)
{
    // C = A.*B
    const int id = get_global_id(0);

    if(id < cMNP)
    {
        float el_norm = fft_PSF[2*id]  *fft_PSF[2*id]
                      + fft_PSF[2*id+1]*fft_PSF[2*id + 1];
        el_norm = sqrt(el_norm);
        if(el_norm < value[0])
        {
            fft_image[2*id]=0.0;
            fft_image[2*id+1]=0.0;
        }
    }
}
