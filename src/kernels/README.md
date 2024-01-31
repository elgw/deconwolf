This folder contains opencl kernels (programs to be run on the GPU)
for different steps in the deconvolution algorithms.

## For complex arrays

The complex algorithms processes `cMNP` elements, defined at compile
time.


### cl_complex_mul.c

`C = A.*B`

```
__kernel void cl_complex_mul(const __global float * A,
                            const __global float * B,
                            __global float * C)
```

### cl_complex_mul_conj.c

`C = A.*conj(B)`

```
__kernel void cl_complex_mul_conj(const __global float *B,
                                 const __global float * A,
                                 __global float * C)
```

### cl_complex_mul_conj_inplace.c
In-place means that the first argument is modified in this case.

`A = A.*conj(B)`

```
__kernel void cl_complex_mul_conj_inplace(const __global float *A,
                                         __global float * B)
```

### cl_complex_mul_inplace.c

`A = A.*B`
```
__kernel void cl_complex_mul_inplace(const __global float * A,
                                    __global float * B)
```

### cl_preprocess_image.c

 `||fft_PSF(kk)|| < value[0] ? fft_image(kk) = 0 : 0`

```
__kernel void cl_preprocess_image(__global float *fft_image,
                                const __global float * fft_PSF,
                                const __global float * value)
```


## For real arrays

Processes `wMNP` numbers, also knows `wM`, `wN` and `wP`. Does not
know if the first dimension is padded i.e. assumes normal strides.

Also knows the size of the original image, `M`, `N`, `P`.

Besides the `cl_idiv_kernel` these kernels can be tricked to work with
padded data as long as we set `wM` and `wMNP` to the padded
values. Doing so will also update the pixels in the padding region,
but that does not matter since those values will not be propagated or
used.

### cl_positivity.c

`A = max(A, threshold)`

```
__kernel void cl_positivity(__global float * A,
                           const __global float * threshold);
```

### cl_idiv_kernel.c
Calculates the iDiv between the forward projection and the input image.
```
kernel void idiv_kernel( global const float * forward,
                         global const float * image,
                          local float *wgBuff, /* For this work group only */
                          global float *partialSums );
```

## cl_real_mul_inplace.c
`B = A.*B`

```
__kernel void cl_real_mul_inplace(const __global float * A,
                                  __global float * B);
```

## cl_shb_update.c
`P = x + alpha*(x-xp)`

```
__kernel void cl_shb_update(__global float * P,
                            const __global float * x,
                            const __global float * xp,
                            const __global float * alpha);
```

## cl_update_y_kernel.c

`y=y/image` for pixels where y > 1e-6, for the other pixels `y = 0`;

```
kernel void update_y_kernel( global float * y,
                             global const float * image);
```
