# Usage
Deconwolf has a command line interface (CLI), i.e., you would typically run
it from a terminal. However, here is also a
[GUI](https://github.com/elgw/dw_gui) that might be handy.

## Command line usage:
To deconvolve an image you need an approximation of the PSF. If you don't
have one you can create one with `dw_bw` that is shipped with deconwolf.
The basic information that you need to know is
 * The numerical aperture, NA.
 * The refractive index of the immersion, ni. Common media are air (ni = 1) and oil (ni = 1.51) and silicon oil (ni = 1.405).
 * The size of the pixels in the image.
 * The distance between the images/planes in z.

For example, if you have NA = 1.45, ni = 1.515 pixel size = 130 nm,
distance between planes = 250 nm, generate a PSF (`PSF_DAPI.tif`) with:

``` shell
dw_bw --resxy 130 --resz 250 --lambda 461 --NA 1.45 --ni 1.515 PSF_DAPI.tif
```

To deconvolve the images, type:

``` shell
dw dapi_001.tif PSF_dapi.tif
```

that will produce a new image called `dw_dapi_001.tif` along with a log file
called `dw_dapi_001.tif.log.txt`. Since you deserve better than the
default settings, see other options:

``` shell
dw_bw --help
dw --help
```

Please note that deconwolf requires that the pixel size is the same for
both the PSF and the input image and does not read that from any metadata.

## Test data
No special test data has been prepared, but you can get some images from the [DeconvolutionLab2](http://bigwww.epfl.ch/deconvolution/deconvolutionlab2/) web page.

## Memory considerations
The peak memory usage is written at the end of the log file. If you have 16 GB of RAM, images up to [1024x1024x60] pixels should work without tiling.

## PSF considerations
deconwolf requires that the PSF is centered, i.e.,
that the largest value is in the middle.
Consequently it prefers PFSs that
have an odd size in each dimension.
If you generate the PSF with some
other program you might have to center it first.


## Supported image formats
Currently deconwolf does only support tif images,
specifically: multipage, 16-bit unsigned or 32-bit floats,
written in strip mode. The output is either 16-bit unsigned or 32-bit
floating point and can be read by Matlab, ImageJ, etc.
If you use 16-bit output note that images will be scaled in order to not be
saturated. The scaling value can be found at the end of the log files.

## Log files and output
 * The reported error is the mean square error between the input image
   and the current guess convolved with the PSF.
 * If 16-bit output is used, the scaling is reported in the log file like:
   ```
   $ cat dw_dapi_001.tif.log.txt | grep scaling
   scaling: 0.425100
   ```
   That means, in order to go back to absolute intensities, divide by
   that number.
 * If deconwolf finished normally you will find that it ends with something
   like:
   ```
   Iteration 19/20, MSE=2.481362e+04
   Iteration 20/20, MSE=2.235636e+04
   1.957566% pixels at bg level in the output image.
   scaling: 0.425100
   Took: 28.767448 s
   peakMemory: 4420728 kiB
   Fri Jan 22 16:42:22 2021
   ```
   If it doesn't, chances are that deconwolf run out of memory and crashed.


# Notes
 * FFTW is self tuning and will perform some tuning every time it presented for a new problem size. The result of this tuning is called wisdom and is stored in files like `fftw_wisdom_float_threads_16.dat` by deconwolf (in `~/config/deonwolf/`). Do not transfer that file to other machines and expect the tuning to take some time.
