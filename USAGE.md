# Usage
Deconwolf has a command line interface (CLI), i.e., you would typically run it from a terminal. If you prefer, there is also a [GUI](https://github.com/elgw/dw_gui) that has to be installed separately..

## Command line usage:
To deconvolve an image, you need to have a correponding point spread function (PSF), and that is also the basic needs of deconwolf. For example, to deconvolve the image `dapi_001.tif` by the PSF in `PSF_dapi.tif`, just type:

``` shell
dw dapi_001.tif PSF_dapi.tif
```
which will produce a new image called `dw_dapi_001.tif` along with a log file
called `dw_dapi_001.tif.log.txt`. For the other options, see:

``` shell
dw --help
```

Please note that deconwolf requires that the pixel size is the same for both the PSF and the input image.

## Test data
No special test data has been prepared, but you can get some images from the [DeconvolutionLab2](http://bigwww.epfl.ch/deconvolution/deconvolutionlab2/) web page.

## Memory considerations
The peak memory usage is written at the end of the log file. If you have 16 GB of RAM, images up to [1024x1024x60] pixels should work without tiling.

## Point Spread Function (PSF)
Theoretical PSFs according to the "Born-Wolf" model can be generated with `dw_bw`. As an example:

``` shell
dw_bw --lambda 466 --resxy 65 --resz 250 PSF_dapi.tif
```

In the current release deconwolf requires that the PSF is centered, i.e., that the largest value is in the middle. Consequently it prefers PFSs that have an odd size in each dimension.

## Supported image formats
Currently deconwolf does only support tif images, specifically: multipage, 16-bit unsigned or 32-bit floats, written in strip mode. The output is 16-bit unsigned and can be read by Matlab, ImageJ, etc.

## Log files and output
The reported error is the mean square error between the input image and the current guess convolved with the PSF.

# Notes
 * FFTW is self tuning and will perform some tuning every time it presented for a new problem size. The result of this tuning is called wisdom and is stored in files like `fftw_wisdom_float_threads_16.dat` by deconwolf (in `~/config/deonwolf/`). Do not transfer that file to other machines and expect the tuning to take some time.
