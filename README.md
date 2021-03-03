# deconwolf

`deconwolf` is a program for deconvolution of fluorescent wide-field image stacks:
 - The deconvolved images shows very mild boundary effects which means that you can crop and deconvolve small regions of interest.
 - RAM usage can be reduced at the cost of slightly longer computation times by tiling. That makes it possible to deconvolve large images on small machines.
 - It can make use of all precious cores of your "big" machine since the critical parts run on separate threads (as many as you would like).
 - Deconwolf is tiny! The installation size of less than 0.1 MB it would fit most floppy drives (if you are fortunate enough to own one of those antiquities).

Except for this README there is also a [CHANGELOG](CHANGELOG.md) and a [TODO](TODO.md) list.

It does not:
 - Show your images, for that, try [ImageJ](https://imagej.net/Welcome)
 - Diagnose your imaging system.
 - Estimate your PSF based on real images.
 ...

There is also a [GUI](https://github.com/elgw/dw_gui) that could be handy.

At the moment we don't maintain any packages so you will have to build the program from the source code.

## Command line usage:
deconwolf has a simple command line interface (CLI) and only need to know which image that you want to deconvolve and what PSF that should be used. For example, to deconvolve the image `dapi_001.tif` by the PSF in `PSF_dapi.tif`, just type:
``` shell
dw dapi_001.tif PSF_dapi.tif
```
which will produce a new image called `dw_dapi_001.tif` along with a log file called `dw_dapi_001.tif.log.txt`. For the other options, see:
``` shell
dw --help
```

### Test data
No special test data has been prepared, but you can get some images from the [DeconvolutionLab2](http://bigwww.epfl.ch/deconvolution/deconvolutionlab2/) web page.

### Memory considerations
The peak memory usage is written at the end of the log file.

### Point Spread Function (PSF)
Theoretical PSFs according to the "Born-Wolf" model can be generated with `dw_bw`. As an example:

``` shell
dw_bw --lambda 466 --resxy 65 --resz 250 PSF_dapi.tif
```

In the current release deconwolf requires that the PSF is centered, i.e., that the largest value is in the middle. Consequently it prefers PFSs that have an odd size in each dimension.

### Supported image formats
Currently deconwolf does only support tif images, specifically: multipage, 16-bit unsigned or 32-bit floats, written in strip mode. The output is 16-bit unsigned and can be read by Matlab, ImageJ, etc.

### Log files and output
The reported error is the mean square error between the input image and the current guess convolved with the PSF.

## Building and installing
To build deconwolf on Ubuntu 20.04 you will need these packages: `fftw3f`, `fftw3f_threads` `openmp` and `tiff-5`, `libgsl-dev`.

Typical installation procedure:
```
make -B
sudo make install
```

On OSX you will need XCode from the App Store and [brew](https://brew.sh/). After that you need to install the following packages. Please install them one by one so that you will notice any errors.
```
xcode-select --install
brew install libopenmpt # Not sure if this is needed
brew install libomp
brew install libtiff
brew install fftw
brew install gsl
```

On Ubuntu 19.10 you might have to:
```
sudo apt-get update

# sudo apt-cache search libtiff
sudo apt-get install libtiff5-dev

# sudo apt-cache search libfftw3
sudo apt-get install libfftw3-dev
```

If you need/want to build fftw3 from source, that was not too tricky:
```
# download source first ...
./configure --help
./configure --enable-threads --enable-single
make
sudo make install
```

## Notes
 * FFTW is self tuning and will perform some tuning every time it presented for a new problem size. The result of this tuning is called wisdom and is stored in files like `fftw_wisdom_float_threads_16.dat` by deconwolf (in `~/config/deonwolf/`). Do not transfer that file to other machines and expect the tuning to take some time.


## Resources and references
 * [fftw3 documentation](http://www.fftw.org/fftw3_doc/).
 * [libtiff repository](https://gitlab.com/libtiff/libtiff)

The algoritm is based on these papers:

 * Richardson, William Hadley (1972). "Bayesian-Based Iterative Method of Image Restoration". JOSA. 62 (1): 55–59. [doi](https://doi.org/10.1364/JOSA.62.000055)
 * Lucy, L. B. (1974). "An iterative technique for the rectification of observed distributions". Astronomical Journal. 79 (6): 745–754. [doi](https://doi.org/10.1086%2F111605)
 * Biggs, D.S.C. “Acceleration of Iterative Image Restoration Algorithms.” Applied Optics. Vol. 36. Number 8, 1997, pp. 1766–1775.
 * M. Bertero and P. Boccacci, A simple method for the reduction of boundary effects in the Richardson-Lucy approach to image deconvolution,
A&A 437, 369-374 (2005), [doi](https://doi.org/10.1051/0004-6361:20052717)
 * Lee, Ji-Yeon & Lee, Nam-Yong. (2014). Cause Analysis and Removal of Boundary Artifacts in Image Deconvolution. Journal of Korea Multimedia Society. 17. 838-848. [doi](https://doi.org/10.9717/kmms.2014.17.7.838).
