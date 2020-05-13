# deconwolf

`deconwolf` is a program for deconvolution of fluorescent wide-field image stacks, the highlights are:
 - almost no most boundary effects compared to ordinary linear deconvolution due to clever algorithm.
 - optional low RAM usage at the cost of slightly longer computation times which makes it possible to use low-cost hardware.
 - highly parallelised so that it can use all the juice in your fancy multi-core computer.
 - extremely tiny: with an installation size of less than 0.1 MB it would fit most floppy drives (if you are fortunate enough to own one of those antiquities).


Except for this readme there is also a [change log](CHANGELOG.md) and a [to do list](TODO.md).

## Usage:
deconwolf has a simple command line interface and only need to know which image that you want to deconvolve and what PSF that should be used. For example, to deconvolve the image `dapi_001.tif` by the PSF in `PSF_dapi.tif`, just type:
```
deconwolf dapi_001.tif PSF_dapi.tif
```
which will produce a new image called `dcw_dapi_001.tif` along with a log file called `dcw_dapi_001.tif.log.txt`. For the other options, see
```
deonwolf --help
```

In the source directory there is also a small script that might be useful for batch processing, example:
```
$ python3 deconwolf_batch.py iMERULA91_20200125_001/ ../PSF/ '--iter 50 --tilesize 512'
deconwolf --iter 50 --tilesize 512 iMERULA91_20200125_001/Cy5_012.tif ../PSF/PSF_Cy5.tif
deconwolf --iter 50 --tilesize 512 iMERULA91_20200125_001/Cy5_013.tif ../PSF/PSF_Cy5.tif
deconwolf --iter 50 --tilesize 512 iMERULA91_20200125_001/dapi_013.tif ../PSF/PSF_dapi.tif
...
```
the output can be run directly with `| bash` or piped to a file and executed later (possibly using `parallel`). Since deconwolf does not overwrite output files such list of jobs to be run can be aborted and restarted.

### Test data

### Memory considerations
The peak memory usage is written at the end of the log file.

### PSF
PSFs can be generate from ImageJ with a [plugin](http://bigwww.epfl.ch/algorithms/psfgenerator/). If the image has N slices, it is recommended that the PSF has 2xN-1 slices. Unfortunately that will require a lot of memory and processing out of your machine. If the PSF is larger than needed it will be cropped automagically.

### Supported image formats
Currently deconwolf does only support tif images, specifically: multipage, 16-bit unsigned or 32-bit floats, written in strip mode. The output is 16-bit unsigned and can be read by Matlab, ImageJ, etc.

### Log files and output
The reported error is the mean square error between the input image and the current guess convolved with the PSF.

## Building and installing
deconwolf requires `fftw3f`, `fftw3f_threads` and `tiff-5` and can be built with [meson](https://mesonbuild.com/). 

Typical installation procedure:
```
meson builddir
cd builddir
ninja 
# to install deconwolf to a standard location, use
sudo ninja install
# if you for some reason don't want it anymore, use
# sudo ninja uninstall
```

On OSX, if you have [homebrew](https://brew.sh/), then you can install meson with
```
brew install meson
```

On Ubuntu 19.10 I managed to get meson and the other libraries by:
```
sudo apt-get update

# sudo apt-cache search libtiff 
sudo apt-get install libtiff5-dev

# sudo apt-cache search libfftw3
sudo apt-get install libfftw3-dev

sudo apt-get install meson
```

If you need to build fftw3 from source, that was not too tricky:
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

