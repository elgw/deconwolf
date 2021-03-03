# deconwolf

`deconwolf` is a program for deconvolution of fluorescent wide-field image stacks:
 - The deconvolved images shows very mild boundary effects which means that you can crop and deconvolve small regions of interest.
 - RAM usage can be reduced at the cost of slightly longer computation times by tiling. That makes it possible to deconvolve large images on small machines.
 - It can make use of all precious cores of your "big" machine since the critical parts run on separate threads (as many as you would like).
 - Deconwolf is tiny! It could even fit on a floppy drive (if you are fortunate enough to own one of those antiquities).

Except for this README there is also a short [manual](USAGE.md), a [CHANGELOG](CHANGELOG.md) and a [TODO](TODO.md) list.

It does not:
 - Show your images, for that, try [ImageJ](https://imagej.net/Welcome)
 - Diagnose your imaging system.
 - Estimate your PSF based on real images.
 ...

There is also a [GUI](https://github.com/elgw/dw_gui) that could be handy.

At the moment no pre-built packages are maintained so you will have to build deconwolf from the source code.

## Building and installing
Deconwolf runs on 64-bit machines with x86_64 architecture, it has been built and installed on Ubuntu 20.04 and macOS Big Sur.

### Ubuntu 20.04
Ensure that the required packages are installed

``` shell
sudo apt-get update
sudo apt-get install fftw3f
sudo apt-get install fftw3f_threads
sudo apt-get install openmp
sudo apt-get install tiff-5
sudo apt-get install libgsl-dev
```

To build and install:
``` shell
make -B
sudo make install
```

### macOS Big Sur
You will need XCode from the App Store and [brew](https://brew.sh/).

First install the required libraries:
``` shell
xcode-select --install
brew install libopenmpt # Not sure if this is needed
brew install libomp
brew install libtiff
brew install fftw
brew install gsl
```

To build and install:
``` shell
make -B
sudo make install
```

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
