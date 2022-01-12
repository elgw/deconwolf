# Installation notes

 - [Ubuntu 16.04](#ubuntu-1604)
 - [FreeBSD](#FreeBSD)
 - [MKF FFT Backend](#MKL)

## Ubuntu 16.04
``` shell
sudo apt-get update
# find out actual names with command like
# sudo apt-cache search fftw
sudo apt-get install libfftw3-single3
sudo apt-get install libfftw3-dev
sudo apt-get install libgsl-dev
sudo apt-get install libomp-dev
sudo apt-get install libtiff-dev
```
## FreeBSD
Deconwolf has been built on FreeBSD. The makefile does not
work out of the box and has to changed slightly. At least `pkg-config` has
to be replaced by `pkgconf` and `gmake`, not `make` should be used.

The following packages were required:
``` shell
pkg install gmake
pkg install fftw3
pkg install tiff
pkg install gsl
```

## MKL
FFTW3 is the default FFT backend for deconwolf but it is also possible to use
Intel MKL. This option is only tested on Ubuntu so far.

### Installation
Install the required package(s):

``` shel
sudo apt install intel-mkl
```

Build using

``` shell
make MKL=1 -B
```

Make will find the MKL libraries using
`pkg-config mkl-static-lp64-seq --cflags --libs`, the procedure
might be different on other platforms.

### Usage
To set the number of threads, set the environmental variable
`MKL_NUM_THREADS`, for example:
``` shell
export MKL_NUM_THREADS=8
dw ...
```
