# Installation notes

 - [Ubuntu 16.04](#ubuntu-1604)
 - [CentOS](#CentOS)
 - [FreeBSD](#FreeBSD)
 - [MKF FFT Backend](#MKL)

## CentOS
Tested on CentOS Linux Release 7.8.2009 (Core).

``` shell
# Install dependencies
sudo yum install gcc gsl-devel libtiff-devel fftw-devel
```

Should work with the
``` shell
make
sudo make install
```

``` shell
# To build with meson
sudo yum python3
sudo python3 -m pip install meson ninja
# Then follow the general meson instructions in README.md
```

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
Differences to building on Linux.
 - `pkgconf` replaces `pkg-config`
 - `gmake` replaces `make`
 - clang is the default compiler.

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
