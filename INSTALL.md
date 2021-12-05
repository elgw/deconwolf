# Installation notes

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


## Using Intel MKL for the FFT
We've build deconwolf against the Intel MKL library (for the FFT) on Ubuntu,
although FFTW3 is the default. In initial performance tests MKL is found
to be faster than FFTW3.

### Installation
Install the required package(s), on Ubuntu:

``` shel
sudo apt install intel-mkl
```

Build using

``` shell
make MKL=1 -B
```
that will find the MKL libraries using
`pkg-config mkl-static-lp64-seq --cflags --libs`, the procedure
might be different on other platforms.

### Usage
To set the number of threads, set the environmental variable
`MKL_NUM_THREADS`, for example:
``` shell
export MKL_NUM_THREADS=8
dw ...
```
