# deconwolf v0.3.7

**deconwolf** is a software for 3-D deconvolution of fluorescent wide-field
images:
 - The deconvolved images shows very mild boundary effects which means
   that relatively small regions of interest can be used.
 - RAM usage can be reduced at the cost of slightly longer computation times by
   tiling. That makes it possible to deconvolve large images on small machines.
 - Critical parts run on separate threads (as many as you would
   like). However, for maximal throughput (and if you have enough
   RAM), run several instances of dw in parallel.
 - It is shipped with program to generate Point Spread Functions (PSFs)
   according to the Born and Wolf model. Our program is the only one we
   know of that actually integrate the PSF over each pixel.
 - The programs are developed under Linux but can be compiled on Mac
   and Windows as well.
 - Fully open source. And we embrace [contributions and
   suggestions](CONTRIBUTING.md).


Deconwolf does not:
 - Have a full-featured Graphical User Interface (GUI), however, if you
   like to click buttons there is a limited
   [GUI](https://github.com/elgw/deconwolf-gui).
 - Show your images, for that, try [ImageJ](https://imagej.net/Welcome)
 - Diagnose your imaging system and correct for typical CMOS or CCD
   artifacts.
 - Estimate your PSF based on real images.
 - Tell you how many iterations that you should use although there are
   three ways to make it stop i) after a fixed number of iterations
   ii) when the progress is slow, or iii) at a fixed difference
   between the input image and the current guess (convolved with the PSF).
 - If you miss one or more of these features,
   or just want something else, there is an (incomplete)
   list of [alternative](#aternatives) software.

Except for this README.me there is also a short [USAGE.md](USAGE.md),
a [CHANGELOG.md](CHANGELOG.md) and a [TODO.md](TODO.md).

After building and installing you will find two binaries:
 - **dw** -- for deconvolution, [man page](doc/dw.txt)
 - **dw_bw** -- to create PSFs using the Born-Wolf model [man page](doc/dw_bw.txt)

At the moment we don't provide pre-built packages. You will have to build
deconwolf from the source code, instructions follows below.

## Build and install
Deconwolf runs on 64-bit machines with, both aarch64 and x86_64, and
require no special hardware. To compile and install deconwolf should
take less than a minute on a Linux machine but might be more
cumbersome on MacOS and Windows. For platform specific build
instructions, see [INSTALL.md](INSTALL.MD).  We hope to provide
pre-compiled binaries in the future.

### Dependencies
Deconwolf uses:

 * [FFTW](http://www.fftw.org/fftw3_doc/) for Fourier Transforms. Or
   alternatively Intel MKL.
 * [libtiff](https://gitlab.com/libtiff/libtiff) to read and write TIFF files.
 * [GNU Scientific Library](https://www.gnu.org/software/gsl/) for
   integration and special functions.
 * [OpenMP](https://www.openmp.org/) for _automatic_ parallelization
   of code, used in **dw**.

If these libraries are available for your platform, chances are that
deconwolf can be built as well.

To enable GPU acceleration, deconwolf also requires a GPU and OpenCL
drivers installed.

### Installation on Ubuntu 22.04 (or Windows via WSL)
For other platforms, see [INSTALL.md](INSTALL.md).

First install the required packages:

``` shell
sudo apt-get update
sudo apt-get install \
 gcc                 \
 libfftw3-single3    \
 libfftw3-dev        \
 libgsl-dev          \
 libomp-dev          \
 libpng-dev          \
 libtiff-dev         \
 pkg-config
```

Build, to build and install deconwolf:
``` shell
make -B
./makedeb-ubuntu_2204
sudo apt-install ./deconwolf_*.deb
# to remove
# sudo apt remove deconwolf
```

If OpenCL is not installed, use `make -B VKFFT=0` to build without GPU
acceleration.

## Minimal usage example
To generate a parametric PSF and deconvolve an image, all you need is
something like this:

``` shell
# generate PSF.tif
dw_bw --resxy 130 \
--resz 250 \
--NA 1.46 \
--ni 1.518 \
--lambda 460 PSF.tiff
# Deconvolve image.tiff -> dw_image.tiff
dw --iter 50 image.tiff PSF.tiff
```
For available options, please see

``` shell
dw --help
dw_bw --help
```

To validate that **dw** does not create random garbage, run it on
`/demo/dapi_001.tif`


``` shell
cd demo
make
imagej dapi_001.tif dw_dapi_001.tif
```

The run time on an AMD 3700x was 8s. To use GPU acceleration. To use
GPU accelerated deconvolution, test

``` shell
dw --iter 20 dapi_001.tif PSF_dapi.tif --gpu
```

For more documentation see the short [usage guide](USAGE.md), and the manual
pages for both binaries, [man dw](doc/dw.txt) [man dw_bw](doc/dw_bw.txt).

### Performance hints
By default deconwolf use a method to avoid the creation of boundary
artifacts [^1] which isn't available in other software. To be able to
test raw numerical performance we turn that option of for this section.

Benchmarking is performed on the [microtubules
image](https://bigwww.epfl.ch/deconvolution/data/microtubules/) using
the accompanying PSF.

Please note that it does not simulate "real" wide
field data very well since it is created by periodic convolution. In
this case it is necessary to use the options **--periodic** to
indicate that periodic deconvolution should be used, and also
**--xyfactor 0** do disable any cropping of the PSF.

System: Ubuntu 22.04.4 LTS, AMD Ryzen 7 3700X 8-Core Processor, 64 GB
RAM, 12 GB RX 6700 XT GPU.

| software                  | time (s) | self-mem (Mb) | sys-mem (Mb) |
| ------------------------  | -------- | --------      | -------      |
| DeconvolutionLab2         | 1025     | 1582          | 48511        |
| DeconvolutionLab2 + FFTW2 | 862      | 1353          | 47387        |
| MATLAB/deconvlucy         | 104      |               | 5270         |
| dw 1.3.7 --threads 1      | 52       |               |  344         |
| dw 1.3.7 --threads 2      | 32       |               |  419         |
| dw 1.3.7 --threads 4      | 21       |               |  566         |
| dw 1.3.7 --threads 8      | 18       |               | 1124         |
| dw 1.3.7 --gpu            |  2.6     |               | 5085         |

Notes:

- sys-mem is measured by parsing the **VmPeak** value from
  `/proc/pid/status`. In the case of DeconvolutionLab2 the values does
  not necessarily reflect the required memory since it is written in
  Java which is garbage collected. For MATLAB/deconvlucy the memory
  includes the full MATLAB environment.

- self-mem is the memory usage reported by the software if available.

- DeconlutionLab2 use vanilla RL, i.e. without any acceleration which
means that more iterations will be needed before convergence.

- MATLAB/deconvlucy use "Biggs" acceleration. Matlab version R2020b
  was used in this case.

- For "real" data, when **--periodic** is not used, the input image is
  padded automagically during processing and the relevant part is
  cropped out at the end.

### Bugs

Most likely there are bugs and they can only be fixed when they are known.
Please open a [new ticket](https://github.com/elgw/deconwolf/issues) if you
have any issues with the program.


## References

There is a [pseudo code](PSEUDOCODE.md) description of what the binaries does.

The deconvolution algorithm is based on the following papers:

 * Richardson, William Hadley (1972). "Bayesian-Based Iterative Method of Image
   Restoration". JOSA. 62 (1): 55–59.
   [doi](https://doi.org/10.1364/JOSA.62.000055)
 * Lucy, L. B. (1974). "An iterative technique for the rectification of observed
   distributions". Astronomical Journal. 79 (6): 745–754.
   [doi](https://doi.org/10.1086%2F111605)

   These together are referred to as Richardson-Lucy
   (RL).

* Wang H, et al. Scaled Heavy-Ball Acceleration of the
   Richardson-Lucy Algorithm for 3D Microscopy Image Restoration. IEEE
   Trans Image Process. 2014 [doi](https://doi.org/10.1109/TIP.2013.2291324).

   This is the default acceleration method as it use less memory than
   the other alternatives below, and seems to generate less shot noise.

[^1]: M. Bertero and P. Boccacci, A simple method for the reduction of boundary
   effects in the Richardson-Lucy approach to image deconvolution,
   A&A 437, 369-374 (2005).
   [doi](https://doi.org/10.1051/0004-6361:20052717)

   The default boundary handling method (corresponds to **\--bq 2**),
   although extended to 3D. To disable boundary handling or to use
   already padded data, set **\--bq 0** which turns off this feature
   completely, i.e., lets dw use periodic boundary handling. The
   option **\--bq 1** is a compromise of speed and memory vs quality.

The PSF calculations in **dw_bw** use these:

 * Max Born. Principles of optics : electromagnetic theory of propagation, interference,
   and diffraction of light. Cambridge: Cambridge University Press, 2019.
   ISBN: 978-1-108-47743-7.

 * F. Aguet. “Super-Resolution Fluorescence Microscopy Based on Physical
   Models”. EPFL Thesis no. 4418 (2009), 209 p. Swiss Federal Institute of
   Technology Lausanne (EPFL), May 2009
   [url](http://bigwww.epfl.ch/publications/aguet0903.html)

 * Jizhou Li, Feng Xue, and Thierry Blu. “Fast and accurate three-dimensional
   point spread function computation for fluorescence microscopy”. In: Journal
   of the Optical Society of America A 34.6 (May 2017), p. 1029.
   [doi](https://doi.org/10.1364/josaa.34001029)

   Enable by **\--li**.

 * [VkFFT](https://github.com/DTolm/VkFFT) is used for FFT transform
   on GPUs (via OpenCL). Until version 0.3.6
   [libclfft2](https://github.com/clMathLibraries/clFFT) was used for
   this purpose.


## Alternatives
This is a non-complete list of alternative deconvolution software:

### Deconvolution
Free and open source:
 - [Deconvolution Lab2](http://bigwww.epfl.ch/deconvolution/deconvolutionlab2/)

Commercial:
 - [Huygens](https://svi.nl/HomePage)
 - [Microvolution](https://www.microvolution.com/)

### Point spread functions

- [PSF Generator](http://bigwww.epfl.ch/algorithms/psfgenerator/)
