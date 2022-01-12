# deconwolf

 * [Introduction](#Introduction)
 * [Build and Install](#install)
   * [Dependencies](#dependencies)
   * [Ubuntu 20.04](#linux)
   * [macOS Big Sur](#osx)
   * [Windows 10](#win10)
   * [FreeBSD 13.0](#freebsd)
 * [Usage](#use)
   * [Bugs](#bugs)
 * [References](#references)
 * [Alternatives](#alternatives)

## Introduction
**deconwolf** is a program for 3-D deconvolution of fluorescent wide-field
images:
 - The deconvolved images shows very mild boundary effects which means that you
   can crop and deconvolve small regions of interest.
 - RAM usage can be reduced at the cost of slightly longer computation times by
   tiling. That makes it possible to deconvolve large images on small machines.
 - It can make use of all precious cores of your "big" machine since the
   critical parts run on separate threads (as many as you would like).
 - Deconwolf is tiny! The binaries could even fit on a floppy drive
   (if you are fortunate enough to own one of those antiquities).
 - It is shipped with program to generate Point Spread Functions (PSFs)
   according to the Born and Wolf model. Our program is the only one we
   know of that actually integrate the PSF over each pixel.
 - Fully open source. And we embrace [contributions and suggestions](CONTRIBUTING.md).

Except for this README.me there is also a short [USAGE.md](USAGE.md),
a [CHANGELOG.md](CHANGELOG.md) and a [TODO.md](TODO.md).

This repository provides two binaries,
 - `dw` -- for deconvolution, [man page](doc/dw.1.txt)
 - `dw_bw` -- to create PSFs using the Born-Wolf model [man page](doc/dw_bw.1.txt)

Deconwolf does not:
 - Have a full-featured Graphical User Interface (GUI), however, if you
   like to click buttons there is a limited
   [GUI](https://github.com/elgw/deconwolf-gui).
 - Show your images, for that, try [ImageJ](https://imagej.net/Welcome)
 - Diagnose your imaging system.
 - Estimate your PSF based on real images.
 - If you miss one or more of these features,
   or just want something else, there is an (incomplete)
   list of [alternatives](#aternatives).

At the moment we don't provide pre-built packages. You will have to build
deconwolf from the source code, instructions follows below.

<a name="install" />

## Build and install
Deconwolf runs on 64-bit machines with x86_64 architecture. Jump directly to
installation instructions for [Ubuntu](#linux) or [macOS](#osx),
[Windows 10](#win10) or [FreeBSD 13](#freebsd). Instruction for other
systems will be collected in [INSTALL.md](INSTALL.md).


### Dependencies
Deconwolf uses:

 * [FFTW](http://www.fftw.org/fftw3_doc/) for Fourier Transforms.
 * [libtiff](https://gitlab.com/libtiff/libtiff) to read and write TIFF files.
 * [GNU Scientific Library](https://www.gnu.org/software/gsl/) for
   integration and special functions.
 * [OpenMP](https://www.openmp.org/) for _automatic_ parallelization of code.
 * [POSIX Threads](https://en.wikipedia.org/wiki/Pthreads) for explicit paralellization.

If these libraries are available for your platform, chances are that that it can
be built.

<a name="linux" />

### Ubuntu 20.04

Install required packages:

``` shell
sudo apt-get update
# find out actual names with command like
# sudo apt-cache search fftw
sudo apt-get install libfftw3-single3
sudo apt-get install libfftw3-dev
sudo apt-get install openmp
sudo apt-get install tiff-5     # or possibly the one on the next line
sudo apt-get install libtiff-dev
sudo apt-get install libgsl-dev
```

Build and install deconwolf:
``` shell
make -B
sudo make install
```


<a name="osx" />

### macOS Big Sur

For building you will need XCode from the App Store and [brew](https://brew.sh/).

Then set up XCode and install the required packages:
``` shell
xcode-select --install
brew install libopenmpt # Not sure if this is needed
brew install libomp
brew install libtiff
brew install fftw
brew install gsl
```

Build and install deconwolf
``` shell
make -B
sudo make install
```

<a name="win10" />

### Windows 10

The simplest way to build native windows binaries seems to be using msys2.
Follow all steps of the [msys2](https://www.msys2.org/) installation guide,
then install the dependencies in the 'MSYS2 MinGW 64-bit` terminal:

``` shell
pacman -S mingw-w64-x86_64-fftw
pacman -S mingw-w64-x86_64-libtiff
pacman -S mingw-w64-x86_64-msmpi
pacman -S mingw-w64-cross-winpthreads-git
pacman -S git
```

Then get deconwolf and build:

``` shell
git clone https://www.github.com/elgw/deconwolf
cd deconwolf
make WINDOWS=1 -B
```

binaries will end up in the `bin` sub folder. All dependencies (DLLs) can be
found in `/mingw64/bin/` under the directory where msys64 is installed. To use
deconwolf outside of msys2 you will need to copy the binaries and the following
DLL files to the same folder.

```
libdeflate.dll libfftw3f-3.dll     libgcc_s_seh-1.dll libgomp-1.dll
libgsl-25.dll  libgslcblas-0.dll   libjbig-0.dll      libjpeg-8.dll
libLerc.dll    liblzma-5.dll       libstdc++-6.dll    libtiff-5.dll
libwebp-7.dll  libwinpthread-1.dll libzstd.dll        zlib1.dll
```

At least one person has build deconwolf using Windows Subsystem for Linux but
beware, there might be some [performance penalty](https://www.phoronix.com/scan.php?page=article&item=wsl-wsl2-tr3970x&num=1).


### FreeBSD

The following packages were required:
``` shell
pkg install git
pkg install gmake
pkg install fftw3
pkg install tiff
pkg install gsl
pkg install sudo
```

To build and install deconwolf:
``` shell
gmake -f makefile-freebsd
sudo gmake install
```


<a name="use" />

## Minimal usage example
To generate a parametric PSF and deconvolve an image, all you need is something
like this:
``` shell
# generate PSF.tif
dw_bw --resxy 130 --resz 250 --NA 1.46 --ni 1.518 --lambda 460 PSF.tiff
# Deconvolve image.tiff -> dw_image.tiff
dw --iter 50 image.tiff PSF.tiff
```
For available options, please see

``` shell
dw --help
dw_bw --help
```

At the moment the documentation is limited, but there is a short
[usage guide](USAGE.md).

### Bugs

Most likely there are bugs and they can only be fixed when they are known.
Please open a [new ticket](https://github.com/elgw/deconwolf/issues) if you
have any issues with the program.

## References

The deconvolution algorithm is based on these papers:

 * Richardson, William Hadley (1972). "Bayesian-Based Iterative Method of Image
   Restoration". JOSA. 62 (1): 55–59.
   [doi](https://doi.org/10.1364/JOSA.62.000055)
 * Lucy, L. B. (1974). "An iterative technique for the rectification of observed
   distributions". Astronomical Journal. 79 (6): 745–754.
   [doi](https://doi.org/10.1086%2F111605)
 * Biggs, D.S.C. “Acceleration of Iterative Image Restoration Algorithms.”
   Applied Optics. Vol. 36. Number 8, 1997, pp. 1766–1775.
   [doi](https://doi.org/10.1364/AO.36.001766)
 * Biggs, D.S.C “Accelerated iterative blind deconvolution”. PhD thesis.
   University of Auckland, New Zealand, 1998.
 * M. Bertero and P. Boccacci, A simple method for the reduction of boundary
   effects in the Richardson-Lucy approach to image deconvolution,
   A&A 437, 369-374 (2005).
   [doi](https://doi.org/10.1051/0004-6361:20052717)
 * Lee, Ji-Yeon & Lee, Nam-Yong. (2014). Cause Analysis and Removal of Boundary
   Artifacts in Image Deconvolution. Journal of Korea Multimedia Society. 17.
   838-848.
   [doi](https://doi.org/10.9717/kmms.2014.17.7.838).

The PSF generation is based on these:

 * Max Born. Principles of optics : electromagnetic theory of propagation, interference,
   and diffraction of light. Cambridge: Cambridge University Press, 2019.
   ISBN: 978-1-108-47743-7.

 * F. Aguet. “Super-Resolution Fluorescence Microscopy Based on Physical
   Models”. EPFL Thesis no. 4418 (2009), 209 p. Swiss Federal Institute of
   Technology Lausanne (EPFL), May 2009
   [url](http://bigwww.epfl.ch/publications/aguet0903.html)

 * The `--li` option uses:
   Jizhou Li, Feng Xue, and Thierry Blu. “Fast and accurate three-dimensional
   point spread function computation for fluorescence microscopy”. In: Journal
   of the Optical Society of America A 34.6 (May 2017), p. 1029.
   [doi](https://doi.org/10.1364/josaa.34001029)


## Alternatives
This is a non-complete list of alternative deconvolution software:

Free and open source:
 - [Deconvolution Lab2](http://bigwww.epfl.ch/deconvolution/deconvolutionlab2/)

Commercial:
 - [Huygens](https://svi.nl/HomePage)
 - [Microvolution](https://www.microvolution.com/)
