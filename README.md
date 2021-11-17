# deconwolf

 * [Build and Install](#install)
   * [Linux](#linux)
   * [macOS](#osx)
   * [Windows 10](#win10)
 * [References](#ref)

`deconwolf` is a program for deconvolution of 3-D fluorescent wide-field
image stacks:
 - The deconvolved images shows very mild boundary effects which means that you
   can crop and deconvolve small regions of interest.
 - RAM usage can be reduced at the cost of slightly longer computation times by
   tiling. That makes it possible to deconvolve large images on small machines.
 - It can make use of all precious cores of your "big" machine since the
   critical parts run on separate threads (as many as you would like).
 - Deconwolf is tiny! The binaries could even fit on a floppy drive (if you are fortunate
   enough to own one of those antiquities).

Except for this README there is also a short [manual](USAGE.md),
a [CHANGELOG](CHANGELOG.md) and a [TODO](TODO.md) list.

Does does not:
 - Have a full-featured Graphical User Interface (GUI), however, if you
 like to click buttons there is a tiny [GUI](https://github.com/elgw/dw_gui).
 - Show your images, for that, try [ImageJ](https://imagej.net/Welcome)
 - Diagnose your imaging system.
 - Estimate your PSF based on real images.

At the moment no pre-built packages are maintained so you will have to build
deconwolf from the source code, instructions follows below.

<a name="install"/>

## Build and install
Deconwolf runs on 64-bit machines with x86_64 architecture, it has been built
and installed on Ubuntu 20.04 and macOS Big Sur and Windows 10.
Instruction for other systems will be collected [here](INSTALL.md). Go directly to [Windows 10](#win10), [Linux](#linux) or [macOS](#osx).

<a name="deps">

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

<a name="linux"/>

### Ubuntu 20.04
Ensure that the required packages are installed

``` shell
sudo apt-get update
# find out actual names with command like
# sudo apt-cache search fftw
sudo apt-get install libfftw3-single3
sudo apt-get install libfftw3-dev
sudo apt-get install openmp
sudo apt-get install tiff-5
sudo apt-get install libgsl-dev
```

To build and install:
``` shell
make -B
sudo make install
```

<a name="osx"/>

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

<a name="win10"/>

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

<a name="ref">

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
