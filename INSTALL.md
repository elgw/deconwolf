# Installation notes

 - [Meson](#Meson)
 - [Windows 10](#Windows10)
 - [macOS Big Sur](#macOS-Big-Sur)
 - [Free BSD](#freebsd)
 - [CentOS](#CentOS)
 - [Ubuntu 16.04](#ubuntu-1604)
 - [FreeBSD](#FreeBSD)
 - [MKL FFT Backend](#MKL)


## Meson
deconwolf can also be installed using [meson](https://mesonbuild.com/),
tested only Ubuntu and MacOS

To build and install:
``` shell
meson setup builddir --buildtype release
cd builddir
meson compile
meson install # Note: only works if meson was installed by sudo
```

To uninstall:

``` shell
sudo ninja -C builddir uninstall
```
## macOS Big Sur

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

## Windows 10
Although it might be possible to build native windows executables, it
is suggested that WSL is used.

With the following instructions it was possible to build deconwolf
0.1.0 on Windows 10 (the --inplace option might not work on the
current version).

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
beware, there might be a
[performance penalty](https://www.phoronix.com/scan.php?page=article&item=wsl-wsl2-tr3970x&num=1).



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
sudo yum install python3
sudo python3 -m pip install meson ninja
# Then follow the general meson instructions in README.md
```

## Ubuntu 16.04
``` shell
sudo apt-get update
sudo apt-get install gcc
sudo apt-get install pkg-config
sudo apt-get install libfftw3-single3
sudo apt-get install libfftw3-dev
sudo apt-get install openmp
sudo apt-get install libtiff-dev # only difference to 20.04
sudo apt-get install libgsl-dev
sudo apt-get install libomp-dev
sudo apt-get install libpng-dev

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

## Arch/ Manjaro

``` shell
# remember to update system
sudo pacman -Suuyy
# install dependencies
sudo pacman -S fftw, gsl, openmp, libtiff
make
sudo make install

```
