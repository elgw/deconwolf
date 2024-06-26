# Installation notes

## General requirements on Linux/Unix

- A compiler, tested with [gcc](https://gcc.gnu.org/) and
  [clang](https://clang.llvm.org/).
- openmp, typically bundled with the compiler.
- [cmake](https://cmake.org/)
- A package manager to get the dependencies.
- [https://www.gnu.org/software/make/](make)
- [fftw3](https://www.fftw.org/)
- [libtiff](https://libtiff.gitlab.io/libtiff/)
- [GSL - GNU Scientific Library](https://www.gnu.org/software/gsl/)
- [libpng](http://libpng.org/pub/png/libpng.html)

## Normal build procedure with CMake

To build:

``` shell
mkdir builddir
cd builddir
cmake ..
cmake --build .
```

### Without GPU acceleration:

By default **dw** will be compiled with OpenCL/GPU support enabled. If
 that does not work for you, turn that option off by passing
 `-DENABLE_GPU=OFF` to **cmake**, i.e. use:

 ``` shell
mkdir builddir
cd builddir
cmake -DENABLE_GPU=OFF ..
cmake --build .
 ```

### Machine specific optimizations

To enable optimizations for your particular machine, pass the option
 `-DENABLE_NATIVE_OPTIMIZATION=ON`. Doing that might create binaries
 that run faster on your machine. However it might not be possible to
 run the same binaries on another machine.


### To install

A direct install (files will be copied to default paths), use:
``` shell
sudo make install
```

To use the system package manager is probably the better way since it
offers a cleaner way to uninstall:

``` shell
mkdir builddir
cd builddir
cmake ..
cpack -G DEB # see cpack --help for other types
sudo apt-get install ./deconwolf-0.3.8-Linux.deb # on Ubuntu
```


## GPU Acceleration via OpenCL

- If you have more than one OpenCL device, please use the `--cldevice
  n` argument to tell deconwolf which device to use. The devices are
  ordered with 0 being the first. If available use the command
  `clinfo` to list the OpenCL devices available on your system.

- Under WSL2 and with the current version of pocl it does not seem to
work, see [issue
#56](https://github.com/elgw/deconwolf/issues/56). Hopefully it can
work in the future.

- If cmake can't find OpenCL and you are on a cluster, check out:
  [https://github.com/elgw/deconwolf/issues/55#issuecomment-2137579066]

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

## macOS Sonoma, Apple Silicon

For building you will need XCode from the App Store and [brew](https://brew.sh/).

Then set up XCode and install the required packages:
``` shell
xcode-select --install
brew install cmake
brew install gcc
brew install pkg-config
brew install libomp
brew install libtiff
brew install fftw
brew install gsl
```

To build, make sure you use the Homebrew-installed version of gcc, in my case gcc-14:

```shell
mkdir builddir
CC=gcc-14 CXX=g++-14 cmake ..
CC=gcc-14 CXX=g++-14 cmake --build .
make install
```

## Windows 10/11
Deconwolf can be built several different ways under Windows:
1. Using WSL, then follow the instructions for Ubuntu. Most likely there will be a
[performance
penalty](https://www.phoronix.com/scan.php?page=article&item=wsl-wsl2-tr3970x&num=1),
and it will not be possible to enable GPU acceleration, see [issue
#56](https://github.com/elgw/deconwolf/issues/56).
2.  [msys2](https://www.msys2.org/) or
[cygwin](https://www.cygwin.com/), however those options will be
slower since OpenMP will be using an pthreads emulation on top of
windows threads. It might be possible to get OpenCL working.
3. As a native windows program. This is the preferred way since it should work with OpenCL without any problems..

To build native windows programs, at least the following software is
required:

- [git](https://git-scm.com/download)
- [cmake](https://cmake.org/download/)
- Visual studio with clang.
- [vcpkg](https://github.com/microsoft/vcpkg?tab=readme-ov-file#quick-start-windows)

The dependencies can get retrieved by vcpkg:
``` shell
git clone https://github.com/microsoft/vcpkg
.\vcpkg\bootstrap-vcpkg.bat
.\vcpkg\vcpkg.exe install fftw3[threads]
.\vcpkg\vcpkg.exe install tiff
.\vcpkg\vcpkg.exe install gsl
.\vcpkg\vcpkg.exe install opencl
.\vcpkg\vcpkg.exe install getopt
.\vcpkg\vcpkg integrate install
```
Please note the value of the `CMAKE_TOOLCHAIN_FILE` as it will be used later.

A visual studio project can be created by

``` shell
cd deconwolf
mkdir build
cd build
cmake "-DCMAKE_TOOLCHAIN_FILE=C:/YOUR/OWN/PATH/vcpkg/scripts/buildsystems/vcpkg.cmake" -T ClangCL -A x64 ../
```
Important: Please use the correct path to `vcpkg.cmake`.

Open the visual studio solution and oo to Linker->Input->Additional
dependencies and add:

``` shell
libomp.lib
```

Change the build type from debug to release and compile.

If you are a windows developer and reading this, please help us out
to make the build process smoother!

## FreeBSD
- Use `gmake`, not `make`.

Packages:
``` shell
pkg install git
pkg install gmake
pkg install fftw3
pkg install tiff
pkg install gsl
pkg install sudo
```

## CentOS
Tested on CentOS Linux Release 7.8.2009 (Core).

Dependencies:
``` shell
sudo yum install gcc gsl-devel libtiff-devel fftw-devel
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

## Ubuntu 22.04

First install the required packages:

``` shell
sudo apt-get update
sudo apt-get install \
 cmake               \
 pkg-config          \
 gcc                 \
 libfftw3-single3    \
 libfftw3-dev        \
 libgsl-dev          \
 libomp-dev          \
 libpng-dev          \
 libtiff-dev
```


## Ubuntu 23.04
Same as Ubuntu 22.04. Possibly also

``` shell
apt-get install opencl-headers
```

## Arch/ Manjaro

``` shell
# remember to update system
sudo pacman -Suuyy
# install dependencies
sudo pacman -S fftw gsl openmp libtiff
```

## Rapsberry PI (64-bit Debian bookworm)
``` shell
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install libfftw3-dev \
libtiff-dev \
libgsl-dev
```

## MKL
FFTW3 is the default FFT backend for deconwolf but it is also possible to use
Intel MKL. At some point it was possible to choose MKL via

``` shel
sudo apt install intel-mkl
make MKL=1 -B
```

To set the number of threads, set the environmental variable
`MKL_NUM_THREADS`, for example:
``` shell
export MKL_NUM_THREADS=8
dw ...
```
If you are interested in using MKL please open a new issue.
