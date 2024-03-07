# Installation notes

## OpenCL support
To use GPU acceleration **dw** needs to be compiled with that option
enabled. Using the makefile that corresponds to:

``` shell
make kernels
make -B VKFFT=1
```

## Build systems
Under linux the makefile should work on most distributions. For those
that prefer something else there is also configuration files for Meson
and CMake.

### Meson
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

### CMake
To build:

``` shell
mkdir builddir
cd builddir
cmake ..
cmake --build .
```

## Platform specific notes
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

## Windows 11
Deconwolf can be built and run using WSL. Most likely there will be a [performance
penalty](https://www.phoronix.com/scan.php?page=article&item=wsl-wsl2-tr3970x&num=1), and it will not be possible to enable GPU acceleration.

It can also be built using [msys2](https://www.msys2.org/) or
[cygwin](https://www.cygwin.com/), however those options will be
slower since OpenMP will be using an pthreads emulation on top of
windows threads. It might be possible to get OpenCL working.

To build native windows programs, at least the following software is
needed:

- [git](https://git-scm.com/download)
- [cmake](https://cmake.org/download/)
- Visual studio with clang.
- [vcpkg](https://github.com/microsoft/vcpkg?tab=readme-ov-file#quick-start-windows)

The dependencies can get retrieved by vcpkg:
``` shell
git clone https://github.com/microsoft/vcpkg
.\vcpkg\bootstrap-vcpkg.bat
.\vcpkg integrate install
.\vcpkg\vcpkg.exe install fftw3
.\vcpkg\vcpkg.exe install tiff
.\vcpkg\vcpkg.exe install gsl
.\vcpkg\vcpkg.exe install opencl
```

A visual studio project can be created by

``` shell
cd deconwolf
mkdir winbuild
cd winbuild
cmake .. -G "Visual Studio 17 2022" -T ClangCL -A x64
```

### FreeBSD
- Use `gmake`, not `make`.
- Default compiler: `clang`
- `pkgconf` not `pkg-config`.

Packages:
``` shell
pkg install git
pkg install gmake
pkg install fftw3
pkg install tiff
pkg install gsl
pkg install sudo
```

### CentOS
Tested on CentOS Linux Release 7.8.2009 (Core).

``` shell
# Install dependencies
sudo yum install gcc gsl-devel libtiff-devel fftw-devel
```

### Ubuntu 16.04
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

### Ubuntu 23.04
Same as Ubuntu 22.04.

``` shell
apt-get install opencl-headers
```

### Arch/ Manjaro

``` shell
# remember to update system
sudo pacman -Suuyy
# install dependencies
sudo pacman -S fftw gsl openmp libtiff
make
sudo make install
```

### Rapsberry PI (64-bit Debian bookworm)
``` shell
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install libfftw3-dev \
libtiff-dev \
libgsl-dev
```

## Library options

### MKL
FFTW3 is the default FFT backend for deconwolf but it is also possible to use
Intel MKL. This option is only tested on Ubuntu so far.

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
