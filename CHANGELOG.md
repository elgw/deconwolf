# CHANGELOG

## 0.3.8
- For systems with multiple GPUs or OpenCL compatible devices it is
  now possible to select which to use with **--cldevice**. To figure
  out which that are available it is simplest to use
  `clinfo`. Alternatively call dw with **--verbose 2** or above.

- Removed depreciated makefile for freebsd as it is no longer needed.

- Removed anything related to meson as it is no longer needed.

- Using the `PRI*` macros from `inttype.h`,
  especially `PRIu64` for `uint64_t` and `PRId64` for `int64_t` to get
  rid of some warnings under MacOS.

- Minor changes.

## 0.3.7
- Deconwolf compiles as a native windows program using clang. So far
  the binaries are only smoke tested since the main target is linux.

- **dw_bw** use OpenMP and does not rely on pthreads any more (for
  portability reasons).

- Removed the AVE and EVE methods since they don't add anything over SHB.

- Added a `CMakeList.txt` for building with cmake.

- Added `--gpu` which at the moment is equivalent to `--method shbcl2`
  but a little more mnemonic.

- Added the `--periodic` option which turns on periodic boundary conditions,
  i.e. is equivalent to `--bq 0`.

## 0.3.6
- The GPU code path uses in-place transformations as much as possible
  to save a little on the memory usage.

- Switched to [VkFFT](https://github.com/DTolm/VkFFT) (v1.3.3) as the
  default FFT backend on the GPU. Unless a big regression is found,
  the clFFT code path will most likely not be maintained in future
  versions and be removed.

  To build with GPU acceleration use:

  ```
  make kernels
  make -B VKFFT=1
  ```

  As before, you need also to choose `--method shbcl2` to use it over
  the CPU implementation.

  Initial tests show a speed up of about 10-30% depending on the image
  size. As a bonus VkFFT will process any sizes while clFFT simply
  refuse to process the tricky ones.

- Identified that `cl_idiv_kernel.c` took a substantial amount of the
  iteration time and rewrote it.

- Removed the "CUDA" backend since it does not make sense any more.

- Checks that the min value of the image > 0. Aborts if not.

- Checks that the max value of the image >= 1. Aborts if not.

## 0.3.5
- **dw maxproj** There were problems reading the output in
  MATLAB. Updated so that the output image will be written as a single
  strip.

## 0.3.4
- Minor bug fixes which gives a clean build with `-fanalyzer`.

## 0.3.3
- Writing pixel size to output file also when tiling is used.

## 0.3.2
- Tested on raspberry pi 4 using 64-bit bookworm.
- Found a bug in `fft.c` where `memcpy` was used wrongly (replaced by
  `memmove`). Strangely that bug never manifested under
  Ubuntu/x86_64.
- Added **fim_realloc** for aligned reallocs. This function could be
  branched depending on the OS since there are platform specific
  aligned reallocation functions.
- Header files: Using `#pragma once` instead of the `#ifndef file_h_`
  pattern.

## 0.3.1
- Introduced **fim_malloc** for all allocations that might benefit
  from a stricter alignment than malloc provides by default. Tested
  with `MADV_HUGEPAGE` for the allocations but the results are
  inconclusive (but it uses more RAM when enabled). Cleared all uses
  of `fftw_free` and `fftw_malloc`.

## 0.3.0
- Respects the NO_COLOR environmental variable in accord with https://no-color.org/.
- Fixed correct capping of pixel values when **--scaling** is used.

## v 0.2.9
 - Added the command line option **--scaling** for setting bypassing the
   automatic image scaling in 16-bit output mode.

## v 0.2.8
 - Switched from `fftw3f_threads` to `fftw3f_omp`. This reduced the
   run time by about 10% on a Intel i7-6700K. Can be reverted by
   commenting in/out the corresponding lines in the makefile.
 - Cleaned up the output of `dw --version`

## v 0.2.7
- Converted a few minor code paths to execute in parallel by OpenMP
  directives.

## v 0.2.6
- Using ISO 8601 in log files, e.g., `2023-02-14T11:14:14`.

## v 0.2.5
- Added the **--xyz** option to **dw maxproj**, for creating max
  projections along the three axes and collecting them on a single 2D
  image.

## v 0.2.4
- **dw --help** now shows the additional commands/modules available.
- Reading 16-bit tif files with **TIFFReadEncodedStrip** instead of
  **TIFFReadRawStrip**. Some programs saves tiff files in other ways
  :)
- Added the command psf-STED for 3D STED PSFs. Use at your own risk.
- Building with meson is temporarily broken and to be fixed.
- Fixed dw chashing when combining --method rl with --iterdump
- Setting the background level automatically to min(image) unless
  specified with **--bg**.

## v 0.2.3
- Fixed some errors introduced in v 0.2.2, especially the **dw
  maxproj** was broken.
- added the subcommand **dw merge**. To be used to merge single z-planes
  into a 3D volume.

## v 0.2.2
- Can deconvolve using clFFT, when compiled with **OPENCL=1** two new
  methods appear, **--method shbcl** and **--shbcl2**, the first using
  clFFT only for the Fourier transforms, the latter using OpenCL for
  the whole deconvolution procedure. Uses quite much GPU memory which
  is something to improve upon in future version, possibly by
  switching to vkFFT.

## v 0.1.1
- Added experimental **dw imshift** for shifting images, also shift
  estimation using normalized cross correlation with **dw imshift
  --ref file.tif**. Might be extended to basic tiling etc.

## v 0.1.0
- Implements the 'Scaled Heavy Ball'. More
memory efficient than eve and about the same speed and image
quality. Might become the default method.
- Reorganization of code with one file per deconvolution method, RL is
  now separated to an own file which improves readability.
- The **--method** argument can be used to switch between several
methods, see **--help** or the man page.
- Showing Idivergence after each iteration, switch back to MSE with
**--mse**
- Cleaned up the text written to the terminal, notably any warnings
from libtiff now go to the log file.
- OMP is set to use as many cores as FFTW.
- Added OMP directives to a few more loops.
- Using static OMP schedule.
- Introduces the **--tsv** argument to save information per
iteration to a separate tsv file for easier plotting and analysis.
- Three different stopping criteria: Relative error (default) Fixed
  number of iterations or at an absolute error.

## v 0.0.26
- **dw maxproj** works with file that are not in the current folder.
- Fixed **--iterdump** not always working.

## v. 0.0.25
- Builds with cuFFT on Linux, use `make CUFFT=1 -B`, requires a CUDA
compatible GPU and of course the cuFFT library installed.

## v. 0.0.24
- Tested on CentOS, install both with make and meson.
- Fixed a memory leak with the **--tilesize** option causing
crashed sometimes.

## v. 0.0.23
- Added 'meson.build' files in order for deconwolf to be built by
[The Meson Build system](https://mesonbuild.com/), tested to work
on both Ubuntu 21.10 and MacOS (on x86_64 hardware).
- Added a small test image under `demo/` together with a **makefile**
to deconvolve it.
- Added [pseudo code](PSEUDOCODE.md) for the binaries hoping to
planning to replace this by a properly typeset and more detailed
document.
- Aborting if the number of threads is set < 1.
- The algorithm is still unchanged since v 0.0.20.

## v. 0.0.22
- Fixed double free-bug in tiling mode.

## v. 0.0.21
- Updated documentation and man-pages based on markdown files
for easier updating.
- Provides `makefile-freebsd` for building on FreeBSD 13.0
- Changed behavior when too few input arguments are given to
only give a two-line message.

## v. 0.0.20
- Changing acceleration technique to use
'Exponential Vector Extrapolation' (EVE) described in Biggs PhD thesis.
Deconvolved images get higher MSE but much lower I-div.
- '--xyfactor 0' does not crash dw anymore.
- Frees the PSF as soon as not needed to save some memory.
- Changing the behavior of the progress dots to appear more linear
in time
- Changing the non-negative condition to strictly positive in order for
pixel not to get stuck at 0.
- Adding the option to turn off Biggs acceleration, i.e. run normal
Richardson-Lucy with --biggs 0.
- Will load PSFs that don't have an odd number of pixels in each dimension
however that is not recommended.
- Can be built against Intel MKL (`make MKL=1 ...`), consider that an
experimental option. 14 percent faster on a small test image, varied
results on larger images.

## v. 0.0.19
- Using lanczos5 instead of lanczos3 for the PSF generation. As a result
GSL_EROUND is not raised for the test cases.
- Faster PSF generation, using more symmetries.
- dw_bw can now use more than one thread (wrongly disabled in v 0.0.18).

## v. 0.0.18
- Provided install instructions for Windows 10.
- Fixed some mismatching fftwf_malloc/fftwf_free where they were
mixed up with malloc/free causing crashes on Windows.
- Added an experimental src/CMakeLists.txt that can be used when
building with cmake. It is also possible to cross compile for Windows
on Linux although it takes some effort to collect the DLL files for the
dependencies.

## v. 0.0.17
- Fixed some bugs in the PSF generation code that did affect the accuracy
of the pixels in the PSF.
- Stared to use GSL for numerical integration. It remains to change the
double integral over x-y into something more dynamic.
