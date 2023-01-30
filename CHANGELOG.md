# CHANGELOG

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
