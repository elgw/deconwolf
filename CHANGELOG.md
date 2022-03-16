# CHANGELOG

- v 0.1.0
  - Reorganization of code with one file per method.
  - Separate RL implementation for clarity.
  - The **--method** argument can be used to switch between several
    methods, see **--help** or the man page.
  - Implements the 'Scaled Heavy Ball'. More
    memory efficient than eve and about the same speed and image
    quality. Might become the default method.
  - OMP is set to use as many cores as FFTW.
  - Introduces the **--tsv** argument to save information per
    iteration to a separate tsv file.

- v 0.0.26
 - **dw maxproj** works with file that are not in the current folder.
 - Fixed **--iterdump** not always working.

- v. 0.0.25
  - Builds with cuFFT on Linux, use `make CUFFT=1 -B`, requires a CUDA
  compatible GPU and of course the cuFFT library installed.

- v. 0.0.24
  - Tested on CentOS, install both with make and meson.
  - Fixed a memory leak with the **--tilesize** option causing
    crashed sometimes.

- v. 0.0.23
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

- v. 0.0.22
  - Fixed double free-bug in tiling mode.

- v. 0.0.21
   - Updated documentation and man-pages based on markdown files
     for easier updating.
   - Provides `makefile-freebsd` for building on FreeBSD 13.0
   - Changed behavior when too few input arguments are given to
     only give a two-line message.

- v. 0.0.20
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

- v. 0.0.19
    - Using lanczos5 instead of lanczos3 for the PSF generation. As a result
      GSL_EROUND is not raised for the test cases.
    - Faster PSF generation, using more symmetries.
    - dw_bw can now use more than one thread (wrongly disabled in v 0.0.18).

- v. 0.0.18
    - Provided install instructions for Windows 10.
    - Fixed some mismatching fftwf_malloc/fftwf_free where they were
      mixed up with malloc/free causing crashes on Windows.
    - Added an experimental src/CMakeLists.txt that can be used when
      building with cmake. It is also possible to cross compile for Windows
      on Linux although it takes some effort to collect the DLL files for the
      dependencies.

- v. 0.0.17
   - Fixed some bugs in the PSF generation code that did affect the accuracy
     of the pixels in the PSF.
   - Stared to use GSL for numerical integration. It remains to change the
     double integral over x-y into something more dynamic.
