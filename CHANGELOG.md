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
