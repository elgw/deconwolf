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
