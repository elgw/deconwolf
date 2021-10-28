This `src/`-folder contains:

 * `deconwolf.c` Entry point for `dw`
 * `deconwolf_tif_max.c/h` Code to do max-projections over z in tif images.
 * `dw.c/h` Main functionality
 * `dw_bw.c/h` The PSF generator, compiles to `dw_bw`
 * `dw_version.h` contains the version number, is updated manually.
 * `fim.c/h` Manipulation of 'float' images.
 * `fim_tiff.c/h` Read/write tiff files.
 * `makefile` generates unit tests for various components.
 * `tiling.c/h` Process images in tiles
 * `fft.c/h` Some FFT operations on the images.
