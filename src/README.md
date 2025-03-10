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
 * `cl_util.c/h` GPU abstraction code.

Implemented methods:
 * `method_rl.c/h` Standard Richardson Lucy
 * `method_shb.c/h` RL accelerated by the Scaled Heavy Ball algorithm
 * `method_shb_cl.c/h` as above, using clFFT for the FFTs.
 * `method_shb_cl2.c/h` as above, but with as much as possible of the
   calculations on the GPU. If there is enough memory this should be
   the fastest alternative. (TODO).

Extra stuff under development:
 * `dw_dots.c/h` for dot detection.
 * `dw_maxproj.c/h` for max projections of tif files.

Tiling:

When the tiling option is used the following functions are relevant:
- `fim_tiff_to_raw` -- converts a tif (u16 or f32) file to a contiguous raw (f32) data file
- `tiling_get_tile` -- read a tile from a raw file
- `tiling_put_tile` -- write a tile to a raw file
- `fim_tiff_from_raw` -- convert a raw file to tif (f32 -> u16)

Status:
 - No support for `npy` files.
