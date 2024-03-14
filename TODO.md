# To do

This document is complementary to the [issues page on
github](https://github.com/elgw/deconwolf/issues).

## Planned

- [ ] c2r and c2r plans needed any more? Probably only the in-place versions.

- [ ] dynamic loading of OpenCL, for examples using [https://github.com/yugr/Implib.so]

- [ ] Complete switch to Vulkan from OpenCL. Doing this will ensure
      better portability, for example on new macs with
      [MoltenVK](https://github.com/KhronosGroup/MoltenVK). The switch
      from clFFT to VkFFT in v 0.3.6 was a first step towards this goal.

- [ ] Crash-safe writing of output images. I.e, follow the standard
  procedure and first, write to temporary file and move to final
  destination when the write has finished. Avoids some bad luck.

- [ ] Add OpenCL-RL method (as a low mem option).

- [ ] Get documentation up to date and write a general introduction.

- [ ] Set default PSF size in `dw_bw` based on the geometry of the PSF
and expected sample thickness.

- [ ] Prepare binary packages, should include Windows 10/11, MacOS and
      linux.

- [ ] Crop the PSF also in the axial direction when it vanishes, this
      should be good for confocal, STED etc.

- [ ] `--xyfactor` is a quite awkward name.

- [ ] Still one kernel to write in order to speedup the initialization
      in `method_shbcl2` when bq > 0.

- [ ] Wrap image pointers in some abstraction, would be especially
      useful for Windows where there are multiple **free**-functions.

- [ ] Make the **--lookahead** option work again, and possibly enable
      it by default, at least on the CPU side.

## Nice to have
- [ ] Skip libtiff, especially for writing.
- [ ] See if Pinned memory allocations can improve the performance for shbcl (don't expect anything drastic from shbcl2).
[LFAT](https://ltfat.github.io/notes/ltfatnote017.pdf),
[code](https://github.com/ltfat/ltfat/blob/master/fourier/nextfastfft.m)
- [ ] Wrap malloc to catch allocation errors to the log file?
- [ ] Decide about default tiling settings based on the PSF.
- [ ] Demos, for example on the effect of the tiling.
- [ ] Implement the Gibson-Lanni PSF model.
- [ ] Double check the tiling weights using gradients and other patterns.
- [ ] Add information on the input image to the log -- especially
detect saturated pixels.
- [ ] Fail gracefully when a pyramidal tif is supplied.
- [ ] Fail gracefully when a multi color/time/... tif is supplied.
- [ ] Include CCD corrections.
- [ ] Make the background correction work also in tiled mode.
- [ ] Automatic image cropping in the axial direction based on
      gradient magnitude or similar.
- [ ] Do the padding/unpadding for in-place transformations on the GPU
      side. This should not affect the memory usage since it is only
      performed at the beginning and at the end of the algorithms.
