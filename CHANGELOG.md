# CHANGELOG

### 2020-11-10 v0.0.9
 * Updated the PSF generator to use Simpson's rule with 9x9 valuations in X-Y per pixel. This should produce more accurate PSF's when the pixel size is large. It is future work to also integrate in Z.
 * Fixed some bugs/issues related to vector aliasing that caused dw to crash when compiled with the `O3` flag.
 * General cleanup including warning free compilation.

### 2020-06-09 v0.0.6
 * Added the binary `dw_bw` which generates PSFs according to the Born-Wolf model. It is based on, and compared to [https://github.com/Biomedical-Imaging-Group/PSFGenerator/blob/master/src/psf/bornwolf/KirchhoffDiffractionSimpson.java](https://github.com/Biomedical-Imaging-Group/PSFGenerator/blob/master/src/psf/bornwolf/KirchhoffDiffractionSimpson.java)

### 2020-06-05 v0.0.5
 * Changing name of the main binary to `dw` and prefixing the utility binaries with `dw_`.
 * Does not load the full images neither for reading or writing in tiling mode.
 * Added the binary `dw_tiffmax` for fast maximum intensity projections.

### 2020-05-17 v0.0.3
 * Further memory savings
 * Everything in the main loop is parallelized with OpenMP
 * Automatic cropping of redundant layers of the PSF, controlled by the `--xyfactor` setting.
 * Less clutter in the output to the terminal.

### 2020-05-13 v0.0.2
 * Memory savings of about 20%
 * Refactorization
 * Added the `--relax` flag in order to scale the central value of the PSF.

### 2020-05-12 v0.0.1
 * First release!
