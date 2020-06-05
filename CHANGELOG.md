# CHANGELOG

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
