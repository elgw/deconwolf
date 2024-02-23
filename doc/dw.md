% DECONWOLF(1) dw ABCDE
% Erik Wernersson
% 2022

# NAME
deconwolf **dw** is a tool for deconvolution of wide field microscopy stacks.

# SYNOPSIS
For deconvolution:

**dw** [*OPTIONS*] file.tif psf.tif

or, for max projections over z:

**dw** maxproj file1.tif file1.tif ...

# OPTIONS
**\--threads t**
: Set the number of computational threads to use.

**\--out file.tif**
: Explicitly set the name of the output file. By default the output
file name is the name of the input file prefixed by `dw_`.

**\--prefix str**
: Set the prefix to use for the output file. An extra `_` will be appended
to the str.

**\--iterdump**
: Save each iteration to disk, not only the final.

**\--niterdump n**
: Save every nth iteration to disk, including the final image.

**\--iter n**
: Set the number of fixed iterations.

**\--relerror r**
: Tell **dw** to stop when the relative error between to iterations is
  below **r**. This can be combined with **\--maxiter** to set a limit
  on the number of iterations to perform.

**\--abserror e**
: Tell **dw** to stop when an absolute error is reached. This can be
  combined with **\--maxiter** to set a limit on the number of
  iterations to perform.

**\--maxiter m**
: If a relative or abolute error is used as stopping criteria, this
  limits the number of iterations to perform.

**\--tilesize s**
: Set the size (axial side length, in pixels) of the largest portion that
can be deconvolved at a time. E.g., if s is 2048 any image larger than 2048
will be deconvolved in tiles. See separate section on tiling below.

**\--tilepad p**
: Set how many pixels the tiles should overlap.  See separate section
  on tiling below.

**\--overwrite**
: Overwrite the target if it exists.

**\--bq n**
: Border/edge handling. n=0: no padding. n=1: a compromise between memory/speed
and quality. n=2 default.

**\--float**
: Save the output using 32-bit floats, also turns off any
scaling. Please note that without this flag, images will be saved as
16-bit unsigned integers and will be scaled to make most use of that
format, see the log files for scaling factor used.

**\--version**, **-V**
: Show version information and quit.

**\--help**
: Show a short help message and quit.

**\--method METHOD**
: Select what method to use. Valid options are
  i/ **id** identity transform, i.e. nothing. Useful to see if images
   loads, scales and saves correctly.
  ii/ **rl** classical Richardson-Lucy.
  iii/ **shb** Scaled Heavy Ball.
  iv/ **shbcl** Scaled Heavy Ball using OpenCl (not compiled by default)

**\--mse**
: Show the Mean Square Error between the input image and the current
  guess convolved with the PSF at each iteration. If not set,
  I-divergence is shown.

**\--psigma s**
: Pre-process the image and the PSF by convolving it by a 3-D isotropic
Gaussian with sigma=s. This acts as a low pass filter.
Have been found useful for very noisy image.

**\--noplan**
: disable FFTW3 planning. This means that FFTW3 uses the default plan
  for the given problem size.

**\--no-inplace**
: disable the use of in-place FFT transformations with fftw3. This
  will cause fftw3 to use more memory and is typically slower. This
  option is mostly for debugging. On some platforms it is expected
  that deconwolf will only work with this flag.

**maxproj**
: With *maxproj* as the first argument deconwolf will create max
projections of all following tif files. Output will be prefixed with `max_`.

# While running
At normal verbosity deconwolf will put one green dot per FFT. After
each iteration the Idiv or MSE (with **\--mse**) is shown, not that
this is not a measurement on the final image quality.  Warnings are
prefixed with ' ! '. More information can be found in the log file.

# INPUT
**dw** reads 16-bit (unsigned integers) or 32-bit (floating point) tif
files. It does not understand compressed files or any dimensions
except x,y,z, so only one color per image etc. To test if **dw** reads
and writes your files, try the **\-\-method id** option. In general, if
an image is saved by ImageJ dw should be able to read it.

# OUTPUT
Without specifying the output file name, the output file will
be prefixed with `dw_`, e.g, if you deconvolve `file.tif`
the output file will be called `dw_file.tif`. Please pay attention
to the log file that also will be created (as `dw_file.tif.log.txt`).

Unless **\--float** was specified the output images will be scaled
to use the full dynamic range of the 16 bit format. The scaling factor
can be found in the log file.

# PERFORMANCE
Without specifying any arguments to **dw** it will use one thread per
cpu core in the host machine. This is typically the fastest way to
deconvolve one image. For maximal throughput it is better to run
several copies of **dw** using fewer cores each. E.g. on an 8-core
machine the maximal throughput can be reached by deconvolving eight
images in parallel using one core each (if enough RAM is available).

# Tiling
In order to use less RAM and deconvolve really large scans deconwolf
can process images in a memory efficient way by dividing them into
smaller portions, or tiles. This process is completely transparent to
the user although the **--tilesize T** parameter has to be set
explicitly.

Tile processing would typically be slow and introduce artifacts at the
boundaries. To reduce the boundary artifacts tiles are overlapped by
up to **--tilepad p** pixels. This artifact remedy does,
unfortunately, slow down the processing even more.

Internally the tile processing performs the following steps:

 1. The tif input image is written to disk as raw float data, without
    loading the full image to RAM.
 2. A tiling grid is set up which divides the lateral domain of the
    image into tiles of size at most $TxT$.
 3. Each tile then is loaded from disk, including extra padding $p$
    where it isn't in contact with the edge. The tile is then
    deconvolved and the data is written to disk.
 4. Where the padding is overlapping another tile, the
    image data is weighted linearly to reduce artifacts.
 5. The raw output images is converted to tif, again without loading
    the full image to RAM.

Tiling is enabled only when **--tilesize** is specified.

# SEE ALSO
**dw_bw** for generation of point spread functions according to
the Born-Wolf model.

# WEB PAGE
[https://github.com/elgw/deconwolf/](https://github.com/elgw/deconwolf/)

# REPORTING BUGS
Please report bugs at
[https://github.com/elgw/deconwolf/issues/](https://github.com/elgw/deconwolf/issues/)

# COPYRIGHT
Copyright © 2022 Erik Wernersson.  License GPLv3+: GNU GPL version 3 or later
<https://gnu.org/licenses/gpl.html>.
This is free software: you are free to change and redistribute it.  There is NO WARRANTY, to the
extent permitted by law.
