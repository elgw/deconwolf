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
: Explicitly set the name of the output file. By default the name of the
input file prefixed by `dw_` will be used.

**\--prefix str**
: Set the prefix to use for the output file. An extra `_` will be appended
to the str.

**\--iterdump**
: Save each iteration to disk, not only the final.

**\--niterdump n**
: Save every nth iteration to disk, including the final image.

**\--iter n**
: Set the number of iterations.

**\--tilesize s**
: Set the size (axial side length, in pixels) of the largest portion that
can be deconvolved at a time. E.g., if s is 2048 any image larger than 2048
will be deconvolved in tiles.

**\--tilepad p**
: Set how many pixels the tiles should overlap.

**\--overwrite**
: Overwrite the target if it exists.

**\--bq n**
: Border/edge handling. n=0: no padding. n=1: a compromise between memory/speed
and quality. n=2 default.

**\--float**
: Save the output using 32-bit floats, also turns off any scaling.

**\--version**, **-V**
: Show version information and quit.

**\--help**
: Show a short help message and quit.

**\--method METHOD**
: Select what method to use. Valid options are
 - `id` identity transform, i.e. nothing. Useful to see if images
   loads, scales and saves correctly.
 - `rl` classical Richardson-Lucy.
 - `ave` Additive Vector Extrapolation by Biggs and Andrews.
 - `eve` Exponential Vector Extrapolation by Biggs.

**\--psigma s**
: Pre-process the image and the PSF by convolving it by a 3-D isotropic
Gaussian with sigma=s. This acts as a low pass filter.
Have been found useful for very noisy image.

**maxproj**
: With *maxproj* as the first argument deconwolf will create max
projections of all following tif files. Output will be prefixed with `max_`.

# While running
At normal verbosity deconwolf will put one green dot per FFT. The fMSE value
reported at each iteration is the Mean Squared Error between the input image
and the forward projected current guess (i.e. the PSF convoved with the
current guess).

# OUTPUT
Without specifying the output file name, the output file will
be prefixed with `dw_`, e.g, if you deconvolve `file.tif`
the output file will be called `dw_file.tif`. Please pay attention
to the log file that also will be created (as `dw_file.tif.log.txt`).

Unless **\--float** was specified the output images will be scaled
to use the full dynamic range of the 16 bit format. The scaling factor
can be found in the log file.

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
