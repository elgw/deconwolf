% DECONWOLF(1) dw_bw ABCDE
% Erik Wernersson
% 2022

# NAME
**dw_bw** is a tool to generate Point Spread Functions according to the
Born-Wolf model.

# SYNOPSIS
**dw_bw** [*OPTIONS*] psf.tif

# OPTIONS

**\--ni #**
: Set the refractive index of the immersion media.

**\--NA #**
: Set the numerical aperture of the objective.

**\--lambda #**
: Set the emission wavelength [nm].

**\--resxy #**
: Set the pixel size [nm] in the lateral plane.

**\--resz P**
: Set the pixel size/distance between planes in the axial direction.

**\--size N**
: Set output size to N x N x P [pixels].
This should be set big enough so that the PSF cone
isn't cropped at the first and last slice of the PSF.

**\--nslice P**
:  Set output size to N x N x P [pixels].
P has to be an odd number. Please not that P should be at
least (2Z-1) where Z is the number of slices in the image
to be deconvolved

**\--threads t**
: Set the number of computational threads to use.

**\--overwrite**
: Overwrite the target if it exists.

**\--li**
: Use Li's method (our implementation) for the BW integral.
About 2x faster for the default PSF size. Not thoroughly tested.

**\--version**, **-V**
: Show version information and quit.

**\--verbose L**
: Set the verbosity level to L where 0=quiet, 1=default, 2=full.

**\--help**
: Show a short help message and quit.

# OUTPUT
A 32-bit float tif file along with a log file.

# EXAMPLE

``` shell
dw_bw --resxy 65 --resz 200 --NA 1.45 --ni 1.515 --lambda 466 PSF_DAPI.tif
```
saves `PSF_DAPI.tif` and `PSF_DAPI.tif.log.txt`.

# SEE ALSO
**dw** for deconvolution.

# WEB PAGE
[https://github.com/elgw/deconwolf/](https://github.com/elgw/deconwolf/)

# REPORTING BUGS
Please report bugs at
[https://github.com/elgw/deconwolf/issues/](https://github.com/elgw/deconwolf/issues/)

# COPYRIGHT
Copyright © 2022 Erik Wernersson.  License GPLv3+: GNU GPL version 3 or later
<https://gnu.org/licenses/gpl.html>.
This is free software: you are free to change and redistribute it.
There is NO WARRANTY, to the extent permitted by law.
