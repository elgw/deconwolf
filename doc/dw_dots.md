% DECONWOLF(1) dw ABCDE
% Erik Wernersson
% 2022

# NAME
deconwolf **dw** is a tool for deconvolution of wide field microscopy
stacks.

**dw dots** is the module that can be used to detect dots in tif
images and export them to a tsv file.

# SYNOPSIS

**dw dots** [*OPTIONS*] file1.tif file2.tif ...

# OPTIONS
**\--LoG**
: Pre-filter the image with a Laplacian of Gaussian
  filter. Recommended for non-deconvolved images. Sometimes good for
  deconvolved as well.

**\--lsigma**
: Lateral sigma for a Gaussian low pass filter or the LoG filter if enabled.

**\--asigma**
: Axial sigma for a Gaussian low pass filter or the LoG filter if enabled.

**\--overwrite**
: If specified existing files will be overwritten.

**\--threads t**
: Set the number of computational threads to use.

**\--fout file.tif**
: Write out the filtered images (e.g. the LoG filtered image). Only
  for debugging and only to be used when there is just one input file.

**\--ndots**
: Specify the max amount of dots that should be exported.

**\--help**
: Show a brief help message

**\--nthreads**
: Specify the number of threads to use.

**\--verbose**
: Set the verbosity level.

# Method
**dw dots** will read one input image at a time and will perform the following:

 - Filter the image with either a Gaussian low pass filter or a LoG
  filter using the values specified by **\--lsigma** and **\--asigma**.
 - Extract all local maxima.
 - Sort the maxima according to the dot intensity.
 - Calculate the FWHM value for each dot in the lateral plane (by
   intersection).
 - Write out at most **\--ndots** dots to a TSV file.

# INPUT
**dw dots** can read the same kind of files as **dw**. It does look
for log files from deconwolf, and if one is found the input image is
rescaled according to the log file.

For deconvolved images **\--lsigma** and **\--asigma** can typically
be left untouched, i.e. at 0. For non-deconvolved images it makes
sense to match the filters to the airy pattern.

# OUTPUT
A TSV file will be generated per input file (`.dots.tsv`) as well as a
log file (`.dots.log.txt`). The TSV file will be sorted by the dot intensities.

# SEE ALSO
**dw**, **dw_bw**

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
