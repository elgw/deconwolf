 - [Usage](#usage)
 - [Question and Answers](#qa)


<a name="usage" />

# Usage
Deconwolf has a command line interface (CLI), i.e., you would typically run
it from a terminal. However, here is also a
[GUI](https://github.com/elgw/deconwolf-gui) that might be handy.

## Command line usage:
To deconvolve an image (say `dapi_001.tif`) you need an approximation
of the PSF of your particular microscope. If you don't have one you
can create one according to the Born-Wolf model with `dw_bw` that is
shipped with deconwolf.

The basic information that you need to know is
 * The numerical aperture, NA.
 * The refractive index of the immersion, ni. Common media are air (ni = 1)
   and oil (ni = 1.51) and silicon oil (ni = 1.405).
 * The size of the pixels in the image.
 * The distance between the images/planes in z.

For example, if you have NA = 1.45, ni = 1.515 pixel size = 130 nm,
distance between planes = 250 nm, generate a PSF (`PSF_DAPI.tif`) with:

``` shell
dw_bw --resxy 130 --resz 250 --lambda 461 --NA 1.45 --ni 1.515 psf_dapi.tif
```

To deconvolve the image, type:

``` shell
dw --iter 40 dapi_001.tif psf_dapi.tif
```

that will produce a new image called `dw_dapi_001.tif` along with a
log file called `dw_dapi_001.tif.log.txt` which in particular say if
the output image was scaled. We have tried to make the default options
as sane as possible, but for more options see:

``` shell
dw_bw --help
dw --help
```
for a brief overview or

``` shell
man dw
man dw_bw
```
for more details.

Please note that deconwolf requires that the pixel size (nm per pixel)
is the same for both the PSF and the input image and does not attempt
to read any metadata from the tif files.

## Test data
There is a single nuclei in the `demo` subfolder (descibed in
README.md).  You can get more images from the
[DeconvolutionLab2](http://bigwww.epfl.ch/deconvolution/deconvolutionlab2/)
web page. Please note that in some cases the PSFs do not contain
enough z-planes for a proper deconvolution.

## Memory considerations
If you have 16 GB of RAM, images up to [1024x1024x60] pixels should
work without tiling. If deconwolf crashes because it runs of out
memory you can:

 1. Use the *--inplace** flag (FFTs will be just slightly slower while
    memory will be reduced).
 2. Use fewer threads, i.e. specify **--threads N**.
 3. Consider switching to a lower quality on the boundaries with
    **--bq 1**. The option **--bq 0** is not recommended for normal
    images (since it assumes that the image repeats itself around the
    edges which is bad particularly in the axial direction).

Accordingly the lowest memory usage will be when you combine these
things i.e. calling dw with **--inplace --threads 1 --bq 1**. If that
still uses too much memory you have to tell dw to process the image in
tiles with the **--tilesize N** option.

The peak memory usage is written at the end of the log file.

## PSF considerations
deconwolf requires that the PSF is centered, i.e., that the largest
value is in the middle.  If you generate the PSF with some other
program you might have to center it first. In our experience it is
better to have a PSF with odd sizes, i.e. that the PSF is centered at
a certain pixel and not between them.

## Supported image formats
Currently deconwolf does only support tif images, specifically:
multi-page, 8-bit unsigned, 16-bit unsigned or 32-bit floats. The data
should be written in strip mode (typically the default output
format). The output is either 16-bit unsigned or 32-bit floating point
(with the **--float** flag) and can be read by Matlab, ImageJ, etc.
If you use 16-bit output note that images will be scaled in order to
not be saturated. The scaling value can be found at the end of the log
files.

## Log files and output
 * The reported error is the mean square error between the input image
   and the current guess convolved with the PSF. Please note that deconwolf
   doesn't know how the deconvolved images should look like so this isn't
   in any way a measurement of how good image quality you get out of the
   program.
 * If the 16-bit output format is used, the scaling is reported in the
   log file like:
   ```
   $ cat dw_dapi_001.tif.log.txt | grep scaling
   scaling: 0.425100
   ```
   That means, in order to go back to absolute intensities, divide by
   that number.
 * If deconwolf finished normally you will find that it ends with something
   like:
   ```
   Iteration 19/20, MSE=2.481362e+04
   Iteration 20/20, MSE=2.235636e+04
   1.957566% pixels at bg level in the output image.
   scaling: 0.425100
   Took: 28.767448 s
   peakMemory: 4420728 kiB
   Fri Jan 22 16:42:22 2021
   ```
   If it doesn't, chances are that deconwolf run out of memory and crashed.

## Notes
 * FFTW is self tuning and will perform some tuning every time it
   presented for a new problem size. The result of this tuning is called
   wisdom and is stored in files like `fftw_wisdom_float_threads_16.dat`
   by deconwolf in `~/config/deconwolf/`. Those files can in general
   not be transferred to other machines.

   This self-tuning can take considerable time but should only be
   needed once per problem size. To turn of the self-tuning use the
   **--noplan** flag.

<a name="qa" />

# Questions and Answers

## How can I prepare tif files for deconwolf?
Deconwolf only accept tif files with one channel per file. ImageJ/FIJI
is typically helpful for conversions. For ND2 files check out
[randiantkit](https://github.com/ggirelli/radiantkit) or
[nd2tool](https://github.com/elgw/nd2tool).

## How do I get metadata from ND2-files?
You could use either [nd2tool](https://github.com/elgw/nd2tool) or the
more mature [showinf](https://www.openmicroscopy.org/bio-formats/downloads/) from
openmicroscopy. With showinf you could do something like:

```
$ bftools/showinf -nopix -omexml iAM337_20190830_001.nd2 > omexml

$ cat omexml | grep -B 1 -A 2 dObjectiveNA
            <OriginalMetadata>
               <Key>dObjectiveNA</Key>
               <Value>1.45</Value>
            </OriginalMetadata>

$ cat omexml | grep EmissionWavelength | grep Channel:0:
         <Channel Color="-16711681" EmissionWavelength="810.0" EmissionWavelengthUnit="nm" ID="Channel:0:0" Name="ir800" SamplesPerPixel="1">
         <Channel Color="-2147467009" EmissionWavelength="700.0" EmissionWavelengthUnit="nm" ID="Channel:0:1" Name="a700" SamplesPerPixel="1">
         <Channel Color="16056319" EmissionWavelength="488.0" EmissionWavelengthUnit="nm" ID="Channel:0:2" Name="a488" SamplesPerPixel="1">
         <Channel Color="-16776961" EmissionWavelength="695.0" EmissionWavelengthUnit="nm" ID="Channel:0:3" Name="Cy5" SamplesPerPixel="1">
         <Channel Color="-1241579265" EmissionWavelength="542.0" EmissionWavelengthUnit="nm" ID="Channel:0:4" Name="tmr" SamplesPerPixel="1">
         <Channel Color="-9109249" EmissionWavelength="590.0" EmissionWavelengthUnit="nm" ID="Channel:0:5" Name="a594" SamplesPerPixel="1">
         <Channel Color="570490879" EmissionWavelength="432.0" EmissionWavelengthUnit="nm" ID="Channel:0:6" Name="dapi" SamplesPerPixel="1">

$ cat omexml | grep PhysicalSizeZ= | head -n 1
      <Pixels BigEndian="false" DimensionOrder="XYCZT" ID="Pixels:0" Interleaved="false" PhysicalSizeX="0.129780110998775" PhysicalSizeXUnit="µm" PhysicalSizeY="0.129780110998775" PhysicalSizeYUnit="µm" PhysicalSizeZ="0.2" PhysicalSizeZUnit="µm" SignificantBits="16" SizeC="7" SizeT="1" SizeX="1024" SizeY="1024" SizeZ="81" Type="uint16">
...
```

## How can I run deconwolf for many files?
Either use the [GUI](https://github.com/elgw/deconwolf-gui), or write your own
script in your favorite language.

## Don't you have more PSF models?
No, not at this time. Please check out the
[PSF Generator](http://bigwww.epfl.ch/algorithms/psfgenerator/) from EPFL.


## What about image dimensions?

**dw** report image dimensions in the same way as ImageJ (although dw does
not understand time and channel dimensions). For example the file
`dapi_001.tif` in the `demo/` folder is reported to have the shape
**101x201x40** by dw and ImageJ.

Other tools might report it the other way around, for example:

skimage reports the size **40x201x101** by
``` python
>>> from skimage import io
>>> I = io.imread('dapi_001.tif')
>>> I.shape
```

In general this is nothing to worry about, see
[https://en.wikipedia.org/wiki/Row-_and_column-major_order]
for a discussion around this topic.

## Maximizing throughput
To get a high throughput, i.e., many processed images / hour it is
typically a little faster to process many images in parallel (if
enough RAM).  Here are some example timings based on images:
[2048x2048x20], psf: [117x117x39] **\--iter 20**, **\--bq 2**. Using
dw version 0.1.0 on an AMD Ryzen 3700X:

 - One dw (using **\--threads 8**) took 88 s using 7294 MB RAM.

 -> 40 images / hour or 88 s / image

 - Two dw in parallel (using **\--threads 4**) took 2 m 44 s using 2x7106 =
   14212 MB RAM

 -> 43 images / hour or 82 s / image

 - Four dw in parallel (using **\--threads 2**) took 297 s using 6795x4 =
   27181 MB RAM

 -> 48 images / hour or 74 s / image

 - Eight dw in parallel (using **\--threads 1**) took 9 m 12 s, using
   6712x8 = 53696 MB RAM

-> 52 images / hour or 69 s / image


It is also possible to mix a few instances of dw using the CPU and a
few offloading calculations to the GPU (with **--method
shbcl**). Results would depends on the GPU/CPU combination.
