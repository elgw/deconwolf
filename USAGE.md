 - [Usage](#usage)
 - [Question and Answers](#qa)


<a name="usage" />

# Usage
Deconwolf has a command line interface (CLI), i.e., you would typically run
it from a terminal. However, here is also a
[GUI](https://github.com/elgw/deconwolf-gui) that might be handy.

## Command line usage:
To deconvolve an image (say `dapi_001.tif`) you need an approximation of
the PSF of your particular microscope. If you don't
have one you can create one according to the Born-Wolf model
with `dw_bw` that is shipped with deconwolf.

The basic information that you need to know is
 * The numerical aperture, NA.
 * The refractive index of the immersion, ni. Common media are air (ni = 1)
   and oil (ni = 1.51) and silicon oil (ni = 1.405).
 * The size of the pixels in the image.
 * The distance between the images/planes in z.

For example, if you have NA = 1.45, ni = 1.515 pixel size = 130 nm,
distance between planes = 250 nm, generate a PSF (`PSF_DAPI.tif`) with:

``` shell
dw_bw --resxy 130 --resz 250 --lambda 461 --NA 1.45 --ni 1.515 PSF_DAPI.tif
```

To deconvolve the image, type:

``` shell
dw --iter 40 dapi_001.tif PSF_dapi.tif
```

that will produce a new image called `dw_dapi_001.tif` along with
a log file
called `dw_dapi_001.tif.log.txt` which in particular say if the output
image was scaled. Since you deserve better than the
default settings, see other options:

``` shell
dw_bw --help
dw --help
```
or

``` shell
man dw
man dw_bw
```

Please note that deconwolf requires that the pixel size is the same for
both the PSF and the input image and does not attempt to read any metadata
from the tif files.

## Test data
There is a single nuclei in the `demo` subfolder (descibed in README.md).
You can get more images from the
[DeconvolutionLab2](http://bigwww.epfl.ch/deconvolution/deconvolutionlab2/)
web page. Please note that in some cases the PSFs do not contain enough z-planes.

## Memory considerations
The peak memory usage is written at the end of the log file. If you have 16
GB of RAM, images up to [1024x1024x60] pixels should work without tiling.

## PSF considerations
deconwolf requires that the PSF is centered, i.e.,
that the largest value is in the middle.
If you generate the PSF with some
other program you might have to center it first.

## Supported image formats
Currently deconwolf does only support tif images,
specifically: multipage, 16-bit unsigned or 32-bit floats,
written in strip mode. The output is either 16-bit unsigned or 32-bit
floating point and can be read by Matlab, ImageJ, etc.
If you use 16-bit output note that images will be scaled in order to not be
saturated. The scaling value can be found at the end of the log files.

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

   This self-tuning can take considerable time but should only be needed
   once per problem size.

<a name="qa" />

# Questions and Answers

## How can I prepare tif files for deconwolf?
Deconwolf only accept tif files with one channel per file. For ND2 files
check out [randiantkit](https://github.com/ggirelli/radiantkit).

## How do I get metadata from ND2-files?
You could use the
[command line tools](https://www.openmicroscopy.org/bio-formats/downloads/) from
openmicroscopy. Then you could do something like:

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
