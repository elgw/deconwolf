# Usage

## Don't expect too much

When working with real widefield images, in contrast to synthetic
ones, it is important to know that the success or failure of
deconvolution depends on several factors and that there are
limitations that can not be overcome.

Most important, the results will always be limited by mismatches
between the PSF of the microscopy and what is given as input to the
software (PSF mismatch).

The "Born and Wolf" PSF model that is included with deconwolf should
be a good, or at least simple, starting point since it only requires
very few parameters. Please be aware there is a compromise between
simplicity and accuracy, especially the Born and Wolf model relies on
several assumptions that typically are met. For instance you can see
that such PSFs are symmetric both around the lateral plane and around
the z-axis. That is something I've never seen in real images.

It is likely that better results can be obtained by using more
advanced models (like the Gibson-Lanni) or by experimentally
estimating the PSF based on image data (of beads). However, to find a
better PSF for a specific microscope is typically quite
complicated.

Richardson-Lucy based methods assume that the image formation can be
described as a linear process, something which is only approximate
since any interesting sample will both absorb and refract
light. Saturated pixel in the input image clearly violates the
linearity assumption so please check for that before deconvolution.

Saturated pixels, sensor noise, sensor artifacts, PSF mismatches, an
excessive number of iterations, ..., can cause artifacts. Ideally one
should acquire samples also with super resolution, to be in a position
to recognize and distinguish interesting features from potential
artefacts.

> [!CAUTION]
> Deconvolution is not guaranteed to improve microscopy images
> and might introduce artefacts.


## Running deconwolf
Deconwolf has a command line interface (CLI), and can be run from a
terminal or script. There is also a limited
[GUI](https://github.com/elgw/deconwolf-gui) that might be handy, at
least initially.


To deconvolve an image (let's call it `dapi_001.tif`) you need an
approximation of the PSF of your particular microscope at the given
emission wave length.

If you don't have one you can create one according to the Born-Wolf
model with the `dw_bw` program that is shipped with deconwolf.

The basic information that you need to know is:
- The numerical aperture, NA.
- The refractive index of the immersion, ni. Common media are air (ni = 1)
and oil (ni = 1.51) and silicon oil (ni = 1.405).
- The lateral size of the pixels in the __image__ -- NOT the size of
the pixels in the sensor.
- The axial size of the pixels in the image, i.e. distance between the
images/planes.

For example: one of our 100x objectives has NA = 1.45, ni = 1.515. At
the given magnification that gives us image pixel of size 130 nm. The
distance between planes was set to 250 nm. To generate a PSF with this
information (we will call it `PSF_dapi.tif`) the following command can
be used:

``` shell
dw_bw --resxy 130 --resz 250 --lambda 461 --NA 1.45 --ni 1.515 PSF_dapi.tif
```

> [!TIP]
> You can of course re-use the PSF with multiple images. Please make sure
> to keep the `.log.txt` file, generated along with the PSF, to know what
> parameters were used.

> [!NOTE]
> Please note that deconwolf does not read any meta data from the tif
> files regarding pixel size etc. Hence it is the PSF file that determines
> those factors implicitly.

With a PSF at hand the image can be deconvolved with a command line like this:

``` shell
dw --iter 40 dapi_001.tif psf_dapi.tif
```

When done, **dw** will write a new image called `dw_dapi_001.tif`
along with a log file called `dw_dapi_001.tif.log.txt`.

> [!CAUTION]
> By default, the pixel values of the output image will be scaled to
> make most out of the 16-bit image format. If absolute intensities matter,
> please read about the various output options below.

For a list of the available command line options, please see:

``` shell
dw_bw --help
dw --help
```

On Linux/Mac there should also be a manual installed, which can be
opened by:

``` shell
man dw
man dw_bw
```

If the manuals are not available on your system, they can also be
found at this link
[https://github.com/elgw/deconwolf/tree/master/doc].


## Test data

If you want to test that **dw** works, there is a small image of a single nuclei
in the `demo` subfolder. The image is small enough that it can be
deconvolved with very modest hardware.

More test images can be found on the
[DeconvolutionLab2](http://bigwww.epfl.ch/deconvolution/deconvolutionlab2/)
web page. Please note that in some cases the PSFs do not contain
enough z-planes for a proper deconvolution.

## Memory considerations
If you have 16 GB of RAM, images up to [1024x1024x60] pixels should
work without tiling. If deconwolf crashes because it runs of out
memory you can:

1. Consider removing top and bottom slices from the image. The
boundary handling in the axial direction is quite good so often
there is no need for a lot of out of focus images planes above and
below the objects of interest.
2. Use the **--tilesize N** option. To process the image as slightly
overlapping tiles.
3. Use fewer threads, i.e. specify **--threads N**.
4. Consider switching to a lower quality on the boundaries with
**--bq 1**. The option **--bq 0** is not recommended for normal
images (since it assumes that the image repeats itself around the
edges which is bad particularly in the axial direction).

The peak memory usage is written at the end of the log file (at least
if run on Linux).

## PSF considerations
The **dw_bw** program will generate PSFs that can be used with
deconwolf, but it is also possible to use PSF from other sources.

deconwolf requires that the PSF is centered, i.e., that the largest
value is in the middle.  If you generate the PSF with some other
program you might have to center it first. In our experience it is
better to have a PSF with odd sizes, i.e. that the PSF is centered at
a certain pixel and not between them. Ideally the PSF should have at
least $2P-1$ planes where $P$ is the number of planes of the input image.

If the PSF has more planes than needed, it will be cropped
automatically when loaded. If it has fewer than the ideal number of
planes, an warning will be issued.

The PSF might also be cropped in the lateral plane when loaded (to
save some memory). To prevent that, add the command line argument
**--xyfactor 0**.

## Supported image formats
Currently deconwolf can only read tif images, specifically:
multi-page, 8-bit unsigned, 16-bit unsigned or 32-bit floats. The data
should be written in strip mode (typically the default output
format). The output is either 16-bit unsigned or 32-bit floating point
(with the **--float** flag) and can be read by Matlab, ImageJ, etc.
If you use 16-bit output note that images will be scaled in order to
not be saturated. The scaling value can be found at the end of the log
files.

It is undefined behavior to use pyramidal (multiple resolution) or
multi-color images.

## Log files
- The reported error, Idiv, is the I-divergence between the input image
  and the current guess convolved with the PSF. Please note that deconwolf
  doesn't know how the deconvolved images should look like so this isn't
  in any way a measurement of how good image quality you get out of the
  program.

- If deconwolf finished normally you will find that it ends with something
  like:
  ``` shell
  Iteration  24/ 25, Idiv=5.539379e-01
  Iteration  25/ 25, Idiv=5.426050e-01
  Deconvolution took 15.52 s
  2.452803% pixels at bg level in the output image.
  scaling: 2.975924
  Took: 17.723509 s
  peakMemory: 18486288 kiB
  2024-03-08T11:33:34
  ```
  If it doesn't, chances are that deconwolf run out of memory and crashed.

## Output format

Internally deconwolf use 32-bit floating point precision but when it
is time to write the output image to disk 16-bit unsigned tif images
are used by default since that is more memory efficient, and 16-bits
per pixel should be good enough in most cases.

To prevent clipping/saturation and to make full use of the 16-bit
dynamic range all pixel values are scaled by a scaling factor, $s$,
prior to writing. The value of the scaling factor is based on the
deconvolved image, $Y$, and set as $s=(2^{16}-1)/\max(Y)$.

The scaling is reported in the
log file like:
``` shell
$ cat dw_dapi_001.tif.log.txt | grep scaling
scaling: 0.425100
```

If you care about the absolute intensities there are three things you
can do:

1. Use the **--float** argument to save images as 32-bit floating
point tif files. The downside is that the output images will be twice
as large compared to the default and that the support for reading
such images is not widespread.

2. Set a scaling factor manually by using **--scale s** where $v$ is
some value that you have to figure out yourself. $s==1$ is a
natural choice to start with but might cause saturated pixels. If
that is the case reduce $s$ until you are satisfied. For example,
if you use **--scale 0.5** the output pixels will have only half of
the computed value, if you load the images you can then multiply
the pixel values by $2$ to get the original values back.

3. Read the scaling value from the log files and divide by that value
to get the computed intensities back.

## Questions and Answers

### I have a GPU that I want to use, what can I do?
Congratulations, that will make deconwolf run much faster. In general
you only need to supply the command line argument **--gpu** to enable
GPU acceleration. However, there are some things to take into consideration:

- Was deconwolf compiled with GPU support? If you run `dw --version`
the output should include:

  ``` shell
  OpenCL: YES
  VkFFT: YES
  ```
  if not, you need to re-build with the GPU options turned on.

- Does OpenCL work on you machine? On linux you can check with the
command line program `clinfo`. In some cases you need to download
specific drivers/software from the GPU vendor.

- Do you have enough memory? The GPU based pipeline use about as much
GPU memory as the non-gpu pipeline use RAM. If you first deconvolve
an image with without the **--gpu** option you can check the log
file (at least under linux) to figure out how much RAM was used. If
there isn't enough memory on the GPU there is the option to
deconvolve using tiles, see the info about the **--tilesize** option.

- The results will not be pixel identical to images deconvolved
without the **--gpu** option. However differences should be
negligable (if not, please report a bug). This is because the FFT
libraries have a limited precision.

- Do you have more than one GPU or "cldevice"? Then you can tell
  **dw** which one to use with the `--cldevice n` argument, where $n$
  is the number of the device to use, the first one has number $0$.

### Do I really need to specify the number of iterations?

RL-based deconvolution methods are iterative and don't know when to
stop. While some programs have automatic stop conditions we argue that
it is better if the user specify the number of iterations for two reasons:

1. How many iterations that are wanted/needed/optimal depends on
downstream analysis. And who knows what your purposes are?
2. For consistency. Given a dataset with multiple FOV is makes sense
that each FOV is deconvolved with the same number of
iterations. With an automatic stop functionality, that is probably
not the case.

### It takes a long time before deconwolf starts iterating

FFTW is self tuning and perform some optimizations whenever
presented to a new problem size. The result of this tuning is called
wisdom and is stored in files like `fftw_wisdom_float_threads_16.dat`
by deconwolf in `~/config/deconwolf/`. Those files can in general not
be transferred to other machines.

This self-tuning can take considerable time but should only be needed
once per problem size. To turn of the self-tuning use the **--noplan**
flag.


### How can I prepare tif files for deconwolf?

Deconwolf only accept tif files with one channel per file. ImageJ/FIJI
is typically helpful for conversions. For ND2 files check out
[randiantkit](https://github.com/ggirelli/radiantkit) or
[nd2tool](https://github.com/elgw/nd2tool).

### How do I get metadata from ND2-files?

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

### How can I run deconwolf for many files?

Either use the [GUI](https://github.com/elgw/deconwolf-gui), or write your own
script in your favorite language.

### Don't you have more PSF models?

No, not at this time. Please check out the
[PSF Generator](http://bigwww.epfl.ch/algorithms/psfgenerator/) from EPFL.


### What about image dimensions?

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


### Maximizing throughput

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

### 2D images?
No, the option to deconvolve 2D images is deliberately turned off.

### My question was not in this list!

Then you are welcome to ask it at [https://github.com/elgw/deconwolf/issues].
