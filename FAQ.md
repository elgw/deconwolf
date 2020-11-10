### How do I get metadata from ND2-files?
One option would be to download the [command line tools](https://www.openmicroscopy.org/bio-formats/downloads/) from openmicroscopy. Then you could do something like:

```
$ bftools/showinf -nopix -omexml iAM337_20190830_001.nd2 > omexml
$ cat omexml | grep sObjective
$ cat omexml | grep wsObjective

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

### How can I run deconwolf for many files at the same time
Try something like
```
find dapi_* -print0 | xargs -p -0 -I{} deconwolf --threads 32 --iter 50 {} ../PSF/PSF_dapi.tif
```
remove the -p when you know what you are doing.

Also check out the python script `deconwolf_batch.py`.

### Which theoretical PSF model is best?
Don't know.

### How to generate a bunch of PSFs ?
Have a program to loop through your settings and create `config.txt` for the PSFGenerator
```
java -cp PSFGenerator.jar PSFGenerator config.txt
```
Would-like-to-have interface:
 1. Go to folder with tiff files
 2. `genPSFs regexp` to identify channel names
 3. Asks for NA, pixel sizes and what model to use. Then ask for emission wavelengths for each channel.
 4. A folder PSF is generated.

Change `PSF-shortname=RW` to edit type, change corresponding `psf-RW-accuracy=Good` to whatever is wanted. 
```
# Number of pixels in output image
NZ=121
NY=121
NX=121
# Emission wavelength
Lambda=470 
# Lateral resolution in nm
ResLateral=130.0
ResAxial=200
```
The output name is `PSF XY.tif` where `XY` is the short name.
