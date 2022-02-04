# dw_bw

```
INPUT:
   Parameters to the Born-Wolf model
   Parameters about he integration
   Parameters about the ouput PSF shape

OUTPUT:
   A 32-bit tiff file with the PSF
   A .log.txt file containing everything to
   reproduce the result

For each image plane, p, in the PSF: # Processed in parallel
    Sample the BW model for discrete radii -> RZ_p
    For each pixel in the output PSF at plane p:
       Integrate RZ_p (using lanczos-5 interpolation) over the pixel

Save the output image and log file.
```

# dw

```
INPUT:
   The image to deconvolve of size [M x N x P]
   The PSF to use of size [m x n x p]
   Deconvolution settings.
OUTPUT:
   A deconvolved image (prefixed by dw_ by default)
   A .log.txt file containing everything to
   reproduce the result.

# Don't use a deeper PSF than necessary
if p > 2 * P - 1:
   crop the PSF so that p = 2 * P - 1

# Don't use a wider PSF than necessary
while the peripheral XZ and YZ planes have a sum > xyfactor:
   remove peripheral XZ and YX planes

if M > tilesize or N > tilesize:
   process the image in tiles.

deconvolve image/tiles

save image to disk as either 16-bit or 32-bit.
```

# deconvolution

```
INPUT:
   Image data to be deconvolved [M' x N' x P']
   PSF to use [m x n x p]

# Determine the work size, i.e. size of padded image.
If bq == 0:
   M = M'
   N = N'
   P = P'

if bq == 1:
   M = M' + (m + 1)/2
   N = N' + (n + 1)/2
   P = P' + (p + 1)/2

if bq == 2:
   M = M' + m - 1
   N = N' + n - 1
   P = P' + p - 1

if bq > 0:
   Set up the weight matrix W according to Bertero

Set up the auxiliary arrays required for Biggs EVE acceleration.

Until max iterations reached:
   Run one RL iteration
   Update arrays according to Biggs

Crop the final image from [M x N x P] to [M' x N' x P'] and return.
```
