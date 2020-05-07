# deconwolf

`deconwolf` is a program written in C for deconvolution of fluorescent wide-field images. It is reasonably fast and memory efficient. Optional tiling enables it to be run on machines with low memory. The most critical parts are parallelised to make use of high-core-count machines.

## Todo
 - [ ] Protect better against misuse.
 - [ ] Documentation, examples and test data.
 - [ ] Utilities for processing of multiple images and handling of PSFs.
 - [ ] Flag to save one image after each iteration (to facilitate setting the number of iterations).
 - [ ] Demos, for example on the effect of the tiling.
 - [ ] Use tif tags to write meta data (also to transfer from input image).
 - [ ] Make sure that it can be compiled on windows and mac...
 - [x] Crop PSF also in x and y (individual crop per tile as well).
 - [x] Faster calculation of tiling weights.
 - [x] Report VmPeak to the log, i.e., the peak memory usage.
 - [x] Automatic cropping of the PSF if it has too many stacks.
 - [x] Block mode for low memory systems, accessible through the options `--tilesize` and `--tilepad`
 - [x] Eliminate one or two arrays in the main loop to save memory.
 - [x] Use some kind of logging
 - [x] Argument parsing 
 - [x] Proper automatic name on output file.
 - [x] Make use of symmetries to save memory?
 - [x] Save FFTW wisdom.
 - [x] Identical results to matlab code.
Will not do:
 - Re-use plans to save some (micro) time.

## Building
deconwolf requires `libtiff` and `fftw3` to run, usually those libraries are already installed. However you might need the header files. It is good to build it on your specific machine to get the last drops of juice out of your CPU.

On Ubuntu this might be enough:
```
sudo apt-get update

# sudo apt-cache search libtiff 
sudo apt-get install libtiff5-dev

# sudo apt-cache search libfftw3
sudo apt-get install libfftw3-dev

make deconwolf -B
```

Optionally you can build fftw3 from source, the you will need to build with support for threads and single precision calculations:
```
# download source first ...
./configure --help
./configure --enable-threads --enable-single
make
sudo make install
```
## Usage:
deconwolf reads tif images (16-bit unsigned and 32-bit floats with some limitations). Usage is as simple as:
```
# for options, see
deonwolf --help
# basic usage
deconwolf dapi_001.tif PSF_dapi.tif
```
however the tricky part is to find a good points spread function (PSF).

### PSF
PSFs can be generate from ImageJ with a [plugin](http://bigwww.epfl.ch/algorithms/psfgenerator/). If the image has N slices, it is recommended that the PSF has 2xN-1 slices. Unfortunately that will require a lot of memory and processing out of your machine. If the PSF is larger than required it will automatically be cropped.

## Notes
 * FFTW is self tuning and will perform some tuning every time it presented for a new problem size. The result of this tuning is called wisdom and is stored in files like `fftw_wisdom_float_threads_16.dat` by deconwolf. Do not transfer that file to other machines and expect the tuning to take some time.
 * Memory consumption is somewhere between 60 and 70 B per voxel, and drops when tiling is enabled. (Some data: 26169630 voxels, VmPeak 1700132 kB gives 65 B per voxel. 339982580 voxels and VmPeak 19498732 gives 58 B per voxel.


## Resources
 * [fftw3 documentation](http://www.fftw.org/fftw3_doc/).
 * [libtiff repository](https://gitlab.com/libtiff/libtiff)

The algoritm is based on these papers:
 * Biggs, D.S.C. “Acceleration of Iterative Image Restoration Algorithms.” Applied Optics. Vol. 36. Number 8, 1997, pp. 1766–1775. 
 * M. Bertero and P. Boccacci, A simple method for the reduction of boundary effects in the Richardson-Lucy approach to image deconvolution, 
A&A 437, 369-374 (2005), [doi](https://doi.org/10.1051/0004-6361:20052717)
 * Lee, Ji-Yeon & Lee, Nam-Yong. (2014). Cause Analysis and Removal of Boundary Artifacts in Image Deconvolution. Journal of Korea Multimedia Society. 17. 838-848. [doi](https://doi.org/10.9717/kmms.2014.17.7.838).


