# deconwolf

`deconwolf` is a program written in C for deconvolution of fluorescent wide-field images. It is reasonably fast and memory efficient and can be run on systems with low memory due to optional tiling.

## Todo
 - [ ] Automatic cropping of the PSF if it is too large.
 - [ ] Protect better against misuse.
 - [ ] Documentation, examples and test data.
 - [ ] Eliminate one or two arrays in the main loop to save memory.
 - [ ] Utilities for processing of multiple images and handling of PSFs.
 - [ ] Flag to save one image after each iteration (to facilitate setting the number of iterations).
 - [ ] Demos, for example on the effect of the tiling.
 - [x] Block mode for low memory systems, accessible through the options `--tilesize` and `--tilepad`
 - [x] Use some kind of logging
 - [x] Argument parsing 
 - [x] Proper automatic name on output file.
 - [x] Make use of symmetries to save memory?
 - [x] Save FFTW wisdom.
 - [x] Identical results to matlab code.
Low priority:
 - [ ] Re-use plans to save some (micro) time.
 - [ ] Faster calculation of tiling weights.

## Building
deconwolf requires `libtiff` and `fftw3` to run, usually those libraries are already installed. However you might need the header files.

On Ubuntu:
```
sudo apt-get update

# sudo apt-cache search libtiff 
sudo apt-get install libtiff5-dev

# sudo apt-cache search libfftw3
sudo apt-get install libfftw3-dev

make deconwolf -B
```

If you decide to build fftw3 from source, these commands worked for me:
```
# download source first ...
./configure --help
./configure --enable-threads --enable-single
make
sudo make install
```
## Usage:
```
# for options, see
deonwolf --help
# basic usage
deconwolf dapi_001.tif PSF_dapi.tif
```

### PSF
PSFs can be generate from ImageJ with a [plugin](http://bigwww.epfl.ch/algorithms/psfgenerator/). It the image has N slices, it is recommended that the PSF has 2xN+1 slices. Unfortunately that will require a lot of memory and processing to handle.

## Notes
 * FFTW is self tuning and will perform some tuning every time it presented for a new problem size. The result of this tuning is called wisdom and is stored in `fftw_wisdom_float_threads.dat` by deconwolf. Do not transfer that file to other machines.
 * Memory consumption: 26169630 voxels and VmPeak 1700132 kB gives 65 B per voxel. 339982580 and VmPeak 19498732 gives 58 B per voxel.

## Resources
 * [fftw3 documentation](http://www.fftw.org/fftw3_doc/).
 * [libtiff repository](https://gitlab.com/libtiff/libtiff)


