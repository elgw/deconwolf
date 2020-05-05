# deconwolf

`deconwolf` is a program written in C for deconvolution of fluorescent wide-field images.

## Todo
 - [ ] Argument parsing 
 - [ ] Proper automatic name on output file.
 - [ ] Automatic cropping of the PSF if it is too large.
 - [ ] Protect better against misuse.
 - [ ] Re-use plans to save some (micro) time.
 - [ ] Use logger
 - [ ] Documentation, examples and test data.
 - [ ] Eliminate one or two arrays in the main loop to save memory.
 - [ ] block mode for low memory systems
 - [x] Make use of symmetries to save memory?
 - [x] Save FFTW wisdom.
 - [x] Identical results to matlab code.

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
bin/deconwolf dapi_001.tif PSF_dapi.tif
```

## Notes
 * FFTW is self tuning and will perform some tuning every time it presented for a new problem size. The result of this tuning is called wisdom and is stored in `fftw_wisdom_float_threads.dat` by deconwolf. Do not transfer that file to other machines.
 * Memory consumption: 26169630 voxels and VmPeak 1700132 kB -> 65 B per voxel.

## Resources
 * [fftw3 documentation](http://www.fftw.org/fftw3_doc/).
 * [libtiff repository](https://gitlab.com/libtiff/libtiff)


