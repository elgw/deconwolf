# deconwolf

`deconwolf` is a program written in C for deconvolution of fluorescent wide-field images. It is reasonably fast and memory efficient. Optional tiling enables it to be run on machines with low memory. The most critical parts are parallelised to make use of high-core-count machines.


## Building
deconwolf requires `fftw3f`, `fftw3f_threads` and `tiff-5` and can be built with meason. 

The installation is typically something like this:
```
meson builddir
ninja -C builddir
# to install deconwolf to a standard location, use
sudo ninja install
# if you for some reason don't want it anymore
sudo ninja uninstall
```

If meson or some of the libraries are missing, use google search! On Ubuntu 19.10 this did the job:
```
sudo apt-get update

# sudo apt-cache search libtiff 
sudo apt-get install libtiff5-dev

# sudo apt-cache search libfftw3
sudo apt-get install libfftw3-dev

sudo apt-get install meson
```

If you need to build fftw3 from source, that was not to tricky:
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
There is a small script that might be useful for batch processing, example:
```
$ deconwolf_batch.py iMERULA91_20200125_001/ ../PSF/ '--iter 50 --tilesize 512'
deconwolf --iter 50 --tilesize 512 iMERULA91_20200125_001/Cy5_012.tif ../PSF/PSF_Cy5.tif
deconwolf --iter 50 --tilesize 512 iMERULA91_20200125_001/Cy5_013.tif ../PSF/PSF_Cy5.tif
deconwolf --iter 50 --tilesize 512 iMERULA91_20200125_001/dapi_013.tif ../PSF/PSF_dapi.tif
...
```
the output can be run directly with `| bash` or piped to a file and executed later (possibly using `parallel`). Since deconwolf does not overwrite output files such list of jobs to be run can be aborted and restarted.

### Memory considerations
Memory consumption is somewhere between 60 and 70 B per voxel, and drops when tiling is enabled. (Some data: 26169630 voxels, VmPeak 1700132 kB gives 65 B per voxel. 339982580 voxels and VmPeak 19498732 gives 58 B per voxel.

### PSF
PSFs can be generate from ImageJ with a [plugin](http://bigwww.epfl.ch/algorithms/psfgenerator/). If the image has N slices, it is recommended that the PSF has 2xN-1 slices. Unfortunately that will require a lot of memory and processing out of your machine. If the PSF is larger than required it will automatically be cropped.

## Notes
 * FFTW is self tuning and will perform some tuning every time it presented for a new problem size. The result of this tuning is called wisdom and is stored in files like `fftw_wisdom_float_threads_16.dat` by deconwolf (in `~/config/deonwolf/`). Do not transfer that file to other machines and expect the tuning to take some time.


## Resources
 * [fftw3 documentation](http://www.fftw.org/fftw3_doc/).
 * [libtiff repository](https://gitlab.com/libtiff/libtiff)

The algoritm is based on these papers:
 * Biggs, D.S.C. “Acceleration of Iterative Image Restoration Algorithms.” Applied Optics. Vol. 36. Number 8, 1997, pp. 1766–1775. 
 * M. Bertero and P. Boccacci, A simple method for the reduction of boundary effects in the Richardson-Lucy approach to image deconvolution, 
A&A 437, 369-374 (2005), [doi](https://doi.org/10.1051/0004-6361:20052717)
 * Lee, Ji-Yeon & Lee, Nam-Yong. (2014). Cause Analysis and Removal of Boundary Artifacts in Image Deconvolution. Journal of Korea Multimedia Society. 17. 838-848. [doi](https://doi.org/10.9717/kmms.2014.17.7.838).

## Todo
 - [ ] Protect better against misuse.
 - [ ] Documentation, examples and test data.
 - [ ] A higher level interface with facilities for handling and generation of PSFs.
 - [ ] Flag to save one image after each iteration (to facilitate setting the number of iterations).
 - [ ] Demos, for example on the effect of the tiling.
 - [ ] Use tif tags to write meta data (also to transfer from input image).
 - [ ] Make sure that it can be compiled on windows and mac...
 - [ ] Use cmake to facilitate cross platform builds
 - [ ] Crash-safe writing, write to temporary file and move when write is complete to avoid bad luck.
 - [ ] Highly parsable log file.
 - [ ] Break out the image processing functions to separate library.
 - [ ] Proper logger.
 - [x] Possible to change the prefix with the `--prefix` flag.
 - [x] Refuse to run if output file already exist (unless `--overwrite`) with status `0`.
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
 - [x] Save FFTW wisdom (saved to `./config/deconwolf/`)
 - [x] Identical results to matlab code.

Will not do:
 - Re-use plans to save some (micro) time.

### Utils (that does not exist yet)

```
# generate a list of things to run
deconwolf --batch (--iter ...) image_folder psf_folder
# run them
./deconwolf_batch
# possibly like this (if you know what you are doing)
parallel --halt-on-error 2 --jobs $nCores --joblog parallel.log < deconwolf_batch
```

