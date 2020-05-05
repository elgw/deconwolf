# deconwolf

`deconwolf` is a program written in C for deconvolution of fluorescent wide-field images.

## Todo
 - [ ] Protect better against misuse.
 - [ ] Re-use plans to save some (micro) time.
 - [ ] Make use of symmetries to save memory?
 - [ ] Use logger
 - [ ] Argparsing
 - [ ] Proper automatic name on output file.

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

## Resources
 * [fftw3 documentation](http://www.fftw.org/fftw3_doc/).
 * [libtiff repository](https://gitlab.com/libtiff/libtiff)


