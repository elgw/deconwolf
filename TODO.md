 - [Prioritized](#prio)
 - [Nice to have](#nice)
 - [Maybes](#maybe)
 - [Done](#done)

# Todo

<a name="prio" />


## Top priority
- [ ] Crash-safe writing of output images, write to temporary file and
       move when write is complete to avoid bad luck.
- [ ] Get documentation up to date.
- [ ] Set default PSF size in `dw_bw` based on the geometry of the PSF
       and expected sample thickness.

<a name="nice" />


## Nice to have

- [ ] Adjust the PSF cropping so that the working image size has at least one power of two. Or even better use `nextfastfft` from [LFAT](https://ltfat.github.io/notes/ltfatnote017.pdf), [code](https://github.com/ltfat/ltfat/blob/master/fourier/nextfastfft.m)
- [ ] Wrap malloc to catch allocation errors to the log file.
- [ ] Provide test data and demos
- [ ] Decide about default tiling settings based on the PSF.
- [ ] Demos, for example on the effect of the tiling.
- [ ] Implement the Gibson-Lanni PSF model.
- [ ] Make `FFTW_ESTIMATE` an option.
- [ ] Double check the tiling weights using gradients and other patterns.
- [ ] Add information on the input image to the log -- especially detect saturated pixels.
- [ ] Make sure that something that makes sense happens when a Pyramidal tif is supplied.

<a name="maybe" />


## Maybes
- [ ] Figure out if there are any performance benefits by buildin fftw3 from source for `-march=native -mtune=native`.
- [ ] Try the FFTW interface to [cuFFT](https://docs.nvidia.com/cuda/cufft/index.html#fftw-supported-interface)

<a name="done" />


## Done

- [x] Include the "normal" `Richardson-Lucy` deconvolution method. Kind of in place with the `--biggs 0` options. However some arrays that are not needed without acceleration are still allocated.
- [x] Intel Math Kernel Library (MKL).
- [x] Write down the image scaling to the log files.
- [x] A higher level interface with facilities for handling and generation of PSFs. -- runPSFGenerator.py
- [x] Reuse more tile sizes in order to keep the number of training sessions low.
- [x] faster BigTIFF writes, examples exist in [matlab](https://github.com/rharkes/Fast_Tiff_Write/blob/master/Fast_BigTiff_Write.m), [In C++](https://github.com/jkriege2/TinyTIFF), ... -- decided that it was not necessary
- [x] Write as BigTIFF if the images don't fit into regular TIFF files.
- [x] Does gomp use more threads than it should? -- It looks like gomp creates more threads than needed, however only the specified number of threads are active at the same time.
- [x] Never load the full image in tile mode.
- [x] Identity included, using `--tilesize 0`
- [x] Reduce memory by using the identity `f(-x) = ifft(conj(fft(x))`
- [x] Break out the image processing functions to separate library.
- [x] Protect better against misuse.
- [x] Support even-sized PSFs (kind of hack, solution by resizing)
- [x] Check that it works on osx
- [x] Cross-platform build-tool (meson).
- [x] Possible to change the prefix with the `--prefix` flag.
- [x] Refuse to run if output file already exist (unless `--overwrite`) with status `0`.
- [x] Crop PSF also in x and y (individual crop per tile as well).
- [x] Faster calculation of tiling weights.
- [x] Report peak memory usage to log file.
- [x] Automatic cropping of the PSF if it has too many stacks.
- [x] Block mode for low memory systems, accessible through the options `--tilesize` and `--tilepad`
- [x] Eliminate one or two arrays in the main loop to save memory.
- [x] Use some kind of logging
- [x] Argument parsing
- [x] Proper automatic name on output file.
- [x] Make use of symmetries to save memory?
- [x] Save FFTW wisdom (saved to `./config/deconwolf/`)
- [x] Identical results to matlab code.
- [x] Custom TIFF warning handle (not to overflow the console). -- disables messages from libtiff at default verbosity level.
- [x] Use tif tags to write meta data (also to transfer from input image). -- Transfers some metadata from ImageJ
- [x] Flag to save one image after each iteration (to facilitate setting the number of iterations). -- There is now the `--iterdump` option for that.
- [x] Put some checks to see what arrays are not properly aligned and causes crashes with `-ftree-vectorize` -- Should work now.
