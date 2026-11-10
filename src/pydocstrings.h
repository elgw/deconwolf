// this file was generated from pydocstrings.txt using genpydocstrings.py
PyDoc_STRVAR(pydeconwolf__doc__,
"Python interface to parts and pieces of deconwolf\n"
"--\n\n"
"EXPERIMENTAL! Please validate the results against the CLI version of deconwolf\n"
);
PyDoc_STRVAR(deconvolve__doc__,
"deconvolve(image, psf, iterations=0, gpu=False, verbose=0, bq=2)\n"
"--\n\n"
"Required Arguments\n"
"------------------\n"
"   image : The 3D images to deconvolve\n"
"   psf : the 3D PSF to used\n"
"   iterations : The number of iterations to use\n"
"\n"
"Optional arguments\n"
"------------------\n"
"   gpu : Use GPU or not\n"
"   verbose: Verbosity level\n"
"   bq : The border handling method. 0=periodic, 1=compromise, 2=best\n"
);
PyDoc_STRVAR(imread__doc__,
"imread(filename, verbose=0)\n"
"--\n\n"
"Read a 2D or 3D tif file and return as a fp32 numpy array\n"
"Using libtiff in the background\n"
);
PyDoc_STRVAR(imwrite__doc__,
"imwrite(array, filename, verbose=0)\n"
"--\n\n"
"Write a 2D or 3D fp32 numpy array as a tif file\n"
"Using libtiff in the background\n"
);
PyDoc_STRVAR(qalign__doc__,
"qalign(X, Y, rel_error, localization_sigma, npoint, verbose)\n"
"--\n\n"
"Calls the low level routine qalign which is a PART of dw align_dots\n"
"\n"
"Required Arguments\n"
"------------------\n"
"   X A [n x 3] numpy array of spots in the first image\n"
"   Y A [n x 3] numpy array of spots in the first image\n"
"\n"
"Optional arguments\n"
"------------------\n"
"   rel_error           The upper bound on relative errors on the spot pairs\n"
"   localization_sigma  The localization precision [pixels]\n"
"   npoint              Most number of points to use from X and Y\n"
"   verbose             The verbosity level\n"
);
PyDoc_STRVAR(gsmooth__doc__,
"gsmooth(image, sigmaxy=None, sigmaz=None, verbose=0)\n"
"--\n\n"
"Normalized anisotropic Gaussian smoothing ...\n"
);
PyDoc_STRVAR(bw__doc__,
"gen_psf_bw(NA=None, ni=None, dx=None, dz=None, emission=None, nplanes=None, size=None, verbose=0)\n"
"--\n\n"
"Generate a 3D PSF using the \"Born and Wolf\" model\n"
"\n"
"Required Arguments\n"
"------------------\n"
"   NA The numerical aperture of the objective, for example 1.45\n"
"   ni The refractive index of the immersion/air, for example 1.515\n"
"   dx The pixel size in [nm]\n"
"   dz The distance between planes [nm]\n"
"   emission The emission wave length [nm]\n"
"\n"
"Optional arguments\n"
"------------------\n"
"   nplanes The number of z planes to calculate [pixels]\n"
"   size    The lateral size of the output [pixels]\n"
"   verbose The verbosity level\n"
"   bq      Border handling, 0=circular, 1=compromise, 2=best\n"
"   gpu     Use GPU acceleration. Please see the notes below.\n"
);
PyDoc_STRVAR(std__doc__,
"std(X)\n"
"--\n\n"
"standard deviation\n"
);
PyDoc_STRVAR(std_masked__doc__,
"std_masked(X, mask)\n"
"--\n\n"
"standard deviation using a u8 mask to select values to use.\n"
);
