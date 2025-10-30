// documentation source as rst, then render as text ...
PyDoc_STRVAR(pydeconwolf__doc__,
             "Python interace to deconwolf\n"
             "--\n\n"
             "Important: This is still EXPERIMENTAL\n"
             "Please validate the results against the CLI version of deconwolf\n"

    );

PyDoc_STRVAR(imread__doc__,
             "imread(filename, verbose=0)\n"
             "--\n\n"
             "Read a 2D or 3D tif file and return as a fp32 numpy array\n"
             "\n"
             "Using libtiff in the background\n");


PyDoc_STRVAR(imwrite__doc__,
             "imwrite(array, filename, verbose=0)\n"
             "--\n\n"
             "Write a 2D or 3D fp32 numpy array as a tif file\n"
             "\n"
             "Using libtiff in the background\n");

PyDoc_STRVAR(bw__doc__,
             "gen_psf_bw(NA=None, ni=None, dx=None, dz=None, emission=None, nplanes=None, size=None, verbose=0)\n"
             "--\n\n"
             "Generate a 3D PSF using the \"Born and Wolf\" model"
             "\n\n"
             "Required Arguments\n"
             "------------------\n"
             "\t NA The numerical aperture of the objective, for example 1.45\n"
             "\t ni The refractive index of the immersion/air, for example 1.515\n"
             "\t dx The pixel size in [nm]\n"
             "\t dz The distance between planes [nm]\n"
             "\t emission The emission wave length [nm]\n"
             "\n"
             "\n"
             "Optional arguments\n"
             "==================\n"
             "\t nplanes The number of z planes to calculate [pixels]\n"
             "\t size    The lateral size of the output [pixels]\n"
             "\t verbose The verbosity level\n"
             "\t bq      Border handling, 0=circular, 1=compromise, 2=best\n"
             "\t gpu     Use GPU acceleration. Please see the notes below.\n"
             "*single*, **double**, ``code`` \n"
             "`Stack Overflow home <https://stackoverflow.com/>`_\n"
    );


PyDoc_STRVAR(deconvolve__doc__,
             "deconvolve(image, psf, iterations=0, gpu=False, verbose=0,  bq=2)\n"
             "--\n\n"
             "Deconvolve ...\n"
             "\n");

PyDoc_STRVAR(gsmooth__doc__,
             "gsmooth(image, sigmaxy=None, sigmaz=None, verbose=0)\n"
             "--\n\n"
             "gsmooth(image, sigmaxy=None, sigmaz=None, verbose=0)\n"
             "Normalized anisotropic Gaussian smoothing ...\n"
             );


