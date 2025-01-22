Spot detection and measurements
===============================

.. note::
   New for version 0.4.4. More documentation planned.


Deconwolf includes a few option to detect diffraction limited dots (or
spots, or signals) via the module `dw dots`. This documentation gives
and overview of what the program does. Please see the man page for up
to date usage and exact command line arguments.

A tiny light emitter on a microscope slide will be seen as a blurry
spot under a microscope. The optical parameters decide how small it
can be. If it is as small as possible it will be denotes as a
*diffraction limited spot*.

The dot extraction pipeline consists of:

1. Detection spot candidates, using DoG or just by intensity combined
   with a local maxima finder. A big part of this is selecting a
   suitable filter to use.

2. Sub pixel localization using 3D gaussian blobs as a template.

3. Extra measurements. Most of the features that are measured are
   side-results of the fitting routine.

4. Exports to disk.

5. Filtering. This step is up to you to implement in you favorite
   programming language.


Spot size
---------

For wide field microscopy the optical parameters gives the spot
size. The Abbe resolution provides the full width half max (fwhm) in
the lateral plane (l) and along the axial direction (a). If
:math:`\Delta x` is the axial pixel sizes and :math:`\Delta z` is the
lateral pixel size, then:

.. math::

   \mbox{fwhm}_l = \frac{1}{\Delta x}\frac{\lambda}{2\mbox{NA}}

   \mbox{fwhm}_a = \frac{1}{\Delta z}\frac{2\lambda}{(\mbox{NA})^2}

And the equivalent size of a Gaussian bell is:

.. math::

   \sigma_{l} = \frac{fwhm_{l}}{2 \sqrt{2 \log(2)}}

   \sigma_{a} = \frac{fwhm_{a}}{2 \sqrt{2 \log(2)}}

which is used for the fitting. For the LoG filter the optimal size is
given by:

.. math::

   \sigma_{l, LoG} = \sigma_l \sqrt{2}

   \sigma_{a, LoG} = \sigma_a \sqrt{2}


These values are calculated automatically by `dw dots` if the relevant
optical parameters are supplied. Alternatively they can also be set
manually.

Dot candidates
--------------

In an fictional world where optics is perfect and there are no noise
sources, we could use the pixel values (after deconvolution) as direct
measurements on the number of photons.

In practice one could detect dots before deconvolution or after
deconvolution. If the images are not deconvolved, or only lightly
deconvolved there will be blur left in the image (due to the PSF). In
that case the pixel values are not very good estimators of the signal
brightness.

A common proxy for signal strength is the Laplacian of Gaussian (LoG)
feature. Intuitively the features is calculated by the difference
between what is in the centre of the filter (signal) to what is in the
surrounding pixels (background). I.e., it is a measure of the
elevation of the spot relative to the surroundings.

If there is a spot with a size that somewhat matches the LoG filter it
can be located as an extrema in the LoG filtered image.

Please note that spot finding method LoG+M will detect spots of any
size (although with decreasing amplitude as the spots size deviates
from the design size) and that it will also pick up noise and non-dot
structures in the image. Hence it is more or less mandatory to filter
the candidate spots by a few image features later on.


Fitting dots
------------

Fitting is performed using a multivariate Gaussian as the model,
restricting the shape to covariance matrices of the form:

.. math::

   \Sigma = \begin{pmatrix}
   \sigma_l^2 & 0 & 0 \\
   0 & \sigma_l^2 & 0 \\
   0 & 0 & \sigma_a^2
   \end{pmatrix}

The spot model is written as:

.. math::
   M(x) = c_0 + c_1G(x ; \mu, \Sigma)

If we collect all the parameters that we optimize over into
:math:`\Theta`, then we have:

.. math::
   \Theta = \mathop{\mathrm{arg\,min}}\nolimits -\sum \left(I(x)\log M(x) - M(x) \right)

Using the assumption of Poisson noise only.

Measurements
------------

The possible column in the output are as follows:

* *x*, *y*, *z* the coordinates of the spot before fitting. Note that
  this is 0-indexed and given in the units of pixels.

* *value* the filter value, i.e. the *LoG* value or the pixel value
  depending on what feature that was used to identify the dot
  candidates.

* *use* A simple (stupid) classification of each dot to say if it is
  relevant or not. This is based on an Otsu threshold on the `log(value)`

* *lat_circularity* the largest over the smallest eigenvalue of the
  covariance matrix of the local pixels.

Properties returned from the fitting routine:


* *f_bg* the background level of the spot.

* *f_signal_count* the total number of photons in the spot (assuming a
  quantum efficiency of 1).

* *f_peak_signal* the largest value of the fitted spot (note: can be
  higher than the maximum pixel value).

* *f_x*, *f_y*, *f_z* the fitted sup pixel location of the spot.

* *f_sigma_lateral*, *f_sigma_axial* the size, i.e., the estimated
  sigma value for the gaussian in the lateral plane and along the
  axial direction.

* *f_status* the status value returned from the fitting function.

* *f_error* :math:`-log(p)` of the spot.

* *f_corr* the correlation between the pixel values and the fitted
  model. This information is really useful to reject spot candidates
  that are found on line- and edge- like structures in images.
