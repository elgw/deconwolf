Spot-based image alignment
==========================

.. note::
   New for version 0.4.4. More documentation planned.


Some microscopes will will have a shift (or translate) images slightly
different when switching imaging channels (most likely dues to mirror
alignments).

With deconwolf it is possible to estimate such shifts using detected
dots, check out the CLI interface with:

.. code:: shell

          dw align-dots --help

The program works by loading the dots from two files, we will denote
those sets as :math:`X = \{X_i | i = 1, .., m\}` and :math:`Y=\{Y_i | i = 1, ..., n\}`.

The method has the following steps:

1. A KD-Tree is constructed, using all offsets, :math:`s_{ij} = X_i - Y_j`
   where :math:`d(X_i,Y_j) < t_0`.

2. A KDE consturcted with the points in the KDTree, with a fixed
   Gaussian kernel determined by parameter X.

3. A grid search is carried out to find the approximate maxima of the KDE.

4. Based on the rough location from the grid search, the position of
   the maxima is refined to give at least a few decimals points.

Analyzing the output
--------------------

The first columns, contain `dx`, `dy`, and `dz`, which are the
components of :math:`\Delta`. They describe how to shift second image
to best overlap the first one. Units are in pixels.

The 4th column give :math:`kde(\Delta)`, the density of overlapping dots for the given
shift. Depending on the input data this can be a high or low
number. For an image with a single bead, a kde value of 1 would mean
that the two beads were correctly shifted to overlap. For more complex
images aim for a kde value of at least 3.

The `goodness` provides the contrast of the maxima,
i.e. :math:`kde(\Delta)/kde(\gamma)` where :math:`\gamma` is the
location having the 2nd higest `kde` value (:math:`d(\Delta,\gamma) >
xxx`).

Trouble shooting
----------------

Not that much should go wrong if this is used for the intended
data. If the shifts are large (the program does not find it) try
increasing the search radius for the initial pairing. This will cost
more processing time and consume more memory. Also, if there are
strong non-linear differences between the images, increasing the KDE
sigma might increase robustness.


Going further
-------------

It is possible to manually tune the magnification and the rotation of
the second point set but typically this kind of adjustments make
little sense to do prior to the alignment. It is better to use this
rigid shift alignment as the input for a more advanced correction
model, like what is described elsewhere [1]_. Ask for that to be be
implemented or send a pull request!

.. [1]
   Kozubek, M. and Matula, P. (2000), An efficient algorithm for measurement and correction of chromatic aberrations in fluorescence microscopy. Journal of Microscopy, 200: 206-217. https://doi.org/10.1046/j.1365-2818.2000.00754.x
