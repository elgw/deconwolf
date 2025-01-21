Ilastic-like segmentation
=========================

Deconwolf is capable of segmenting 2D images.


Workflow
--------

1. Create training data


.. code::

   dw nuclei --init image1.tif image2.tif


This will create `image1.tif.a.png` and `image2.tif.a.png`.

Open these files with your favorite image editor and draw the nuclei
(or any objects of interest) in green, and non-nuclei (background etc)
in red.


2. Create a classifier

.. code::

   dw nuclei --fit model.trf *.tif

This command will tell dw to read the annotated images that you
created above and create a random forest classifier which will be
saved to disk as `model.trf`. For each image it will look of an
associated png image (with file extension `.a.png`). The tif file will
simply be ignored if there is no associated png file.

3. Classify

.. code::

   dw nuclei --predict model.trf file1.tif file2.tif
