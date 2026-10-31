Propagation with segmented apertures
====================================
Creating a model of a segmented aperture doesn't require any special
treatment. A |Plane| works the same with segmented or sparse ``amplitude``,
``opd``, and ``mask`` attributes as with monolithic ones.

Defining a segmented plane
--------------------------
A plane's :attr:`~lentil.Plane.mask` can be a a 2-dimensional global mask or a
3-dimensional stack of segment masks:

.. plot:: _plot/segmask.py
    :scale: 50

With a segment mask, each segment is represented by a separate field cropped
to the segment mask's bounding box. Because wavefronts are monochromatic,
these separate fields cab ne propagate dindependently and added coherently as
needed.

Fitted tilt with segmented apertures
------------------------------------
Using :func:`~lentil.Plane.fit_filt` coupled with 
:func:`~lentil.propagate_dft` can significantly increase propagation
performance when modeling large tilts. This is not universally true since
tightly packed segmented apertures will result in overlapping segment mask
bounding boxes. 

.. EXAMPLE

.. note::

    A single segment's PSF is larger than the full aperture's by about the
    ratio of the aperture diameter to the segment diameter. Be sure to size
    ``prop_shape`` for the segment.