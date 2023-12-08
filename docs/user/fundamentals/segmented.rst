.. _user.fundamentals.segmented:

*************************
Segmented optical systems
*************************

Lentil works with segmented apertures out of the box.




Creating a model of a segmented aperture optical system in Lentil doesn't require any
special treatment. The |Plane| object works the same with sparse or
segmented amplitude, opd, and mask attributes as with monolithic ones.

That being said, it is advantageous from a performance point of view to supply a
3-dimensional `segment mask` when specifying a Plane's :attr:`~lentil.Plane.mask`
attribute rather than a flattened 2-dimensional `global mask` when working
with a segmented aperture, as depicted below:

.. plot:: _img/python/segmask.py
    :scale: 50

This modification is not necessary to achieve accurate propagations, but can
greatly improve performance. 




Since a wavefront is monochromatic, it can be broken into any arbitrary number of
fields and recombined coherently at a later time.
