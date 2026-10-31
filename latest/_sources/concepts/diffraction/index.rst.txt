Diffraction Propagation
=======================

.. toctree::
    :hidden:

    basics
    sampling
    output
    tilt
    segmented
    broadband
    multi
    near-field

Lentil models diffraction by propagating a |Wavefront| object between pupil
and image planes using Fourier transform-based propagators. The pages in this
section describe how the propagators work and how to use them. For worked,
end-to-end examples see the 
:doc:`diffraction propagation guides<../../guides/diffraction/index>`

:doc:`basics`
    An introduction to the propagation workflow and the available propagation
    functions.

:doc:`sampling`
    How input and output sampling set the accuracy and the usable field of
    view of a propagation, and how to choose ``pixelscale``, ``shape``, and
    ``oversample``.

:doc:`output`
    Computing only part of the output plane with ``prop_shape``, ``mask``,
    and ``shift``.

:doc:`tilt`
    How to handle tilts too large to sample in the pupil.

:doc:`segmented`
    Propagating segmented and sparse apertures efficiently.

:doc:`broadband`
    Modeling multi-wavelength propagations.

:doc:`multi`
    How to simulate multi-plane propagations with intermediate pupil and image
    planes.

:doc:`near-field`
    Propagating to planes that are not conjugate to a pupil or an image.