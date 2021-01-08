.. _user-guide.coordinate-system:

*****************
Coordinate System
*****************

.. |Pupil| replace:: :class:`~lentil.Pupil`
.. |Image| replace:: :class:`~lentil.Image`

Lentil adopts the widely used convention of aligning the z-axis along the direction
of light propagation through an optical system. By the right hand rule, it follows that
the remaining axes are oriented as shown in the figure below:

.. image:: /_static/img/coordinate_system.png
    :width: 500px
    :align: center


Lentil also adopts the right hand rule convention for rotations about the coordinate
system defined above. The following rotations are used:

* Rotations in +x rotate the yz plane counter-clockwise about the x-axis
* Rotations in +y rotate the xz plane counter-clockwise about the y-axis
* Rotations in +z rotate the xy plane counter-clockwise about the z-axis

When viewing a plane in 2D, the z-axis is assumed to come out of the screen with the
positive x-axis pointing to the right and the positive y-axis pointing up.

.. image:: /_static/img/coordinate_system_plot.png
    :width: 700px
    :align: center

Note that Matplotlib's ``imshow()`` method (and MATLAB's ``imagesc()`` method) places
the origin in the upper left corner of the plotted image by default. The result is that
the direction of y-axis is flipped relative to Lentil's coordinate system. This doesn't
necessarily present a problem as long as results are consistently plotted "incorrectly",
but to be completely correct (particularly when comparing model-generated images against
intuition or measured data) the origin should be located in the lower left corner of the
displayed image.
