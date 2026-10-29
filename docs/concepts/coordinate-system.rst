Coordinate System
=================

Overview
--------
Lentil adopts the widely used convention of aligning the z-axis along the direction
of light propagation through an optical system. For a right-handed coordinate system,
it follows that the remaining axes are oriented as shown in the figure below:

.. image:: /_static/img/coords_3d.png
    :width: 500px
    :align: center

When viewing a plane in 2D, the z-axis comes out of the screen with the
positive x-axis pointing to the right and the positive y-axis pointing up:

.. image:: /_static/img/coords_2d.png
    :width: 225px
    :align: center

Angles are measured from +x increasing counterclockwise towards +y. Positive
array values indicate more positive z (relative to the z-location of the plane)
and negative values indicate more negative z. 

Array storage
-------------
Lentil stores data in standard NumPy arrays with (0, 0) in the first row and
column. The physical x-axis runs along columns with +x increasing with column
index and the physical y-axis runs along rows with +y increasing with row
index, measured from the center pixel at ``(nrows//2, ncols//2)``. This is the
same convention used by `numpy.meshgrid`.

Functions that accept or return coordinates take an ``indexing`` argument that
controls only the order of those values: ``(x, y)`` for ``'xy'`` and
``(r, c)`` for ``'ij'``. Indexing does not change how data is stored
internally.

.. image:: /_static/img/coords_2d_internal.png
    :width: 225px
    :align: center

Wavefront error sign convention
-------------------------------
The sign convention used for representing wavefront error follows directly
from the convention chosen for harmonic time dependence: 
:math:`\exp\left(\pm i\omega t\right)`. Lentil adopts the negative sign 
convention :math:`\exp\left(-i\omega t\right)` used by popular physical optics
textbooks, which has the following implications:

* Positive optical path difference and wavefront error correspond with a
  longer path length
* Positive focus error corresponds with a longer focal length
* Rotations (i.e. tilts) follows the right hand rule
* Forward propagation is accomplished using the (forward) Fourier transform

A more in depth treatment of this topic is 
:doc:`available in the technical notes<../resources/tech-notes/coords>`.

Plotting
--------
Matplotlib's ``imshow()`` method (and MATLAB's ``imagesc()`` method) place the
origin in the upper left corner of the plotted image by default. This displays
arrays in the standard "matrix indexing" (i.e., ``'ij'``) scheme with (0,0) in
the upper left corner, columns increasing to the right, and rows increasing 
downwards. Because Lentil stores data with +y increasing with row index,
arrays plotted using the Matplotlib default will be flipped top to bottom.
Plotting with ``origin='lower'`` displays arrays in their correct (physical)
orientation:

.. code:: python

   ax.imshow(opd, origin='lower', extent=[-0.5, 0.5, -0.5, 0.5])


.. image:: /_static/img/matplotlib_origin.png
    :width: 500px
    :align: center
