Coordinate System
=================

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

Internally, data arrays are stored using the standard "matrix indexing" 
(i.e., ``ij``) scheme with (0,0) in the upper left corner, columns increasing
to the right, and rows increasing downwards:

.. image:: /_static/img/coords_2d_internal.png
    :width: 225px
    :align: center

Positive array values indicate more positive z (relative to the z-location of
the plane) and negative values indicate more negative z. 

.. note::

    Because Lentil stores physical (x, y) coordinates directly in arrays
    ordered as (r, c), you can't infer physical orientation from the array
    indices where data is stored. 

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
arrays in the standard "matrix indexing" (i.e., ``ij``) scheme with (0,0) in
the upper left corner, columns increasing to the right, and rows increasing 
downwards. Because this is also how Lentil stores its array data, plots made
using this default behavior will be correctly oriented with respect to 
Lentil's coordinate frame.

.. image:: /_static/img/matplotlib_origin.png
    :width: 500px
    :align: center

When using ``imshow()``, it can also be helpful to provide the ``extent``
parameter, which defines the physical coordinates of the array being plotted.
Extent ``should`` be given as ``(xmin, xmax, ymin, ymax)``:

.. plot::
    :context: reset
    :include-source:
    :scale: 50

    fig, ax = plt.subplots()
    ax.imshow(lentil.circle((256,256), 120), extent=[-0.5, 0.5, -0.5, 0.5])
