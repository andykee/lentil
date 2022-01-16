.. _user_guide.coordinate_system:

*****************
Coordinate system
*****************

.. |Pupil| replace:: :class:`~lentil.Pupil`
.. |Image| replace:: :class:`~lentil.Image`

Lentil adopts the widely used convention of aligning the z-axis along the direction
of light propagation through an optical system. By the right hand rule, it follows that
the remaining axes are oriented as shown in the figure below:


.. image:: /_static/img/coordinate_system.png
    :width: 500px
    :align: center

When viewing a plane in 2D, the z-axis comes out of the screen with the
positive x-axis pointing to the right and the positive y-axis pointing up.

Lentil also adopts the right hand rule convention for rotations about the coordinate
system defined above:

* Rotations in +x rotate the yz plane counter-clockwise about the x-axis
* Rotations in +y rotate the xz plane counter-clockwise about the y-axis
* Rotations in +z rotate the xy plane counter-clockwise about the z-axis

.. plot::
    :scale: 50

    import matplotlib.pyplot as plt
    import numpy as np
    import lentil

    import matplotlib as mpl
    mpl.rcParams['figure.figsize'] = (3.5, 3.5)

    amp = lentil.circle((256, 256), 120)
    x_tilt = 8e-6 * lentil.zernike(amp, 3) # +x tilt
    y_tilt = 8e-6 * lentil.zernike(amp, 2) # +y tilt

    py = lentil.Pupil(focal_length=10, pixelscale=1/240, amplitude=amp, phase=y_tilt)
    wy = lentil.Wavefront(650e-9)
    wy *= py
    wy = wy.propagate_image(pixelscale=5e-6, npix=256, oversample=5)

    px = lentil.Pupil(focal_length=10, pixelscale=1/240, amplitude=amp, phase=x_tilt)
    wx = lentil.Wavefront(650e-9)
    wx *= px
    wx = wx.propagate_image(pixelscale=5e-6, npix=256, oversample=5)

    plt.subplot(2, 2, 1)
    plt.imshow(x_tilt, origin='lower')
    plt.title('Pupil plane [$+R_x$]')
    plt.xticks(np.linspace(0, 256, 5), labels=np.linspace(-1, 1, 5))
    plt.yticks(np.linspace(0, 256, 5), labels=np.linspace(-1, 1, 5))
    #plt.grid()

    plt.subplot(2, 2, 2)
    plt.imshow(y_tilt, origin='lower')
    plt.title('Pupil plane [$+R_y$]')
    plt.xticks(np.linspace(0, 256, 5), labels=np.linspace(-1, 1, 5))
    plt.yticks(np.linspace(0, 256, 5), labels=np.linspace(-1, 1, 5))
    #plt.grid()

    plt.subplot(2, 2, 3)
    plt.imshow(wx.intensity**0.2, origin='lower')
    plt.title('Image plane [$+R_x$]')
    plt.xticks(np.linspace(0, 256*5, 5), labels=np.linspace(-1, 1, 5))
    plt.yticks(np.linspace(0, 256*5, 5), labels=np.linspace(-1, 1, 5))
    plt.subplot(2, 2, 4)
    plt.imshow(wy.intensity**0.2, origin='lower')
    plt.title('Image plane [$+R_y$]')
    plt.xticks(np.linspace(0, 256*5, 5), labels=np.linspace(-1, 1, 5))
    plt.yticks(np.linspace(0, 256*5, 5), labels=np.linspace(-1, 1, 5))

    plt.tight_layout()

.. note::

    Matplotlib's ``imshow()`` method (and MATLAB's ``imagesc()`` method) place
    the origin in the upper left corner of the plotted image by default. This presents
    arrays in the standard (row, column) ordering. The result is that the direction of
    y-axis is flipped relative to Lentil's coordinate system. This doesn't necessarily
    present a problem as long as results are consistently plotted "incorrectly", but
    to be completely correct (particularly when comparing model-generated images against
    intuition or measured data) the origin should be located in the lower left corner
    of the displayed image.

    .. image:: /_static/img/coordinate_system_plot.png
        :width: 700px
        :align: center
