.. _user.fundamentals.apertures:

*******************
Apertures and masks
*******************

Lentil expects apertures and masks to be described in a discretely sampled 
array-like format - most commonly in the form of a NumPy array. Portions of
the array with values greater than one are assued to be transmissive while 
portions of the array that are zero are assumed to be opaque.

Aperture or mask shapes can be loaded from a file, created manually, or 
constructed using one or more of Lentil's functions for drawing common 
shapes in an array.

Core shapes
===========
Lentil provides a number of functions for drawing basic shapes in arrays.
Multiple arrays can be combined to create more complicated shapes.

.. plot:: user/fundamentals/plots/apertures.py
    :scale: 50

Basic antialiasing is supported (and is enabled by default):

.. plot:: user/fundamentals/plots/shape_antialias.py
    :scale: 50

The core shape functions also support shifting the center of the shape 
arbitrarily relative to the center of the array. The shift is defined 
in terms of (row, col).

.. plot:: user/fundamentals/plots/shape_shift.py
    :scale: 50


Circle
------
Draw a circle with :func:`~lentil.circle`:

.. plot::
    :include-source:
    :scale: 50

    >>> circ = lentil.circle(shape=(256,256), radius=100)
    >>> plt.imshow(circ, cmap='gray')

Rectangle
---------
Draw a rectangle with :func:`~lentil.rectangle`:

.. plot::
    :include-source:
    :scale: 50

    >>> rect = lentil.rectangle(shape=(256,256), width=200, height=100)
    >>> plt.imshow(rect, cmap='gray')

Hexagon
-------
Draw a hexagon with :func:`~lentil.hexagon`:

.. plot::
    :include-source:
    :scale: 50

    >>> hex = lentil.hexagon(shape=(256,256), radius=75)
    >>> plt.imshow(hex, cmap='gray')

Spider
------
Draw a spider with :func:`~lentil.spider`:

.. plot::
    :include-source:
    :scale: 50

    >>> spider = lentil.spider(shape=(256,256), width=3, angle=30)
    >>> plt.imshow(spider, cmap='gray')


Composite apertures
===================
Because the core shape functions simply return NumPy arrays, it is possible
to create more complicated shapes by combining multiple arrays together. Below
is an example of how to draw the Hubble mask:

.. plot::
    :include-source:
    :scale: 50

    # dimensions from Tiny Tim (Krist & Hook 2011)
    outer_diam = 2.4
    central_obsc = .33
    spider = 0.0264
    mount_diam = 0.13
    mount_dist = 0.8921

    shape = (256, 256)
    pixelscale = 0.01

    npix_outer = outer_diam/pixelscale
    npix_inner = (outer_diam * central_obsc)/pixelscale
    npix_spider = spider/pixelscale
    npix_mount = mount_diam/pixelscale
    npix_mount_dist = npix_outer * mount_dist / 2

    # primary mirror
    hubble_outer = lentil.circle(shape, radius=npix_outer/2)
    hubble_inner = lentil.circle(shape, radius=npix_inner/2)
    hubble = hubble_outer - hubble_inner

    # secondary spiders
    for angle in (45, 135, 225, 315):
        hubble *= lentil.spider(shape, width=npix_spider, angle=angle)

    # primary mirror mounting pads
    for angle in (75, 195, 315):
        mount_shift = (npix_mount_dist * -np.sin(np.deg2rad(angle)),
                        npix_mount_dist * np.cos(np.deg2rad(angle)))
        hubble *= 1 - lentil.circle(shape, npix_mount/2, shift=mount_shift)

    plt.imshow(hubble, cmap='gray')


Hex segmented apertures
=======================
The :func:`~lentil.hex_segments` function constructs an aperture made up of
a number of concentric rings of hexagonal segments. An example showing how
to construct the James Webb Space Telescope aperture is below:

.. plot::
    :include-source:
    :scale: 50

    # dimensions from WebbPSF (STScI)
    segment_diam = 1.524
    segment_gap = 0.0075
    spider = 0.083

    # pixelscale is selected to provide at least 2 samples across
    # the smallest feature (segment gap)
    pixelscale = 0.003
    
    npix_seg = segment_diam/pixelscale
    npix_gap = segment_gap/pixelscale
    npix_spider = spider/pixelscale

    jwst = lentil.hex_segments(rings=2, seg_radius=npix_seg/2, 
                               seg_gap=npix_gap, flatten=True)

    # secondary spiders
    for angle in (90, 240, 300):
        jwst *= lentil.spider(jwst.shape, width=npix_spider, angle=angle)
    jwst *= 1-lentil.rectangle(jwst.shape, width=36, height=60, shift=(-750,0))
    jwst *= 1-lentil.rectangle(jwst.shape, width=36, height=20, shift=(-825,0))

    plt.imshow(jwst, cmap='gray')





