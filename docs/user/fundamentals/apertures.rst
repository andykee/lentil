.. _user.fundamentals.apertures:

*******************
Apertures and masks
*******************

Lentil expects apertures and masks to be described in a discretely sampled 
array-like format - most commonly in the form of a NumPy array. Portions of
the array with values greater than one are assued to be transmissive while 
portions of the array that are zero are assumed to be blocking.

Aperture or mask shapes can be loaded from a file, created manually, or 
constructed using one or more of Lentil's functions that draw common shapes
in an array.

Core shapes
===========
Lentil provides a number of functions for drawing basic shapes in arrays.
Multiple arrays can be combined to create more complicated shapes.

The core shape functions support shifting the center of the shape arbitrarily
relative to the center of the array. The shift is defined in terms of 
(row, col).

.. plot:: user/fundamentals/plots/shape_shift.py
    :scale: 50


Basic antialiasing is also supported (and is enabled by default):

.. plot:: user/fundamentals/plots/shape_antialias.py
    :scale: 50

Circle
------


Rectangle
---------

Hexagon
-------


Hex segmented apertures
=======================


