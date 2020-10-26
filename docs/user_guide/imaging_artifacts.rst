*****************
Imaging Artifacts
*****************

.. currentmodule:: lentil

Lentil's :ref:`convolvable module<api-convolvable>` contains methods for applying various imaging 
artifacts via convolution. 

.. note::

    To ensure accuracy and avoid introducing ailiasing artifacts, the input data should 
    be Nyquist sampled.

Smear
=====
Smear is used to represent image motion with a relatively low temporal frequency relative
to integration time. The motion occurs in a slowly varying or fixed direction over one
integration time. Lentil represents smear as a directional blur over some distance (or 
number of pixels):



Jitter
======
Jitter is used to represent image motion with a relatively high temporal frequency relative
to integration time. The motion occurs in all directions randomly during an integration time.
Lentil represents jitter using a Gaussian blur:



Pixel MTF
=========
A focal plane array samples a continuous light field to produce a digital image. Because
Lentil models diffraction numerically by propagating a finite set of points through an
optical system, the discretely sampled image plane intensity must be convolved with the
pixel's aperture function to accurately represent the intensity signal seen by each
pixel. Lentil's :func:`convolvable.pixel` method implements this convolution. 
After convolving the image plane intensity with the pixel MTF, the data should be 
resampled to native detector sampling using :func:`util.rescale`. The
:func:`detector.pixelate` method combines the convolution and resampling operations
into a single method.
