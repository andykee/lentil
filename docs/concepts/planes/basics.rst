Plane basics
============

Overview
--------
Lentil represents optical systems using one or more planes that have some 
influence on a propagating wavefront. All planes in Lentil are described by a
common set of attributes. The most commonly used attributes are:

* ``amplitude`` - defines the electric field amplitude transmission through 
  the plane
* ``opd`` - defines the optical path difference that a wavefront experiences 
  when propagating through the plane
* ``mask`` - defines the binary mask over which the plane data is valid. If 
  ``mask`` is 2-dimensional, the plane is assumed to be monolithic. If 
  ``mask`` is 3-dimensional, the plane is assumed to be segmented with the 
  individual segment masks provided along the first dimension. If mask is not 
  provided, it is automatically created from the nonzero values in 
  ``amplitude``.

.. plot:: _plot/segmask.py
    :scale: 50

* ``pixelscale`` - defines the physical sampling of the above attributes. A 
  simple example of how to calculate the pixelscale for a discretely sampled 
  circular aperture is given below:

  .. image:: /_static/svg/pixelscale.svg
    :width: 450px
    :align: center

* ``ptype`` - specifies the plane type. Additional details are available in 
  the sections describing :ref:`user.wavefront.plane_wavefront` 
  and :ref:`user.wavefront.ptype_rules`

  Valid types are:
  
  ============= ============================
  ptype         Planes with this type
  ============= ============================
  ``none``      ``Plane``
  ``pupil``     ``Pupil``
  ``image``     ``Image``
  ============= ============================

Plane creation
--------------
Create a new plane with

.. plot::
    :context: reset
    :include-source:
    :nofigs:

    p = lentil.Plane(amplitude=lentil.circle((256,256), 120))

.. plot::
    :context: close-figs
    :scale: 50

    fig, ax = plt.subplots()
    ax.imshow(p.amplitude, cmap='gray', origin='lower',
              extent=[-128, 127, -128, 127])

Once a plane is defined, its attributes can be modified at any time:

.. plot::
    :context:
    :include-source:
    :nofigs:

    p = lentil.Plane(amplitude=lentil.circle((256,256), 120))
    p.opd = 2e-6 * lentil.zernike(p.mask, index=7)

.. plot::
    :context: close-figs
    :scale: 50

    fig, ax = plt.subplots()
    ax.imshow(p.opd, cmap='RdBu_r', origin='lower',
              extent=[-128, 127, -128, 127])


Resampling and rescaling
------------------------
Plane sampling can be modified using the :func:`~lentil.Plane.resample` or 
:func:`~lentil.Plane.rescale` methods:

* :func:`~lentil.Plane.resample` modifies plane sampling to the provided 
  ``pixelscale``
* :func:`~lentil.Plane.rescale` adjusts plane sampling by the given scale
  factor

Both methods use intrepolation to resample the ``amplitude``, ``opd``, and
``mask`` attributes and update the plane's ``pixelscale`` accordingly.


Combining planes
----------------


More information
----------------
* Descriptions of :doc:`available planes<types>` and
  :doc:`plane transforms<transforms>`
* Different approaches for
  :doc:`representing wavefront error<wavefront-error>`
* Details of the Plane-Wavefront interaction model
* How to modify default plane behavior and develop custom planes