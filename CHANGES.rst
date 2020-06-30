Changes
=======

v0.2.0
------
Released June 29, 2020

* Collapse ``Detector`` and ``Image`` planes into single ``Image`` plane. The pupil to
  image plane propagation method is now chosen based on whether the ``Image`` plane has
  a defined ``pixelscale`` (propagate via matrix triple product DFT) or if
  ``pixelscale`` is None (propagate via FFT - eventually). ``Detector`` class has been
  deprecated. `#5`_
* Deprecate ``FPA`` and ``BayerFPA``. Some functionality has been retained but converted
  to atomic functions in the :ref:`detector<api-detector>` module. `#6`_
* Completely rework the contents of the :ref:`detector<api-detector>` module. All
  objects have been deprecated. Some functionality has been retained but converted to
  atomic functions instead. `#6`_
* Deprecate ``util.coordinates``
* Change the way ``Rotate`` angle is interpreted to behave more intuitively
* A number of small bugfixes and enhancements
* Updated documentation
* More unit tests

.. _#5: https://github.com/andykee/lentil/issues/5
.. _#6: https://github.com/andykee/lentil/issues/6

v0.1.1
------
Released June 21, 2020

* Update ``propagate`` to support :class:`~lentil.Tilt` planes `#1`_
* Streamline the innards of :func:`~lentil.propagate`
* Update :func:`lentil.wfe.power_spectrum` to return phases with a slightly more correct
  RMS
* Remove unused code
* Increase unit testing coverage
* Set up Travis CI, Coveralls

.. _#1: https://github.com/andykee/lentil/issues/1

v0.1.0
------
Released June 12, 2020

* Initial public release
