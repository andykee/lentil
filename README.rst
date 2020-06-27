Lentil
======
|build status| |coverage| |docs status| |pypi version|

Lentil is a Python library for modeling the imaging chain of an optical system.
It was originally developed at NASA's Jet Propulsion Lab by the Wavefront Sensing and
Control group (383E) to provide an easy to use framework for simulating point spread
functions of segmented aperture telescopes.

Lentil provides classes for representing optical elements with a simple interface for
including effects like wavefront error, radiometric properties, and various noise and
aberration sources. Lentil also provides numerical methods for performing Fraunhofer
(far-field) diffraction calculations. The collection of classes provided by Lentil can
be used to simulate imagery for a wide variety of optical systems.

Lentil is still under active development and new features continue to be added. Until
Lentil reaches version 1.0, the API is not guaranteed to be stable, but changes breaking
backwards compatibility will be noted.

Installing
----------
Install and update using `pip`_:

.. code-block:: text

    pip install lentil

Links
-----
* Documentation: https://lentil.readthedocs.io/
* Releases: https://pypi.org/project/lentil/
* Code: https://github.com/andykee/lentil/
* Issue tracker: https://github.com/andykee/lentil/issues/

.. _pip: https://pip.pypa.io/en/stable/quickstart/

.. |pypi version| image:: https://img.shields.io/pypi/v/lentil.svg
    :target: https://pypi.python.org/pypi/lentil

.. |build status| image:: https://travis-ci.com/andykee/lentil.svg?branch=master
    :target: https://travis-ci.com/andykee/lentil

.. |coverage| image:: https://coveralls.io/repos/github/andykee/lentil/badge.svg
    :target: https://coveralls.io/github/andykee/lentil

.. |docs status| image:: https://readthedocs.org/projects/lentil/badge/?version=latest
    :target: https://lentil.readthedocs.io/en/latest/?badge=latest
