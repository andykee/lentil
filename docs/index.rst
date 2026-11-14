.. _docs:

.. toctree::
    :hidden:

    getting-started/index
    concepts/index
    guides/index
    resources/index
    api/index

Lentil
======
.. plot:: _plot/intro_banner.py
  :width: 100%

Lentil is an Python package for developing high-performance diffraction
simulations. Lentil was originally developed at NASA's Jet Propulsion Library
and is now available as open-source software. Lentil makes it easy to

* Model :doc:`optical systems<concepts/planes/basics>` using a suite of
  pre-built :doc:`planes<concepts/planes/types>` and
  :doc:`optical transformations<concepts/planes/transforms>`
* Numerically simulate :doc:`diffraction propagation<concepts/diffraction/basics>`
  including :doc:`multi-plane propagation<concepts/diffraction/multi>`
* Easily represent optical systems with
  :doc:`segmented apertures<concepts/diffraction/segmented>`
* Work with :doc:`spectral data<concepts/radiometry>` to accurately represent
  the radiometry of an imaging system
* Model common :doc:`focal plane effects<concepts/image-sensors>` and
  :doc:`dynamic imaging artifacts<concepts/imaging-artifacts>`

.. note::
    
    Lentil is under active development and new features continue to be added. 
    Until Lentil reaches version 1.0, the API is not guaranteed to be stable, but
    changes breaking backwards compatibility will be noted in the 
    :doc:`release notes <resources/release>`.


Installation
------------
Install Lentil with pip:

.. code:: bash

    pip install lentil

Additional details are available on the 
:doc:`installation<getting-started/installation>` page.

Learn more
----------
See the :doc:`first steps<getting-started/first-steps>` or jump straight to
the :doc:`guides<guides/index>` to start using Lentil.

