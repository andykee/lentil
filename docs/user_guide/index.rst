.. _user-guide:

**********
User Guide
**********

.. describing optical systems
..    planes
..    coordinate system
..    wavefront errors
.. modeling diffraction
.. radiometry
.. image sensors (incl. image artifacts)
.. model patterns
.. optimizing performance
.. external interfaces
.. algorithm verification

The User Guide provides documentation for all of Lentil's features and capabilities. It
is generally organized by topic area.

Brand new users should start with the :ref:`package-overview` and :ref:`tutorial`.

Detailed information on any specific class or method can be found in the :ref:`api`.

.. note::
    Many of the topics covered in this user guide and the API reference require some
    understanding of the underlying physics that govern the processes modeled
    by Lentil.

    This documentation does not attempt to provide a comprehensive reference to these
    topics but does provide links and references to supporting materials where applicable.


.. toctree::
    :maxdepth: 2

    optical_systems/index
    diffraction
    radiometry
    image_sensors/index
    ../patterns/index
    performance
    external/index
    verification/index
