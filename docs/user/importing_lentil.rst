.. _user.import:

********************************
Importing lentil and subpackages
********************************

Most of the functionality in Lentil resides in the standard library.
Once installed, Lentil's standard library can be imported as follows:

.. code-block:: python3

    >>> import lentil

If you'd prefer to use Lentil *en Français*, you can try

.. code-block:: python3

    >>> import lentil as le

API documentation for Lentil's standard library is available :ref:`here <api>`.

Lentil subpackages
==================
Lentil includes subpackges which provide domain-specific functionality. The 
available subpackages are described below:

==================================  ==========================
Subpackage                          Purpose
==================================  ==========================
:ref:`detector <api.detector>`      Focal plane modeling
:ref:`radiometry <api.radiometry>`  Working with spectral data
==================================  ==========================

Import a subpackage with

.. code-block:: python3

    >>> from lentil import detector
