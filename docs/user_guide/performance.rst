**********************
Optimizing Performance
**********************

.. _caching:

Caching
=======
While the computation that occurs when accessing the ``opd`` attribute in the
example above isn't terribly expensive, that computation (along with many other
similar computations) may be performed tens, hundreds, or even thousands of times
in the process of performing a single propagation.

Lentil provides a caching mechanism to reduce the number of times these
computations are performed. The :class:`~lentil.modeltools.cached_property`
decorator is a drop-in replacement for Python's ``property`` decorator and ensures
that decorated functions have their values computed and cached once at the beginning
of a propagation. The cache is automatically cleared when the propagation is
complete.

.. _performance-image-simulation:

Image simulation
================



Speeding up the FFT
===================


Multiprocessing
===============

Other Performance Tweaks
========================

Faster photon to electron (quantum efficiency) calculations
-----------------------------------------------------------
