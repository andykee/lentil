.. _package-overview:

****************
Package overview
****************

Lentil is a Python library for modeling the imaging chain of an optical system. It was 
originally developed at NASA's Jet Propulsion Lab by the Wavefront Sensing and Control 
group (383E) to provide an easy to use framework for simulating point spread functions 
of segmented aperture telescopes.

Lentil provides a framework for defining optical planes and combining them to model 
diffraction and other imaging effects. It can be used to create models with varying 
levels of complexity; from simple prototypes to extremely high-fidelity models of ground 
or space telescopes. Lentil was developed to be:

* User friendly - a simple, consistent interface optimized for common use cases.
* Modular - models are made by combining configurable building blocks. It is easy to 
  replace or update individual elements without impacting the rest of the model.
* Easy to extend - core objects are designed to be subclassed, modified, extended, 
  broken apart, and put back together in different ways.

Lentil organization
===================
Lentil is organized into two parts: a standard library and a few subpackages 
providing domain-specific tools. The features of the standard library and each
of the subpackages is described in the table below:

============================================== ===============================================
Namespace                                      Purpose
============================================== ===============================================
:ref:`lentil <api>`                            Standard library
:ref:`lentil.detector <api.detector>`          Model focal planes
:ref:`lentil.radiometry <api.radiometry>`      Work with spectral data
============================================== ===============================================

Lentil's standard library can be imported as follows:

.. code-block:: python3

    >>> import lentil

If you'd prefer to use Lentil *en français*, you can try

.. code-block:: python3

    >>> import lentil as le

Public subpackages are automatically imported with Lentil.


Getting help
============
The best place to ask for help on subjects not covered in this documentation or suggest new 
features/ideas is by opening a ticket on `Github <https://github.com/andykee/lentil/issues>`__.

License
=======
Copyright (c) 2021, California Institute of Technology ("Caltech"). U.S. Government sponsorship acknowledged.

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.





