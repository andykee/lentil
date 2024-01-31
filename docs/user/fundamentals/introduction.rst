.. _user.fundamentals.introduction:

************
Introduction
************

.. plot:: user/fundamentals/plots/intro_banner.py

Lentil is a Python library for modeling the imaging chain of an optical 
system. It was originally developed at NASA's Jet Propulsion Lab by the 
Wavefront Sensing and Control group (383E) to provide an easy to use framework 
for simulating point spread functions of segmented aperture telescopes.

Lentil provides a framework for defining optical planes and combining them to 
model diffraction and other imaging effects. It can be used to create models 
with varying levels of complexity - from simple prototypes to extremely 
high-fidelity models of ground or space telescopes.

Lentil uses one or more :ref:`Plane<user.fundamentals.planes>` objects to 
numericlly model an optical system. A 
:ref:`Wavefront<user.fundamentals.wavefront>` object representing a 
discretely sampled monochromatic plane wave is then propagated from plane 
to plane using one of the available :ref:`numerical diffraction propagation 
routines<user.fundamentals.diffraction>`. The resulting complex field can 
be visualized, analyzed, or further propagated through an optical model.
Finally, Lentil includes a number of additional tools for :ref:`representing
common imaging artifacts<user.fundamentals.artifacts>`, :ref:`modeling focal 
plane arrays<user.fundamentals.image_sensors>` including many common noise 
sources, and :ref:`working with radiometric quantities
<user.fundamentals.radiometry>`.
