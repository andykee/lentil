.. _user.fundamentals.introduction:

************
Introduction
************

.. plot:: user/fundamentals/plots/intro_banner.py
    :scale: 48

Lentil is a Python library for modeling the imaging chain of an optical 
system. It was originally developed at NASA's Jet Propulsion Lab by the 
Wavefront Sensing and Control group (383E) to provide an easy to use framework 
for simulating point spread functions of segmented aperture telescopes.

Lentil provides a framework for defining optical planes and combining them to 
model diffraction and other imaging effects. It can be used to create models 
with varying levels of complexity; from simple prototypes to extremely high-
fidelity models of ground or space telescopes.

Lentil uses one or more |Plane| objects to represent an optical system. A
|Wavefront| object 
