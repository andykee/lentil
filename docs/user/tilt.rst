.. _user.tilt:

****************************************
Diffraction propagation with large tilts
****************************************


One of lentil's unique features is its hybrid approach to handling tilt.

The plane's :func:`~lentil.Plane.fit_tilt` method performs a least squares fit 
to estimate and remove tilt from the :attr:`~lentil.Plane.opd` attribute. The 
removed tilt is accounted for by appending an equivalent :class:`~lentil.Tilt` 
object to the plane's :attr:`~lentil.Plane.tilt` attribute. The optical effect 
of the tilt is automatically applied during a propagation step.





.. image:: /_static/img/propagate_tilt_phase.png
    :width: 450px
    :align: center

.. image:: /_static/img/propagate_tilt_phase_wrap.png
    :width: 650px
    :align: center

.. image:: /_static/img/propagate_tilt_angle.png
    :width: 600px
    :align: center

.. image:: /_static/img/propagate_tilt_angle_steps.png
    :width: 600px
    :align: center