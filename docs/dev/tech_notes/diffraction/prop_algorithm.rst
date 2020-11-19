*****************************
General propagation algorithm
*****************************

.. currentmodule:: lentil

The general propagation algorithm implemented by :class:`prop.Propagate` is as follows:

.. code::

    for each plane in planes
        cache plane amplitude, phase, ptt_vector



    for each plane in planes
        clear plane amplitude, phase, ptt_vector cache