.. _user.advanced.extend:

****************
Extending Lentil
****************

Subclassing existing planes
===========================
Lentil planes are designed to be subclassed, allowing users to easily modify
default behavior. A number of special methods are provided for easy operator
overloading:

Customizing attributes
----------------------
The following methods can be defined to customize plane attributes. These 
are commonly used when a plane attribute is stateful and needs to be 
recomputed at the time of attribute access. For example:

.. code-block:: python3

    class TipTiltMirror(lentil.Plane):
        def __init__(self):
            amplitude = lentil.circle((256,256),120)

            self.x = np.zeros(2)
            self._tt_basis = lentil.zernike_basis(amplitude, modes=[2,3])

            super().__init__(amplitude=amplitude)

        def __opd__(self):
            return np.einsum('ijk,i->jk', self._tt_basis, self.x)

.. py:function:: object.__amp__(self)

   Called when the plane's ``amplitude`` attribute is accessed. The 
   ``__amp__()`` method should return a Numpy array.

.. py:function:: object.__mask__(self)

   Called when the plane's ``mask`` attribute is accessed. The 
   ``__mask__()`` method should return a Numpy array.

.. py:function:: object.__opd__(self)

   Called when the plane's ``opd`` attribute is accessed. The 
   ``__opd__()`` method should return a Numpy array.

.. note::

    Numerical diffraction propagation calculations can repeatedly access 
    the opd, amplitude, and mask attributes. Because these attributes 
    remain fixed during a propagation, it is inefficient to recalculate 
    them each time the attribute is accessed. To mitigate this, the plane's
    ``freeze()`` method can be used to temporarily cache these attributes.
    For more information, see :ref:`user.performance.freeze`

Customizing how a plane interacts with a wavefront
--------------------------------------------------
.. py:function:: object.__mul__(self, wavefront)

   Called when a plane is multiplied with a wavefront.


Creating new planes
===================





.. _user.advanced.extend.tiltinterface:

.. Using TiltInterface
.. -------------------
