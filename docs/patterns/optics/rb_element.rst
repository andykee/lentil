Developing an Element with 6DOF Rigid Body Control
==================================================
The example below assumes a ``dwdx`` influence function matrix is available.

.. code-block:: python3

    import lentil

    class RBElement(lentil.Plane):

        def __init__(self):
            # Establish the initial pose state as zeros
            self._x = np.zeros(6)

            # Define a WFE object which will return an OPD when called with a state.
            # dwdx is an influence function matrix mapping pose change (dx) to
            # wavefront (dw).
            self.rb_wfe = lentil.wfe.ParameterizedWFE(basis = dwdx)

        @property
        def opd(self):
            return self.rb_wfe(self.x, reshape=True)

        @property
        def x(self):
            return self._x

        def translate(self, x=0, y=0, z=0):
            self.perturb([0, 0, 0, x, y, z])

        def rotate(self, x=0, y=0, z=0):
            self.perturb([x, y, z, 0, 0, 0])

        def perturb(self, dx):
            self._x = self._x + dx


A slightly more involved model also bookkeeps rigid body actuator stroke, includes
actuator length limits, and introduces a small RBA stepping error every time the
element is perturbed. In addition to ``dwdx``, this example requires a ``dxdu``
influence function matrix.


.. code-block:: python3

    import numpy as np
    import lentil

    class RBElement(mo.Element):

        def __init__(self):
            # Establish the initial pose state as zeros
            self._x = np.zeros(6)

            # Define a WFE object which will return an OPD when called with a state.
            # dwdx is an influence function matrix mapping pose change (dx) to
            # wavefront (dw).
            self.rb_wfe = le.wfe.ParameterizedWFE(basis = dwdx)

            # Establish the initial rigid body actuator stroke as zeros
            self._u = np.zeros(6)

            # Define RBA stroke limits and stepping error
            self._umin = -0.5e-3
            self._umax = 0.5e-3
            self._uerror = 1e-6


            # dxdu is an influence function matrix which maps RBA stroke (du) to
            # pose change (dx)
            self._dxdu = dxdu

        @property
        def opd(self):
            return self.rb_wfe(self.x, reshape=True)

        @property
        def x(self):
            return self._x

        def translate(self, x=0, y=0, z=0):
            self.perturb([0, 0, 0, x, y, z])

        def rotate(self, x=0, y=0, z=0):
            self.perturb([x, y, z, 0, 0, 0])

        def perturb(self, dx):
            # Compute desired command in actuator space
            du = np.dot(np.linang.pinv(self._dxdu), dx)

            # Compute random stepping error
            uerror = np.random.uniform(low=-0.5*self._uerror, high=0.5*self._uerror, size=6)

            # The stepping error only applies to actuators involved in this move
            umove = np.zeros(6)
            umove[np.nonzero(du)] = 1

            # Include the stepping error on the moving actuators
            du = du + (uerror * umove)

            # Enforce stroke limits
            du[du > self._umax] = self._umax
            du[du < self._umin] = self._umin

            # Update u
            self._u = self._u + du

            # Update x
            self._x = self._x + np.dot(self._dxdu, du)
