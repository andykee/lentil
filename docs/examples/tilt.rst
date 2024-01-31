.. _examples.tilt:

**************************************
Diffraction simulation with large tilt
**************************************

.. plot::
    :context: reset
    :include-source:
    :scale: 50

    import lentil
    import numpy as np
    import matplotlib.pyplot as plt

    amp = lentil.circle((256,256), 120)
    opd = lentil.zernike(amp, 2)*15e-6

    p = lentil.Pupil(amplitude=amp, opd=opd, focal_length=10, pixelscale=1/240)
    w = lentil.Wavefront(500e-9)
    w *= p
    w = lentil.propagate_dft(w, pixelscale=5e-6, shape=(256, 256))

    psf = w.intensity/np.max(w.intensity)

    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    ax.imshow(psf**0.1, cmap='inferno', vmin=0.15)
    ax.set_title('Large tilt exposes periodic wraparound')

.. plot::
    :context: close-figs
    :include-source:
    :scale: 50

    p_geom_tilt = p.fit_tilt()

    w = lentil.Wavefront(500e-9)
    w *= p_geom_tilt
    w = lentil.propagate_dft(w, pixelscale=5e-6, shape=(256, 256))

    psf = w.intensity/np.max(w.intensity)

    fig, ax = plt.subplots(figsize=(2.5, 2.5))
    ax.imshow(psf**0.1, cmap='inferno', vmin=0.15)
    ax.set_title('Geometric tilt handling eliminates wraparound')