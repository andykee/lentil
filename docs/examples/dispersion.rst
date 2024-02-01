.. _examples.dispersion:

*******************
Modeling dispersion
*******************

.. plot::
    :context: reset
    :include-source:
    :scale: 50

    import matplotlib.pyplot as plt
    import numpy as np
    import lentil

    amp = lentil.circle(shape=(256,256), radius=120)

    pupil = lentil.Pupil(amplitude=amp, opd=0, pixelscale=1/240,
                        focal_length=10)

    trace = [1, 0]  # y = x
    dispersion = [3e-4, 550e-9]  # center wavelength = 550 nm
    dispersive_element = lentil.DispersiveTilt(trace=trace, dispersion=dispersion)

    shape = (128,128)
    oversample = 3
    out = np.zeros((shape[0]*oversample, shape[1]*oversample))

    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(6,3))

    # dispersed PSF @ 475 nm
    w0 = lentil.Wavefront(475e-9)
    w1 = w0 * pupil
    w2 = w1 * dispersive_element
    w3 = lentil.propagate_dft(w2, shape=shape, prop_shape=(32, 32), 
                            pixelscale=5.5e-6, oversample=oversample)
    ax[0].imshow(w3.intensity, cmap='inferno')
    ax[0].set_title('$\lambda$ = 475 nm')
    ax[0].axis('off')

    # dispersed PSF @ 625 nm
    w0 = lentil.Wavefront(625e-9)
    w1 = w0 * pupil
    w2 = w1 * dispersive_element
    w3 = lentil.propagate_dft(w2, shape=shape, prop_shape=(32, 32), 
                            pixelscale=5.5e-6, oversample=oversample)
    ax[1].imshow(w3.intensity, cmap='inferno')
    ax[1].set_title('$\lambda$ = 625 nm')
    ax[1].axis('off')

    # broadband dispersion
    for wave in np.linspace(475, 625, 150):
        w0 = lentil.Wavefront(wavelength=wave*1e-9)
        w1 = w0 * pupil
        w2 = w1 * dispersive_element
        w3 = lentil.propagate_dft(w2, shape=shape, prop_shape=(32, 32), 
                                pixelscale=5.5e-6, oversample=oversample)

        out = w3.insert(out)


    ax[2].imshow(out, cmap='inferno')
    ax[2].set_title('$\lambda$ = [475 625] nm')
    ax[2].axis('off')