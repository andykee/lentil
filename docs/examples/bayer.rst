.. _examples.bayer:

*************************
Modeling a Bayer detector
*************************

.. plot::
    :context: reset
    :include-source:
    :scale: 50

    import lentil
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.signal

    gaussian = scipy.signal.windows.gaussian(250, 25)
    gaussian = gaussian[25:-25]
    qe_wave = np.arange(350, 750)*1e-9

    qe_red = np.zeros(qe_wave.size)
    qe_green = np.zeros(qe_wave.size)
    qe_blue = np.zeros(qe_wave.size)

    qe_blue[10:210] = gaussian
    qe_green[90:290] = gaussian
    qe_red[170:370] = gaussian


    plt.plot(qe_wave, qe_blue, 'b', label='Blue')
    plt.plot(qe_wave, qe_green, 'g', label='Green')
    plt.plot(qe_wave, qe_red, 'r', label='Red')
    plt.legend()
    plt.title('Quantum Efficiency')
    plt.xlabel('Wavelength [nm]')
    plt.grid()


.. plot::
    :context: close-figs
    :include-source:
    :scale: 50

    amp = lentil.circle((256,256), 120)
    amp -= lentil.circle(shape=(256,256), radius=40)
    for angle in (0, 90, 180, 270):
        amp *= lentil.spider((256,256), 2, angle)
        
    opd = lentil.zernike(amp, 4)*1e-6

    psf = []

    p = lentil.Pupil(amplitude=amp, opd=opd, focal_length=10, pixelscale=1/240)

    for wave in np.arange(350, 750):
        w = lentil.Wavefront(wave*1e-9)
        w *= p
        w = lentil.propagate_dft(w, pixelscale=5e-6, shape=(64, 64), 
                                oversample=2)
        psf.append(w.intensity)

    psf = np.asarray(psf)
    psf = lentil.rebin(psf, 2)

    psf_flat = np.sum(psf, axis=0)

    fig, ax = plt.subplots(figsize=(2.5,2.5))
    ax.imshow(psf_flat, cmap='inferno')
    ax.set_title('Broadband PSF')
    plt.axis('off')


.. plot::
    :context: close-figs
    :include-source:

    img = lentil.detector.collect_charge_bayer(psf, qe_wave, qe_red, qe_green,
                                            qe_blue, 'BGGR')

    red_kernel = np.array([[0,0],[0,1]])
    green_kernel = np.array([[0,1],[1,0]])
    blue_kernel = np.array([[1,0],[0,0]])

    red_mosaic = np.tile(red_kernel, (img.shape[0]//2, img.shape[1]//2))
    green_mosaic = np.tile(green_kernel, (img.shape[0]//2, img.shape[1]//2))
    blue_mosaic = np.tile(blue_kernel, (img.shape[0]//2, img.shape[1]//2))

    img = img/np.max(img)

    r = red_mosaic * img
    g = green_mosaic * img
    b = blue_mosaic * img

    fig, ax = plt.subplots(nrows=1, ncols=4, figsize=(7,3))

    ax[0].imshow(r, cmap='gray')
    ax[0].set_title('Red')
    ax[0].axis('off')

    ax[1].imshow(g, cmap='gray')
    ax[1].set_title('Green')
    ax[1].axis('off')

    ax[2].imshow(b, cmap='gray')
    ax[2].set_title('Blue')
    ax[2].axis('off')

    ax[3].imshow(img, cmap='gray')
    ax[3].set_title('Full Bayer image')
    ax[3].axis('off')
