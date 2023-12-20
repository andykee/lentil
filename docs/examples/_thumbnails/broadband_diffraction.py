import matplotlib.pyplot as plt
import numpy as np
import lentil

amp = lentil.circle(shape=(256,256), radius=120)
coef = [0, 0, 0, 150e-9, -50e-9, 0, 0, 0, 35e-9]
opd = lentil.zernike_compose(mask=amp, coeffs=coef)
pupil = lentil.Pupil(amplitude=amp, opd=opd, pixelscale=1/240,
                     focal_length=20)
psf = 0
for wl in np.arange(450e-9, 650e-9, 5e-9):
    w0 = lentil.Wavefront(wavelength=wl)
    w1 = w0 * pupil
    w2 = lentil.propagate_dft(w1, shape=(128,128), pixelscale=5e-6, oversample=3)
    psf += w2.intensity

fig, ax = plt.subplots()
ax.imshow(psf, norm='log')
ax.axis('off')
fig.savefig('broadband_diffraction.png', dpi=150, bbox_inches='tight')
