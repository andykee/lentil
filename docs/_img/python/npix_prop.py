import matplotlib.pyplot as plt
import lentil

import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (4.5, 4.5)

amp = lentil.circle((256, 256), 120)
opd = lentil.zernike(amp, 4) * 1e-6
pupil = lentil.Pupil(amplitude=amp, phase=opd, pixelscale=1 / 240, focal_length=10)

w1 = lentil.Wavefront(wavelength=500e-9)
w1 *= pupil
w1 = w1.propagate_image(pixelscale=5e-6, npix=128, npix_prop=128, oversample=5)

w2 = lentil.Wavefront(wavelength=500e-9)
w2 *= pupil
w2 = w2.propagate_image(pixelscale=5e-6, npix=128, npix_prop=40, oversample=5)

plt.subplot(1, 2, 1)
plt.imshow(w1.intensity, origin='lower')
plt.title('npix_prop ok')
plt.axis('off')

plt.subplot(1, 2, 2)
plt.imshow(w2.intensity, origin='lower')
plt.title('npix_prop too small')
plt.axis('off')
