import matplotlib.pyplot as plt
import lentil

import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (4.5, 4.5)

amp = lentil.circle((256, 256), 120)
opd = lentil.zernike(amp, 4) * 1e-6
pupil = lentil.Pupil(amplitude=amp, opd=opd, pixelscale=1 / 240, focal_length=10)

w1 = lentil.Wavefront(wavelength=500e-9)
w1 *= pupil
w1 = lentil.propagate_dft(w1, pixelscale=5e-6, shape=128, prop_shape=128, oversample=5)

w2 = lentil.Wavefront(wavelength=500e-9)
w2 *= pupil
w2 = lentil.propagate_dft(w2, pixelscale=5e-6, shape=128, prop_shape=40, oversample=5)


fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(4, 4))

ax1.imshow(w1.intensity, cmap='inferno')
ax1.set_title('prop_shape ok')
ax1.axis('off')

ax2.imshow(w2.intensity, cmap='inferno')
ax2.set_title('prop_shape too small')
ax2.axis('off')
