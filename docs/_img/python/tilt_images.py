import matplotlib.pyplot as plt
import numpy as np
import lentil

import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (3.5, 3.5)

amp = lentil.circle((256, 256), 120)
x_tilt = 8e-6 * lentil.zernike(amp, 3)  # +x tilt
y_tilt = 8e-6 * lentil.zernike(amp, 2)  # +y tilt

py = lentil.Pupil(focal_length=10, pixelscale=1 / 240, amplitude=amp, phase=y_tilt)
wy = lentil.Wavefront(650e-9)
wy *= py
wy = lentil.propagate_dft(wy, pixelscale=5e-6, shape=200, oversample=5)

px = lentil.Pupil(focal_length=10, pixelscale=1 / 240, amplitude=amp, phase=x_tilt)
wx = lentil.Wavefront(650e-9)
wx *= px
wx = lentil.propagate_dft(wx, pixelscale=5e-6, shape=200, oversample=5)

plt.subplot(2, 2, 1)
plt.imshow(x_tilt, origin='lower')
plt.title('Pupil plane ($+R_x$)')
plt.xticks(np.linspace(0, 256, 5), labels=np.linspace(-0.5, 0.5, 5))
plt.yticks(np.linspace(0, 256, 5), labels=np.linspace(-0.5, 0.5, 5))
plt.xlabel('[m]')

plt.subplot(2, 2, 2)
plt.imshow(y_tilt, origin='lower')
plt.title('Pupil plane ($+R_y$)')
plt.xticks(np.linspace(0, 256, 5), labels=np.linspace(-0.5, 0.5, 5))
plt.yticks(np.linspace(0, 256, 5), labels=np.linspace(-0.5, 0.5, 5))
plt.xlabel('[m]')

plt.subplot(2, 2, 3)
plt.imshow(wx.intensity ** 0.2, origin='lower')
plt.title('Image plane ($+R_x$)')
plt.xticks(np.linspace(0, 200 * 5, 5), labels=np.linspace(-1, 1, 5))
plt.yticks(np.linspace(0, 200 * 5, 5), labels=np.linspace(-1, 1, 5))
plt.xlabel('[mm]')

plt.subplot(2, 2, 4)
plt.imshow(wy.intensity ** 0.2, origin='lower')
plt.title('Image plane ($+R_y$)')
plt.xticks(np.linspace(0, 200 * 5, 5), labels=np.linspace(-1, 1, 5))
plt.yticks(np.linspace(0, 200 * 5, 5), labels=np.linspace(-1, 1, 5))
plt.xlabel('[mm]')

plt.tight_layout()
