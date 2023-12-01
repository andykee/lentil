import matplotlib.pyplot as plt
import numpy as np
import lentil

amp = lentil.circle((256, 256), 120)
x_tilt = 8e-6 * lentil.zernike(amp, 3)  # +x tilt
y_tilt = 8e-6 * lentil.zernike(amp, 2)  # +y tilt

py = lentil.Pupil(focal_length=10, pixelscale=1 / 240, amplitude=amp, opd=y_tilt)
wy = lentil.Wavefront(650e-9)
wy *= py
wy = lentil.propagate_dft(wy, pixelscale=5e-6, shape=200, oversample=5)

px = lentil.Pupil(focal_length=10, pixelscale=1 /240, amplitude=amp, opd=x_tilt)
wx = lentil.Wavefront(650e-9)
wx *= px
wx = lentil.propagate_dft(wx, pixelscale=5e-6, shape=200, oversample=5)

fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2, ncols=2, figsize=(3.5, 3.5))

ax1.imshow(x_tilt, origin='lower')
ax1.set_title('Pupil plane ($+R_x$)')
ax1.set_xticks(np.linspace(0, 256, 5), labels=np.linspace(-0.5, 0.5, 5))
ax1.set_yticks(np.linspace(0, 256, 5), labels=np.linspace(-0.5, 0.5, 5))
ax1.set_xlabel('[m]')

ax2.imshow(y_tilt, origin='lower')
ax2.set_title('Pupil plane ($+R_y$)')
ax2.set_xticks(np.linspace(0, 256, 5), labels=np.linspace(-0.5, 0.5, 5))
ax2.set_yticks(np.linspace(0, 256, 5), labels=np.linspace(-0.5, 0.5, 5))
ax2.set_xlabel('[m]')

ax3.imshow(wx.intensity ** 0.2, origin='lower')
ax3.set_title('Image plane ($+R_x$)')
ax3.set_xticks(np.linspace(0, 200 * 5, 5), labels=np.linspace(-1, 1, 5))
ax3.set_yticks(np.linspace(0, 200 * 5, 5), labels=np.linspace(-1, 1, 5))
ax3.set_xlabel('[mm]')

ax4.imshow(wy.intensity ** 0.2, origin='lower')
ax4.set_title('Image plane ($+R_y$)')
ax4.set_xticks(np.linspace(0, 200 * 5, 5), labels=np.linspace(-1, 1, 5))
ax4.set_yticks(np.linspace(0, 200 * 5, 5), labels=np.linspace(-1, 1, 5))
ax4.set_xlabel('[mm]')

fig.tight_layout()
