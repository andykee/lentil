import numpy as np
import matplotlib.pyplot as plt
import lentil

amp = lentil.circle((256, 256), 120)
x_tilt = 8e-6 * lentil.zernike(amp, 3) # +x tilt
y_tilt = 8e-6 * lentil.zernike(amp, 2) # +y tilt

py = lentil.Pupil(focal_length=10, pixelscale=1/240, amplitude=amp, phase=y_tilt)
wy = lentil.Wavefront(650e-9)
wy *= py
wy = wy.propagate_image(pixelscale=5e-6, npix=256, oversample=5)

px = lentil.Pupil(focal_length=10, pixelscale=1/240, amplitude=amp, phase=x_tilt)
wx = lentil.Wavefront(650e-9)
wx *= px
wx = wx.propagate_image(pixelscale=5e-6, npix=256, oversample=5)

plt.subplot(2, 2, 1)
plt.imshow(x_tilt, origin='lower')
plt.title('Pupil plane [$+R_x$]')
plt.xticks(np.linspace(0, 256, 5), labels=np.linspace(-1, 1, 5), fontsize=8)
plt.yticks(np.linspace(0, 256, 5), labels=np.linspace(-1, 1, 5), fontsize=8)
#plt.grid()

plt.subplot(2, 2, 2)
plt.imshow(y_tilt, origin='lower')
plt.title('Pupil plane [$+R_y$]')
plt.xticks(np.linspace(0, 256, 5), labels=np.linspace(-1, 1, 5), fontsize=8)
plt.yticks(np.linspace(0, 256, 5), labels=np.linspace(-1, 1, 5), fontsize=8)
#plt.grid()

plt.subplot(2, 2, 3)
plt.imshow(wx.intensity**0.2, origin='lower')
plt.title('Image plane [$+R_x$]')
plt.xticks(np.linspace(0, 256*5, 5), labels=np.linspace(-1, 1, 5), fontsize=8)
plt.yticks(np.linspace(0, 256*5, 5), labels=np.linspace(-1, 1, 5), fontsize=8)
plt.subplot(2, 2, 4)
plt.imshow(wy.intensity**0.2, origin='lower')
plt.title('Image plane [$+R_y$]')
plt.xticks(np.linspace(0, 256*5, 5), labels=np.linspace(-1, 1, 5), fontsize=8)
plt.yticks(np.linspace(0, 256*5, 5), labels=np.linspace(-1, 1, 5), fontsize=8)

plt.tight_layout()
plt.savefig('../../_static/img/tilt_coords.png', transparent=False, bbox_inches='tight', dpi=150)