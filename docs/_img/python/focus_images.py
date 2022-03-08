import matplotlib.pyplot as plt
import numpy as np
import lentil

import matplotlib as mpl
mpl.rcParams['figure.figsize'] = (4, 4)

amp = np.zeros((256, 256))
amp += lentil.circlemask((256, 256), 64, shift=(32, 8))
amp -= lentil.circlemask((256, 256), 40, shift=(32, 8))
amp[:, 0:128+8] = 0
amp[32:225, 104-16:128-16] = 1
amp[96:120, 128-16:144-8] = 1
amp[201:225, 128-16:144-8] = 1

focus = lentil.zernike(mask=np.ones((256, 256)), index=4)

pupil_neg = lentil.Pupil(amplitude=amp, pixelscale=1/240, focal_length=10)
pupil_neg.phase = -6e-6 * focus
w_neg = lentil.Wavefront(650e-9)
w_neg *= pupil_neg
w_neg = w_neg.propagate_image(pixelscale=5e-6, npix=200, oversample=2)

pupil_pos = lentil.Pupil(amplitude=amp, pixelscale=1/240, focal_length=10)
pupil_pos.phase = 6e-6 * focus
w_pos = lentil.Wavefront(650e-9)
w_pos *= pupil_pos
w_pos = w_pos.propagate_image(pixelscale=5e-6, npix=200, oversample=2)

fig, axs = plt.subplots(1, 2)

ax1 = axs[0]
ax1.imshow(w_pos.intensity, origin='lower')
ax1.set_title('Image plane (+focus)')
ax1.set_xticks(np.linspace(0, 400, 5), labels=np.linspace(-1, 1, 5))
ax1.set_yticks(np.linspace(0, 400, 5), labels=np.linspace(-1, 1, 5))
ax1.set_xlabel('[mm]')

ax2 = axs[1]
ax2.imshow(w_neg.intensity, origin='lower')
ax2.set_title('Image plane (-focus)')
ax2.set_xticks(np.linspace(0, 400, 5), labels=np.linspace(-1, 1, 5))
ax2.set_yticks(np.linspace(0, 400, 5), labels=np.linspace(-1, 1, 5))
ax2.set_xlabel('[mm]')

plt.tight_layout()
