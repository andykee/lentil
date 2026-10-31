import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
import numpy as np
import lentil

amp = lentil.circle((256, 256), 120)
pupil = lentil.Pupil(amplitude=amp, pixelscale=1/240, focal_length=10)

# only compute the region covered by an off-axis 40 x 60 px subarray
mask = np.zeros((1280, 1280))
#mask[150:190, 160:220] = 1
mask[750:950, 800:1100] = 1


w = lentil.Wavefront(650e-9) * pupil
w_full = lentil.propagate_dft(w, pixelscale=2e-6, shape=256, oversample=5)
w_mask = lentil.propagate_dft(w, pixelscale=2e-6, shape=256, mask=mask,
                              oversample=5)

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(4, 2.2))
for a, wf, title in ((ax[0], w_full, 'mask=None'), (ax[1], w_mask, 'mask=mask')):
    a.imshow(wf.intensity, cmap='turbo', norm=AsinhNorm(1e-4), origin='upper')
    a.set_title(title, fontsize=7, fontfamily='monospace')
    a.axis('off')
ax[0].add_patch(plt.Rectangle((800, 750), 60*5, 40*5, fill=False, ec='w', lw=0.6, ls='--'))
fig.tight_layout()
