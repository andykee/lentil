import matplotlib.pyplot as plt
import numpy as np
import lentil

seg = lentil.hexagon((256, 256), 128)
segmask = [np.zeros((768, 768)) for s in np.arange(6)]
segmask[0][0:256, 256:512] += seg
segmask[1][128:384, 37:293] += seg
segmask[2][383:639, 37:293] += seg
segmask[3][512:768, 256:512] += seg
segmask[4][128:384, 475:731] += seg
segmask[5][383:639, 475:731] += seg

segmask = np.asarray(segmask)
mask = np.sum(segmask, axis=0)


fig, ax = plt.subplots(nrows=1, ncols=4, width_ratios=[1,0.3,1,1], figsize=(4,1.5))

plt.figtext(0.23,0.85,"mask.ndim = 2", va="center", ha="center", size=8, fontname='monospace')
plt.figtext(0.68,0.85,"mask.ndim = 3", va="center", ha="center", size=8, fontname='monospace')

ax[0].imshow(mask)
ax[0].tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

ax[1].axis('off')

ax[2].imshow(np.ones_like(mask), cmap='gray', vmin=0, vmax=1)
ax[2].tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
ax[2].set_frame_on(False)
ax[2].set_xlabel('mask', fontname='monospace', size=7)

for n, loc in enumerate(np.flip(np.linspace(0, 0.25, 6))):
    ax_in = ax[2].inset_axes([loc, loc, 0.75, 0.75])
    ax_in.imshow(segmask[5-n])
    ax_in.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

ax[3].imshow(segmask[0])
ax[3].tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
ax[3].set_xlabel('mask[0]', fontname='monospace', size=7)