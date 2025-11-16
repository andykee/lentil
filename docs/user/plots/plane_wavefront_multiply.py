import lentil
import matplotlib.pyplot as plt
import numpy as np

amp = lentil.circle((256,256), 120)
opd = lentil.zernike(amp, 4)*500e-9

plane = lentil.Plane(amp, opd)

w = lentil.Wavefront(500e-9)
print(w.field)

w2 = w * plane

fig, ax = plt.subplots(nrows=2, ncols=5, width_ratios=[1,0.4,1,0.4,1], figsize=(5,3))

ax[0,0].imshow(np.ones_like(amp), vmin=0, vmax=1)
ax[0,0].tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
ax[0,0].set_xlabel('wavefront.amplitude')

ax[1,0].imshow(np.zeros_like(amp), vmin=0, vmax=1)
ax[1,0].tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
ax[1,0].set_xlabel('wavefront.phase')

ax[0,1].annotate('x', xy=(0.45,0.5))
ax[0,1].axis('off')

ax[1,1].annotate('x', xy=(0.45,0.5))
ax[1,1].axis('off')

ax[0,2].imshow(plane.amplitude), ax[0,2].set_title('')
ax[0,2].tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
ax[0,2].set_xlabel('plane.phase')

ax[1,2].imshow(plane.opd)
ax[1,2].tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
ax[1,2].set_xlabel('plane.opd')

ax[0,3].annotate('=', xy=(0.4,0.48))
ax[0,3].axis('off')

ax[1,3].annotate('=', xy=(0.4,0.48))
ax[1,3].axis('off')

ax[0,4].imshow(np.abs(w2.field))
ax[0,4].tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
ax[0,4].set_xlabel('wavefront.amplitude')

ax[1,4].imshow(np.angle(w2.field))
ax[1,4].tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
ax[1,4].set_xlabel('wavefront.phase')