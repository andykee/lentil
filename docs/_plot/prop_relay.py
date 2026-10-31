import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
import numpy as np
import lentil

# f/10 entrance pupil with an off-axis obscuration so the orientation of the
# relayed pupil is visible
amp = lentil.circle((256, 256), 120) - lentil.circle((256, 256), 20, shift=(40, 60))
amp = lentil.normalize_power(np.clip(amp, 0, 1))
pupil = lentil.Pupil(amplitude=amp, pixelscale=1/240, focal_length=10)

# circular field stop at focus with a 4 lambda F/# (26 um) radius
stop = lentil.Image(amplitude=lentil.circle((256, 256), 26), pixelscale=1e-6)

w = lentil.Wavefront(650e-9) * pupil
w = lentil.propagate_dft(w, to=stop)    # sample focus on the stop's grid

# relay with a 5 m collimator -> 0.5 m diameter output pupil
wa = lentil.propagate_dft(w, pixelscale=1/480, shape=256, focal_length=5)
wb = lentil.propagate_dft(w * stop, pixelscale=1/480, shape=256, focal_length=5)

fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(3.5, 3.5))
ax[0,0].imshow(np.abs(amp), cmap='gray', origin='lower',
             extent=np.array([-128, 128, -128, 128])/240)
ax[0,0].set_title('Entrance pupil (D = 1 m)', fontsize=7)
ax[0,0].set_xlabel('[m]')

ax[0,1].imshow(w.intensity, cmap='turbo', norm=AsinhNorm(1e-4), origin='lower',
             extent=[-128, 128, -128, 128])
ax[0,1].add_patch(plt.Circle((0, 0), 26, fill=False, ec='w', lw=0.6, ls='--'))
ax[0,1].set_title('Image plane + field stop', fontsize=7)
ax[0,1].set_xlabel('[µm]')

for a, wf, title in ((ax[1,0], wa, 'Relayed pupil (no stop)'),
                     (ax[1,1], wb, 'Relayed pupil (with stop)')):
    a.imshow(wf.intensity, cmap='gray', origin='lower',
             extent=np.array([-128, 128, -128, 128])/480)
    a.set_title(title, fontsize=7)
    a.set_xlabel('[m]')

#for a in ax:
#    a.tick_params(labelsize=6)
fig.tight_layout()
