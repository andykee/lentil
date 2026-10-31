import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
import lentil

amp = lentil.circle((256, 256), 120)
pupil = lentil.Pupil(amplitude=amp, pixelscale=1/240, focal_length=10)
detector = lentil.Image(pixelscale=2e-6, mask=lentil.rectangle((64, 64), 64, 64))

w = lentil.Wavefront(650e-9) * pupil
w10 = lentil.propagate_dft(w, to=detector, oversample=5)
w20 = lentil.propagate_dft(w, to=detector, focal_length=20, oversample=5)

fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(4, 2.2))
ext = [-64, 64, -64, 64]
for ax, wf, title in ((ax1, w10, 'focal_length=10'),
                      (ax2, w20, 'focal_length=20')):
    ax.imshow(wf.intensity, cmap='turbo', norm=AsinhNorm(), origin='lower',
              extent=ext)
    ax.set_title(title, fontsize=7)
    ax.set_xlabel('[µm]')
ax1.set_ylabel('[µm]')
fig.tight_layout()