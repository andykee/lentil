import lentil
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(nrows=1, ncols=6, figsize=(8,3))



opd_cmap = matplotlib.colormaps.get_cmap('RdBu_r')
opd_cmap.set_bad(color='#222222')

a1 = lentil.circle((256,256), 120)
p1 = lentil.Pupil(amplitude=a1, focal_length=10, pixelscale=1/240)
w1 = lentil.Wavefront(650e-9)
w1 = w1 * p1
f1 = lentil.propagate_dft(w1, pixelscale=5e-6, shape=(32,32), oversample=5)
psf1 = f1.intensity

ax1.imshow(a1, cmap='gray')
ax1.axis('off')

ax2.imshow(psf1, cmap='inferno', norm='log', vmin=1e-3)
ax2.axis('off')




segmask = lentil.hex_segments(rings=2, seg_radius=64, seg_gap=2, flatten=False)
amp = np.squeeze(np.sum(segmask, axis=0, keepdims=True))
opd = 0
for seg in segmask:
    coeff = np.random.uniform(-1, 1, 2)*3e-6
    opd += lentil.zernike_compose(seg, [0, *coeff])

p2 = lentil.Pupil(amplitude=amp, opd=opd, focal_length=10, pixelscale=1/240)
w2 = lentil.Wavefront(650e-9)
w2 = w2 * p2
f2 = lentil.propagate_dft(w2, pixelscale=5e-6, shape=(128, 128), oversample=2)
psf2 = f2.intensity

opd[np.where(opd == 0)] = np.nan
ax3.imshow(opd, cmap=opd_cmap)
ax3.axis('off')

ax4.imshow(psf2, cmap='inferno', norm='log', vmin=5e-3)
ax4.axis('off')




outer_diam = 2.4
central_obsc = .33
spider = 0.0264
mount_diam = 0.13
mount_dist = 0.8921

shape = (1024, 1024) # 256
pixelscale = 0.01/4 # .01

npix_outer = outer_diam/pixelscale
npix_inner = (outer_diam * central_obsc)/pixelscale
npix_spider = spider/pixelscale
npix_mount = mount_diam/pixelscale
npix_mount_dist = npix_outer * mount_dist / 2

# primary mirror
hubble_outer = lentil.circle(shape, radius=npix_outer/2, antialias=False)
hubble_inner = lentil.circle(shape, radius=npix_inner/2)
hubble = hubble_outer - hubble_inner

# secondary spiders
for angle in (0, 90, 180, 270): # (45, 135, 225, 315)
    hubble *= lentil.spider(shape, width=npix_spider, angle=angle)

# primary mirror mounting pads
for angle in (0, 120, 240): # (75, 195, 315)
    mount_shift = (npix_mount_dist * -np.sin(np.deg2rad(angle)),
                    npix_mount_dist * np.cos(np.deg2rad(angle)))
    hubble *= 1 - lentil.circle(shape, npix_mount/2, shift=mount_shift)

hubble_opd = lentil.zernike(hubble, index=11)*700e-9

p3 = lentil.Pupil(amplitude=hubble, opd=hubble_opd, focal_length=67.92, pixelscale=2.4/1024)
psf3 = 0
for wl in np.arange(500e-9, 600e-9, 10e-9):
    w3 = lentil.Wavefront(wl)
    w3 = w3 * p3
    f3 = lentil.propagate_dft(w3, pixelscale=5e-6, shape=(512, 512), oversample=2)
    psf3 += f3.intensity
psf3 /= np.max(psf3)

hubble_plot = hubble_opd
hubble_plot[np.where(hubble_plot == 0)] = np.nan
ax5.imshow(hubble_plot, cmap=opd_cmap)
ax5.axis('off')

ax6.imshow(psf3, cmap='inferno', norm='log', vmin=5e-4)
ax6.axis('off')
