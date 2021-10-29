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

xx, yy = np.meshgrid(np.linspace(0, 1, mask.shape[0]), np.linspace(0, 1, mask.shape[1]))

X = xx
Y = yy

fig = plt.figure()

ax1 = fig.add_subplot(131)
ax1.imshow(mask)
ax1.set_title('mask', fontname='monospace', fontsize=9)
ax1.axis('off')

ax2 = fig.add_subplot(132, projection='3d', proj_type='ortho')
ax2.view_init(elev=25, azim=235)

offset = np.linspace(-0.3,1.2, 6)
alpha = 0.97

xs = [0, 0, 1, 1, 0]
ys = [0, 1, 1, 0, 0]

xss = [1, 1, 0, 0, 0.2]
yss = [0.28, 0, 0, 1, 1]

ax2.contourf(X, Y, segmask[3], 10, zdir='z', offset=offset[5], alpha=alpha, zorder=20)
ax2.plot(xs, ys, offset[5], linewidth=1, color='k', zorder=100)

ax2.contourf(X, Y, segmask[2], 10, zdir='z', offset=offset[4], alpha=alpha, zorder=20)
ax2.plot(xss, yss, offset[4], linewidth=1, color='k', zorder=99)

ax2.contourf(X, Y, segmask[1], 10, zdir='z', offset=offset[3], alpha=alpha, zorder=20)
ax2.plot(xss, yss, offset[3], linewidth=1, color='k', zorder=100)

ax2.contourf(X, Y, segmask[0], 10, zdir='z', offset=offset[2], alpha=alpha, zorder=20)
ax2.plot(xss, yss, offset[2], linewidth=1, color='k', zorder=100)

ax2.contourf(X, Y, segmask[4], 10, zdir='z', offset=offset[1], alpha=alpha, zorder=20)
ax2.plot(xss, yss, offset[1], linewidth=1, color='k', zorder=100)

ax2.contourf(X, Y, segmask[5], 10, zdir='z', offset=offset[0], alpha=alpha, zorder=20)
ax2.plot(xss, yss, offset[0], linewidth=1, color='k', zorder=100)

ax2.axis('off')
ax2.set_title('segmented_mask', fontname='monospace', fontsize=9)

ax3 = fig.add_subplot(133)
ax3.imshow(segmask[0])
ax3.set_title('segmented_mask[0]', fontname='monospace', fontsize=9)
ax3.axis('off')

plt.savefig('../../_static/img/segmask.png', transparent=True, bbox_inches='tight', dpi=300)

