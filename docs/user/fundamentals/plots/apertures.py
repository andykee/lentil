import lentil
import matplotlib.pyplot as plt

fig, (ax1,ax2,ax3, ax4) = plt.subplots(nrows=1, ncols=4, figsize=(5, 3))

ax1.imshow(lentil.circle((256,256), 100), cmap='gray')
ax1.axis('off')

ax2.imshow(lentil.rectangle((256,256), 100, 150), cmap='gray')
ax2.axis('off')

ax3.imshow(lentil.hexagon((256,256), 110), cmap='gray')
ax3.axis('off')

c4 = lentil.circle((256,256), 100) - lentil.circle((256,256), 25)
c4 *= lentil.spider((256,256), 4)
c4 *= lentil.spider((256,256), 4, angle=90)
c4 *= lentil.spider((256,256), 4, angle=180)
c4 *= lentil.spider((256,256), 4, angle=270)
ax4.imshow(c4, cmap='gray')
ax4.axis('off')
