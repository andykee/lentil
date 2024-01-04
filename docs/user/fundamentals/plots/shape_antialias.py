import matplotlib.pyplot as plt
import lentil

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(4.5, 4.5))
ax[0].imshow(lentil.circle(shape=(64,64), radius=128, shift=(85,-85), antialias=True), 
             cmap='gray')
ax[0].set_title('antialias = True')

ax[1].imshow(lentil.circle(shape=(64,64), radius=128, shift=(85, -85), antialias=False), 
             cmap='gray')
ax[1].set_title('antialias = False')

fig.tight_layout()