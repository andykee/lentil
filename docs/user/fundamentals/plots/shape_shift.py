import matplotlib.pyplot as plt
import lentil

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(4.5, 4.5))
ax[0].imshow(lentil.circle(shape=(256,256), radius=32), cmap='gray')
ax[0].plot(128,128, '+', markersize=6, color='#bf616a')
ax[0].set_title('shift = (0, 0)')

ax[1].imshow(lentil.circle(shape=(256,256), radius=32, shift=(-80,50)), cmap='gray')
ax[1].plot(128,128, '+', markersize=6, color='#bf616a')
ax[1].arrow(128, 128, 50, -80, color='#bf616a', linewidth=0.5, head_width=8, length_includes_head=True)
ax[1].set_title('shift = (-80, 50)')

fig.tight_layout()