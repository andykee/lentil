import matplotlib.pyplot as plt
import numpy as np
import lentil

amp = lentil.hex_segments(rings=2, seg_radius=32, seg_gap=2, flatten=True)
fig, ax = plt.subplots()
ax.imshow(amp, cmap='gray')
ax.axis('off')
fig.savefig('segmented_diffraction.png', dpi=150, bbox_inches='tight')
