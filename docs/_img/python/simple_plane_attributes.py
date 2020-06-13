import numpy as np
import matplotlib.pyplot as plt
import monocle as mo

circlemask = mo.util.circlemask((256, 256), 128)
plt.imshow(circlemask)
plt.savefig('../../_static/img/circle_mask.png', transparent=True, bbox_inches='tight', dpi=150)

amplitude = mo.util.circle((256, 256), 128)
plt.imshow(amplitude)
plt.savefig('../../_static/img/circle_amplitude.png', transparent=True, bbox_inches='tight', dpi=150)

z4 = mo.zernike.zernike(circlemask, 4)
plt.imshow(z4)
plt.savefig('../../_static/img/circle_focus.png', transparent=True, bbox_inches='tight', dpi=150)
