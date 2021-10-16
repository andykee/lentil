import numpy as np

import field

f1 = field.Field(data=np.ones((5,4)), pixelscale=1, offset=(-2,-2))
f2 = field.Field(data=np.ones((3,3)), pixelscale=1, offset=(0,-1))

print(f1.extent)
print(f2.extent)
