# 95073220 galaxias

import numpy as np

npontos = input
zmin = 0.8
zmax = 1.
Z = np.random.uniform(zmin,zmax,npontos)
np.savetxt('randZ', np.transpose(Z))