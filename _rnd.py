import numpy as np
import matplotlib.pyplot as plt
import math as m

nsim = 1000

ramin = 10
ramax = 20
decmin = 0
decmax = 90


degtorad = m.pi/180

ramin_rad = ramin*degtorad
ramax_rad = ramax*degtorad

decmin_rad =  decmin*degtorad
decmax_rad = decmax*degtorad

rasim_rad = (ramax_rad-ramin_rad)*np.array(np.random.uniform(0,1,nsim))+ramin_rad;
decsim_rad = np.arcsin(   np.array(np.random.uniform(0,1,nsim))*(np.sin(decmax_rad)-np.sin(decmin_rad))+np.sin(decmin_rad)    );


cosdec = np.cos((m.pi/2. - decsim_rad))
plt.hexbin(rasim_rad,cosdec)

plt.colorbar()
plt.xlabel(r'RA(rad)')
plt.ylabel(r'COS(DEC)')
plt.title('Sphere')

plt.show()
