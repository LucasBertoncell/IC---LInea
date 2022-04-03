import numpy as np
import matplotlib.pyplot as plt
import math as m

nsim = 1000000

ramin = 0
ramax = 90
decmin = 0
decmax = 90


#transformando para deg to rad
degtorad = m.pi/180.
ramin_rad = ramin*degtorad
ramax_rad = ramax*degtorad
decmin_rad =  decmin*degtorad
decmax_rad = decmax*degtorad

#creating the random points
rasim_rad = (ramax_rad-ramin_rad)*np.array(np.random.uniform(0,1,nsim))+ramin_rad;
decsim_rad = np.arcsin(   np.array(np.random.uniform(0,1,nsim))*(np.sin(decmax_rad)-np.sin(decmin_rad))+np.sin(decmin_rad)    );

random = np.zeros(shape = (len(rasim_rad), 4))
random[:,0] = rasim_rad*180/np.pi
random[:,1] = decsim_rad*180/np.pi
random[:,2] = np.array(np.random.uniform(0,1,nsim)) + 0.2
random[:,3] = np.ones(len(rasim_rad))

np.savetxt('random.dat', random)

#Making the plot
#plot in a plane
'''
plt.subplot(1,2,1)
plt.hexbin(rasim_rad,decsim_rad)
plt.colorbar()
plt.xlabel(r'RA[rad]')
plt.ylabel('DEC[rad]')
plt.title('Plane')
#plot in a sphere
plt.subplot(1,2,2)
cosdec = np.cos((m.pi/2. - decsim_rad))
plt.hexbin(rasim_rad,cosdec)
plt.colorbar()
plt.xlabel(r'RA[rad]')
plt.ylabel(r'COS(DEC)')
plt.title('Sphere')
plt.subplots_adjust(wspace=0.5, hspace=0.5)

plt.savefig('footPrint.png')
'''