import pyfits as pf
import healpy as hp
from pylab import *
import numpy as np
import matplotlib.pyplot as plt
import math as m
# GENERATES TWO PICTURES: CATALAOG MASKED  AND   RANDOM POINTS MASKED



CAT = pf.open('/home/lucas.bertoncello/files/catalog_sva.fit') # READING .fits TO GET ARRAYS RA,DEC

degtorad = np.pi/180

RA =  CAT[1].data.field('RA').copy()  
phi = RA*degtorad

DEC =  CAT[1].data.field('DEC').copy()
theta = (90-DEC)*degtorad




mask = hp.read_map('/home/lucas.bertoncello/files/mask_SVA.fits') # gives an array as a map healpix form already, has 1s for pixels with galaxies and 0s otherwise
nside = hp.get_nside(mask) # gets the size (resolution) of the given map
npixel = hp.nside2npix(nside)# gets equivalent ammount of pixels depending on size


ID = np.arange(0,npixel,1) #creating empty array os size npixel---arange([start,]stop,[step,]dtype=None)
#initial_map = np.zeros(npixel) #all IDs are given value 0 (UNSEEN) for this array

ID_ = np.hstack((array(ID-0.5),array(max(ID)+0.5)))
 
ID_cat = hp.ang2pix(nside, theta, phi) # gives each angular coordinate a pixel ID (NOT A COMPLETE MAP, ONLY THE POINTS IN 'PIXEL COORDINATES')

cat_map, bin_edges = np.histogram(ID_cat, bins = ID_)# we need a bin storing galaxies for each point there is a galaxy

#cat_map = ID_cat.filled()

hp.mollview(cat_map, title='Catalog_sva', norm ='hist')
plt.title('catalog sva (without mask')
plt.savefig('cat_map.png')

#MASKING-----


#mask_cat =  np.in1d(ID_cat, ID) may also use np.isin()

ID_cat_mask = ID_cat[where(mask)[0]]

#catm_map = hp.SphericalProjAxes.hist(ID_cat, bins = ID)

cat_masked_map,bin_edges = np.histogram(ID_cat_mask,bins = ID_) # bins =  since





# now we have pixel coordinates where there are galaxies (ID_cat)    and   we have a whole map defining our area of interest (mask, 0s for UNSEEN and 1s for used points)
# therefore, wherever (mask==1 && ID_cat==mask) = point in masked_cat

# trying to use loop   while (i <= len(ID_cat))



#MAP_masked = hp.ma(MAP)

#MAP_masked.mask = np.logical_not(mask)



hp.mollview(cat_masked_map, title='Catalog_masked_sva', norm='hist')
plt.title('Catalog_sva Masked')
plt.savefig('cat_ma.png')






# CREATING RANDOM CATALOG WITH nsim POINTS-------------------------------------------------------------------

nsim = 10*len(RA) # len(RA) = 14842110

# generating points between the max and min of the given catalog and getting respective angular coordinates for each

ramin_rad = np.amin(phi)
ramax_rad = np.amax(phi)
decmin_rad = np.amin(theta)
decmax_rad = np.amax(theta)

rasim_rad = (ramax_rad-ramin_rad)*np.array(np.random.uniform(0,1,nsim))+ramin_rad;
decsim_rad = np.arcsin(   np.array(np.random.uniform(0,1,nsim))*(np.sin(decmax_rad)-np.sin(decmin_rad))+np.sin(decmin_rad)    );
cos_of_dec = np.cos((m.pi/2 - decsim_rad))

plt.hexbin( rasim_rad, cos_of_dec,)
plt.colorbar()
plt.xlabel(r'rad_random')
plt.ylabel(r'cos dec random')
plt.title('Random points at cat_sva area')

plt.savefig('rnd_hexbin.png')

'''mtheta, mphi = hp.pix2ang(nside, mask)

plt.hexbin( mphi, mtheta,)
plt.colorbar()
plt.xlabel(r'mask phi')
plt.ylabel(r'mask theta')
plt.title('Mask of cat_sva area')

plt.savefig('m_hexbin.png')'''

ID_rnd = hp.ang2pix(nside, rasim_rad, cos_of_dec)

rnd_map, xedges, yedges = np.histogram(ID_rnd, bins = ID_) # counts how many galaxies there are for each pixel

rnd_map = array(ma_map, dtype = float64) #changing data type (healpy takes in only floats)

#now we have an array 'rnd_map' with floats that gives: how many galaxies there are in each npixels where there are any

hp.mollview(rnd_map, title = 'Random_masked')
plt.title('Random Points Masked')



'''ID_masked = hp.ma(MAP)

MAP_masked.mask = np.logical_not(mask)

plt.savefig('rand_ma.png')
'''
