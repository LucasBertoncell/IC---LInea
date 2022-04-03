import numpy as np
import matplotlib.pyplot as plt
from math import *
from scipy.integrate import quad
from scipy.interpolate import interp1d
import sympy as sp
import healpy as hp
'''
-------------------------------------------------------Random Catalog for LSS-----------------------------------------------------------------

TASK
	Since the Dark Energy Survey (DES) will observe at different depths, the measured density of
galaxies in the catalog cannot be translated directly into the mean density of galaxies. In
general, we will observe more galaxies in regions where the survey is deeper and less
galaxies where the survey is shallower. In order to generate the random distribution (useful
for the calculation of the Angular Correlation Function), following the true density field, we
must take into account this information.
	In order to generate non-uniform random catalogs according to the DES depths we have
applied the following methodology: 

1- Create a uniform random catalog according to the galaxy catalog footprint. Use the
Mangle mask to take information about magnitude limit (magL) for each object in the
uniform random catalog and each galaxy in the galaxy catalog. 

2 - Use the Mangle mask to take information about the area, magL of each polygon into
the footprint.
(Instead of Mangle, Healpix can be used)

3 - Create a magL bin (For example, 50 bins from magL 20 to 26).

4 - Count the number of galaxies in each magL and calculate the area.

5- Build the density distribution in each magL bin as the number of galaxies over the
area in the given bin and generate the density(m) function, i.e., the density as a function of
magL.

6- Create the Probability Function according to:

Prob[i] = integral(max_i,min_i)of_density[m] / integral(max_total,min_total)of_density[m]

where Max_ i and Min_i are the maximum and the minimum values of magL in the bin i and
Min_total and Max_total are the initial and final magL according to the user binning. For
example, according to the step 3, Max_total = 26 and Min_total = 20. Equation 1 a probability
function, by definition the sum of Prob[i] in all magL bins must be 1.

7- Normalize the Probability function such the maximum of Prob[i] has value 1( ).

8- Generate a uniform number u for each object in the uniform random catalog. U must be
in the range [0,1]

9- Give a Probability number for each object in the uniform random catalog (given by its
magL). For instance, if an object has magL = 22.56 (belonging to the magL bin between 22.5
and 22.6) and a Prob = 0.6, then the object will have this same probability

10- Compare the uniform number u generated on the step 8 with the Prob of each object in
the uniform random catalog. If u < Prob[i], keep the object in the random catalog, otherwise
remove it from the catalog.
~ 



INPUT: mask_file.fromcat

OUTPUT: random_cat_nonUniform.out

You create a random catalog that is uniform and (make it better/more similar to your data catalog) by making it non uniform

'''

####################################################################################################################################################


#reading given mask, storing ID and respective magL values

mask_IDs, magL = hp.read_map('/home/lucas.bertoncello/files/mask_SVA.fits MUST FIND PATH') # DO WE GET ALL IDS AND 0 MAGLS FOR WHEN THEY ARE UNSEEN?
nside = hp.get_nside(mask_IDs) # gets the size (resolution) of the given map (array of IDs)
npixel = hp.nside2npix(nside) # gets equivalent ammount of pixels depending on size

mask_IDs_magLs = [mask_IDs,magL] #2D
mask_copy = mask_IDs

#uniform random points within pixels that are not covered (generating random points on masked sky)

nsim = 10*len(magL)
ramin_rad = 0
ramax_rad = pi*2
decmin_rad = 0
decmax_rad = pi*2

rasim_rad = (ramax_rad-ramin_rad)*np.array(np.random.uniform(0,1,nsim))+ramin_rad;
decsim_rad = np.arcsin( np.array(np.random.uniform(0,1,nsim))*(np.sin(decmax_rad)-np.sin(decmin_rad))+np.sin(decmin_rad) );

#masking uniformly random points

ID_rand = hp.ang2pix(nside, decsim_rad, rasim_rad) #angular coord transformed into healpy pixel IDs

k = 0
while k < len(magL):

	if mask_IDs[k] = hp.UNSEEN():


	k = k+1

ID = np.arange(0,npixel,1)
ID_empty = np.hstack((array(ID-0.5),array(max(ID)+0.5))) #empty map (IDs) with npixels



ID_rand_masked = ID_rand[where(mask_IDs)[0]] #covering UNSEEN IDs,
rand_map_masked ,bin_edges = np.histogram(ID_rand_masked, bins = ID_empty) #stacking ponts in each ID they appeared (bins are IDs)

#rand_map_masked  will be storing IDs that are not UNSEEN and will have a number associated with each, being the randomly generated amount of points there


'''
TEST

hp.mollview(ID_rand_masked, title='MaskedRandom', norm ='hist')
plt.show
theta, phi = hp.pix2ang(nside, ID_rand_masked)
plt.hexbin(theta, phi)
plt.colorbar()
plt.title('Masked Uniform Random Catalog (without mask)')
plt.savefig('rand_map_uniform.png')

'''

#MagBins

mag_max = np.amax(magL)
mag_min = np.amin(magL)

intervalo_magL = ( mag_max - mag_min )
size_bin = intervalo_magL / sqrt(len(magL)) #basic method
num_bins = intervalo_magL / size_bin

#Counting n_galxies in each magL bin

binned_magL, b_edges = np.histogram(magL, num_bins) # now we have num_bins intervals of size size_bin with a number associated with each interval 

#Finding area 
pixel_area = hp.nside2pixarea(nside, degrees = True) #in squared degrees, area of 1 pixel

#Density

area_bin = binned_magL * pixel_area #in squared degrees, array with areas of each binned magnitude. Each time 

D = binned_magL / area_bin #density of each magL interval, num_bins many
'''
remembering:

[i]: num_bins many bins starting at [0]
D[i]: desity associated with each bin [i]
binned_magL[i]: number of points in each bin [i]
b_edges[i+1]: value of magL of each edge of the bins [i+1] many
rand_map_masked[i]: number of galaxies in each bin
ID_rand[10*i]: 10 times more objects, all IDs representing random point on that pixel
'''
def density(bin_number): # using interpolationm find function density(magL)
	
	return interp1d(bin_number, D[bin_number], kind='cubic') #for each bin_number there is a range of magL associated already

#Porbability function

def Prob(bin_number)

	return  integrate.quad(density(i), b_edges[i], b_edges[i+1]) / integrate.quad(density(i), mag_min, mag_max) #prob of being in bin [i], already normalized

'''
TEST

if np.sum(Prob) = 1:
	print("OK")

'''

#Random number for each ID not UNSEEN (the piexel is the object)


u = np.random.uniform(0,1, len(mask_IDs))

#Prob for each object (ID)

def bin_finder(magnitude):#finds the index of the bin of given magL
	
	i = 0
	while i < num_bins: #considering num_bins was given starting at 1, not 0
		if magnitude >= b_edges[i] and magnitude <= b_edges[i+1]:
			return i
		i=i+1

probability_of_ID = Prob[bin_finder(magL)] #everything if ordered following the index [i]

#Comparing and covering

i = 0
while i < len(mask_IDs):

	if u[i] < probability_of_ID[i]:

		mask_copy[i] = UNSEEN()

	i = i+1

hp.mollview(mask_copy, title='Mask for Uniform random points', norm ='hist')

ID_rand_nonUniform = ID_rand_masked[where(mask_copy)[0]] #covering UNSEEN IDs,
rand_map_nonUniform ,bin_edges_nonU = np.histogram(ID_rand_nonUniform, bins = ID_empty)

plt.histogram()


#covering with new mask printing new non-uniform cat



