import numpy as np
import matplotlib.pyplot as plt
from math import *
from scipy.integrate import quad


#################### EX 5


filename = 'bao_z_theta_erro(theta).dat'
z, theta, error_theta = np.loadtxt(filename, unpack=True, usecols=(0,1,2))

# PARAMETERS

mass_dens = 0.3
rad_dens = 0.000084 
vel_sound = 0.02106781

# FUNCTIONS



def E(z):

	return sqrt((1 - mass_dens) + mass_dens*((1 + z)**3) + rad_dens*((1+z)**4) )

def E2(z): #for integral

	return 1/E(z)

def com_dist(z):

	integral, error = (quad(E2, 0, z))

	return integral


def theta_teo(Ez):

	return 

def dist_comoving_squared(theta_teo_list, theta_observed_list):

	j = 0
	while len(theta) > j:


		j+=1

	np.sum(summ)


	










######################## EX 4


#  EXTRACTING DATA

filename = 'Legacy_z_(m-M)_erro(m-M).dat'
z, m_M, error_m_M = np.loadtxt(filename, unpack=True, usecols=(0,1,2))



# PARAMETERS


mass_dens = 0.3
rad_dens = 0.000084 

# FUNCTIONS


def mod_of_dist(z,h):

	return (5*np.log10(dist_lum(z,h)) + 25)




def dist_lum(z,h):

	(3000/h)*(1+z)*com_dist(z)



def values_teo(z,h):


	list_teo = []
	j = 0

	while j < len(z):

		list_teo.append(dist_lum(z[j],h))
		j+=1

	return list_teo


h = 0.5     # initial value

while h <= 1:


	Xsquared_min = lum_dist_comoving_squared(list_teo , m_M)


	lum_dist_comoving_squared()
	if Xsquared_min < 
	