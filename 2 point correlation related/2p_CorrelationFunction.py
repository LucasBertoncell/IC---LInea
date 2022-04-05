'''
Lucas Bertoncello de Oliveira
'''

import numpy as np
import matplotlib.pyplot as plt
from math import *
from scipy.integrate import quad
import sympy as sp

#------------------------------- READING A CATALOG AND COUNTING DISTANCES OF PAIRS OF POINTS AT GIVEN RADIUS------------------------------------
# THE IDEA
# Calculate the Dist(array), distances between all points. Chose your R (radius) and check how many points are there with Dist == R (this will be the value DD)
# Do the same for the Randomly generated points (10x more points than used catalog) and get RR
# ((len(randomCat)/len(Cat))**2)*(DD/RR) - 1 = E (2 point correlation function value for chosen R)
# Get values of E for multiple Rs (bins)
# Then plot a curve: Correlation Function vs R(Mpc) (expecting BAO peak)


'''
INPUT
observed_cat.in (file containing: column1 as RA, column2 as DEC, column3 as Z)
random_cat.in (file containing random values (10x number of points observed_cat.in): column1 as RAr, column2 as DECr, column3 as Zr)
radius_min  = 10 (Mpc)
radius_max = 190 (Mpc)
n_cfs_bins = 30 (would give you for instance 180/30=6, 30 bins with 6 Mpc in between each)
resolution_npix = 2048 or 4096 ...

OUTPUT
chosen_name.out (containing column1 as r, column2 as xi(r))

'''




# ---------RANDOM CATALOG RR (CAN BE OMITTED IF GIVEN A RANDOM CATALOG)-==============================================================

npontos = 200
'''
MASKING --

ramin = 10 #issues if using 0?
ramax = 12airra
decmin = 10
decmax = 20
Area = (ramax - ramin)*(decmax - decmin)

degtorad = pi/180

ramin_rad = ramin*degtorad
ramax_rad = ramax*degtorad

decmin_rad = decmin*degtorad
decmax_rad = decmax*degtorad

'''
def X(z):
	return (299792/70)*1/(sqrt(0.3*(z + 1)**3 + 0.7))

def COStheta(dec,dec_, ra,ra_):
	return sin(dec)*sin(dec_) + cos(dec)*cos(dec_)*cos(ra - ra_)

def r(coseno,z1,z2):
	return sqrt((X(z1))**2 + (X(z2))**2 - 2*X(z1)*X(z2)*coseno)


RA = (ramax_rad-ramin_rad)*np.array(np.random.uniform(0,1,npontos))+ramin_rad;
DEC = np.arcsin(   np.array(np.random.uniform(0,1,npontos))*(np.sin(decmax_rad)-np.sin(decmin_rad))+np.sin(decmin_rad)    );

zmin = 0.6
zmax = 0.8
7
Z = np.random.uniform(zmin,zmax,npontos)

H = 70
m = 0.3
A = 0.7
c = 299792

z = sp.Symbol('z')
RRDist = []

i = 0
while (i < npontos):
	distancia , erro = (quad(X, 0, Z[i]))# function 'quad' returns both numerical result and error associated 

	RRDist.append(distancia)

	i=i+1
#===============================================================================================================================================


# ------ REAL CATALOG DD

RAd,DECd,Zd = np.loadtxt('catalog.dat',unpack=True) #reading to get RA DEC Z

Dist = []
Contado = np.ones(npontos) #1s pra ja contados
z = sp.Symbol('z')

i = 0
j = 0

while (i < npontos):

	#whiles vao diminuindo, o ultimo while nao tem que contar (todos pares ja foram contados) ---> recurcao?
	while (j < npontos):

		if Contado[j] == 1 and i != j: #ponto nao foi contado && nao calculo dist do ponto com ele mesmo

			coseno = COStheta(DECd[i], DECd[j], RAd[i], RAd[j])

			distanciaD = r(coseno, Zd[i], Zd[j])

			Dist.append(distanciaD)

		j=j+1

	Contado[i] = 1 #esse ponto ja teve sua distancia calculada em relacao a todos outros pontos
	i=i+1
	j=0

# Now we have RRDist (list of distances between all random galaxies) and Dist (list of distances between all observed galaxies)
# Note the unit of these distances depends on what CF you are dealing with.
# In this case, we are given RA, DEC and Z and get r in Mpc

#--------------------SEPARATING OUR DISTANCES INTO GROUPS
'''

We have all distances between all pairs of points
We separate these in "groups" based on the bins and ranges:
	from the example, we have 30 bins with 6Mpc in between, staring at 10 until 190,
	we will search for distances at [10,16), [16,22), [22,28), ...., [185, 190]
	the number of pairs distant [10,16) will be DD of our first bin, and so on for the other bins, getting 30 DDs
	We do the same for the random_cat and get 30 RRs 
	With RRs and DDs we calculate our xi(r) following some estimator (for now, (N*(DD[r]/RR[r]) - 1) = xi(r), where N = (len(randomCat)/len(Cat))**2))

'''


