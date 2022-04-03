
################################# 
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import math as m



## JACKNIFE RESAMPLING ##

'''

- generate random catalog
- decide on: number of pieces to divide & number of cfs calculated on each piece
- run through list of (ra, dec) getting values that are within range of each piece
- calculate CF for each sub area of (ra, dec) 
- do this for all pieces, keeping a sum of each CF bin
- by the end of running through lists and calculating cfs you have a list of summed up cfs
- divide the summed up values by the number of pieces to get the MEAN

'''

#--- inputs

Nrd = 1000

#1/8 of the sky
ramin = 0
ramax = 90
decmin = 0
decmax = 90

nside = 1024


raRAND_deg = np.random.uniform(ramin,ramax,Nrd)
decRAND_deg = np.random.uniform(decmin,decmax,Nrd)

phi = np.deg2rad(raRAND_deg)
theta = np.deg2rad(decRAND_deg)


i=0
list_pixels = [] #will contain pixel Identity Numbers for the points from the catalog using the healpy format

while i<len(theta): #points are grouped together into pixels of given nside resolution, to check all pixels that actually contain a point

	pixelID = hp.ang2pix(nside,theta[i],phi[i])  #points in form of pixels
	is_this_pix_not_new =  pixelID in list_pixels #true if piexelID is in list of runned pixels

	if is_this_pix_not_new == False: #if this pix is new, it appends it to the list 

			list_pixels.append(pixelID)
		
	i+=1



k=0
while k<len(list_pixels): #running through list of Ids

	area_ra = []#these will contain ra and dec of each Jackknife area
	area_dec = []

	i=0
	while i<len(theta):

		pixelID = hp.ang2pix(nside,theta[i],phi[i]) #getting the ids of each pixel again

		if pixelID == list_pixels[k]: #if that Id is the same as one of the areas

			area_ra.append(phi[i])
			area_dec.append(theta[i])

		i+=1

	#f = open("jackArea{0}.txt".format(i), "w+")
	#np.savetxt('jackArea{0}.txt'.format(k),np.array([area_ra,area_dec]).transpose(),fmt='%e',delimiter='\t')
	#f.close()
	#numpy savetxt

	plt.scatter(area_ra,area_dec)

	k+=1


#COV matrix

wmean = sum(output)/len(output)

i = 0
j = 0
norm = (N - 1 / N)
Dif =[]

while  i < m.sqrt(N):

	while j < sqrt(N):

		Dif.append((ouput[i] -  wmean)*(output[j] - wmean))

		j+=1

	i+=1




COV = np.matrix('1 2; 3 4')


##





plt.xlabel("RA[deg]")
plt.ylabel("DEC[deg]")
plt.title("Random catalog divided into Jackknife areas NSIDE = 1024")
plt.show()

#---------------------------------------------------------------------------------




nside = 256
data_filename= "/home/anderson/CUTE/MICE/shell_0.dat" # columns RA DEC ZPHOT

RA = np.loadtxt(data_filename, usecols=(0,))
DEC = np.loadtxt(data_filename, usecols=(1,))

phi= np.deg2rad(RA)
theta=  np.deg2rad(DEC)

i=0
list_pixels = []
while i<len(theta):

	pixelID = hp.ang2pix(nside,theta[i],phi[i], lonlat=True)
	is_this_pix_not_new =  pixelID in list_pixels

	if is_this_pix_not_new == False:
			list_pixels.append(pixelID)
	i+=1


k=0
while k<len(list_pixels):

	area_ra = []
	area_dec = []
	i=0

	while i<len(theta):

		pixelID = hp.ang2pix(nside,theta[i],phi[i], lonlat=True)

		if pixelID == list_pixels[k]:

			area_ra.append(phi[i])
			area_dec.append(theta[i])
			

		i+=1

	#f = open("jackArea{0}.txt".format(i), "w+")
	numpy.savetxt('jackArea{0}.txt'.format(i),np.array([area_ra,area_dec]).transpose(),fmt='%e',delimiter='\t')
	f.close()
	

	k+=1


#---------------------------------------------------------------------------------

#cycol = cycle('bgrcmky')
# colors: one of {'b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'}; 8 types
#----generating random points-------------------------------------------------------

degtorad = m.pi/180.
ramin_rad = ramin*degtorad
ramax_rad = ramax*degtorad
decmin_rad = decmin*degtorad
decmax_rad = decmax*degtorad


#rasim_rad = np.random.uniform(ramin_rad,ramax_rad,Nrd)
#decsim_rad = np.random.uniform(decmin_rad,decmax_rad,Nrd)

rasim_rad = (ramax_rad-ramin_rad)*np.array(np.random.uniform(0,1,Nrd))+ramin_rad;
decsim_rad = np.arcsin(   np.array(np.random.uniform(0,1,Nrd))*(np.sin(decmax_rad)-np.sin(decmin_rad))+np.sin(decmin_rad)    );

cosdec = np.cos((m.pi/2. - decsim_rad))



phi = rasim_rad
theta = decsim_rad

mapa = hp.ang2pix(nside, theta, phi, lonlat=True)
hp.mollview(mapa, title="nside=512")


#procuro por quais coordenadas então no mesmo pixel (fazendo uso dos pixels para separar)
i=0
list_pixels = []

while i<len(theta):

	pixelID = hp.ang2pix(nside,theta[i],phi[i], lonlat=True)
	is_this_pix_not_new =  pixelID in list_pixels

	if is_this_pix_not_new == False:

			list_pixels.append(pixelID)
		
	i+=1

	areas.append(pixelID)

	if hp.ang2pix(nside,theta[i],phi[i]) == 


	i+=1



areas = [][]

while k < len(areas):

	plt.scatter(copy_ra, copy_dec)


plt.xlabel("RA[rad]")
plt.ylabel("DEC[rad]")
plt.title("Random catalog divided into Jackknife areas - NSIDE = 1024")
plt.show()


#----getting vertices 

lenght_ra = ramax_rad/n_pieces
lenght_dec = decmax_rad/n_pieces

vertice1 = {(decmax_rad - lenght_dec), ramin_rad}
vertice2 = {decmax_rad, ramin_rad}
vertice3 = {decmax_rad, lenght_ra}
vertice4 = {(decmax_rad - lenght_dec), lenght_ra}


#----getting the map NESTED

map_rnd = hp.ang2pix(nside, theta, phi)

hp.mollview(map_rnd, norm='hist')
plt.show()


#----find coordinates that should be UNSEEN NESTED

copy_ra = rasim_rad #this will store the full catalog with the marked coordinates
copy_dec = decsim_rad
sub_ra = []
sub_dec = []
ids2cover = []

k = 0
i = 0
where_to_remove = 0

dinamic_Nrd = Nrd
while k < n_pieces:

	while i < Nrd:

		if (decmax_rad - lenght_dec) < decsim_rad[i]  and rasim_rad[i] < lenght_ra:


			sub_ra.append(rasim_rad[i])
			sub_dec.append(decsim_rad[i])

			copy_ra = np.delete(copy_ra, where_to_remove)#remove
			copy_dec = np.delete(copy_dec, where_to_remove)

			
			where_to_remove -=1
		

		where_to_remove+=1
		i+=1
	#---- covering each sub_section and calculating CF n_pieces times


	k+=1

n_bins = 30 #30 cfs  for each sub_area
CFsum = [0]*n_bins #certifying all [] are empty





plt.scatter(copy_ra, copy_dec)
plt.scatter(sub_ra,sub_dec)
plt.xlabel("RA")
plt.ylabel("DEC")
plt.grid(color = 'k', linestyle = '--', linewidth = 0.1)
plt.title("Random catalog divided into 9 Jackknife areas (area highlighted)")
plt.show()


plt.scatter(rasim_rad, decsim_rad)
plt.xlabel("RA [rad]")
plt.ylabel("DEC [rad]")
plt.title("Random points plot RA x DEC")
plt.show()

plt.scatter(rasim_rad, cosdec)
plt.xlabel("RA [rad]")
plt.ylabel("cos(pi/2 - DEC) [rad]")
plt.title("Random points plot RA x cos(pi/2 - DEC)")
plt.show()

mask_ra = rasim_rad < lenght_ra 
mask_dec = decsim_rad > (decmax_rad - lenght_dec)

plt.scatter(rasim_rad[mask_ra], cosdec[mask_dec])
plt.xlabel("RA [rad]")
plt.ylabel("cos(pi/2 - DEC) [rad]")
plt.title("Random points plot RA x cos(pi/2 - DEC)")
plt.show()



	### save values in a file1
	### run cfs with sub_ra , sub_dec
	### store found values in file2

# CF_values is an array containing the values CF medians over 10 areas

## vertices-----------------------------araquivo--------------------------mapa com cores


###-----------------------------------------------------------------------------------------------------------------------------for all areas

'''

- transform RA DEC --> theta phi [radians]
- ang2pix with chosen nside -->npix
- nside2pixarea(nside[, degrees])Gives pixel area given nside in square radians or square degrees. --> A_tot of catalog
- A_Tot/n_pieces = A_of_pieces
- find point of maximun RA or DEC (pixels?) --> starting point, a border, pixel1
- get_all_neighbours(nside, theta[, phi, ...])	Return the 8 nearest pixels --> staring at pixel1 --> keep spreading
- each pixel counted adds to a sum_of_subArea -- > once subArea = A_of_pieces you stop spreading and calculate CF for this area
problem: will it cover the whole area no matter the shape? the spread could get pluged


resolução 

areas nao precisam ser iguais

'''

##

