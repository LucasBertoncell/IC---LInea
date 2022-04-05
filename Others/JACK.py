
################################# 

import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import subprocess
import skymapper as skm
#import math as m
#import healpy as hp


#################################

'''

Lucas Bertoncello de Oliveira
Unicamp
l182381@dac.unicamp.br




What this does:

-Jackknifes DATA and RANDOM catalogs into N

-Runs CUTE for each jackarea

-Plots Covariance Matrix


'''


def Jackknife(RA,DEC,ZPHOT,N):


	RA_sub = []
	DEC_sub = []
	Z_sub = [] #necessary for CFs but doesnt involve intervals

	i=0

	while i<len(RA):#running over data catalog

		if RA[i]>=RAintervalMin and RA[i]<=RAintervalMax and DEC[i]>=DECintervalMin and DEC[i]<=DECintervalMax: #contains small area to be removed

			n +=1#n is only used for tests

		else: #saves the values that arent into sub area array

			RA_sub.append(RA[i])
			DEC_sub.append(DEC[i])
			Z_sub.append(ZPHOT[i])

		i+=1

	return RA_sub,DEC_sub,Z_sub



###############################---- JACKKNIFE USING INTERVALS, DIVINDING INTO 64 ----##################################### 

#obs.: this divides a rectangular shape, works with generated catalogs like MICE, that are mostly rectangular and usually divided into around 60 



### INPUTS ###

data_filename = '/home/luke/Desktop/IC/dados/example.cat' # columns RA DEC ZPHOT
randomData_filename = '/home/luke/Desktop/IC/dados/example.cat' #same
ACFini_filename = '/home/luke/Desktop/IC/dados/ACF_ini.ini'

##############


RA = np.loadtxt(data_filename, usecols=(0,)) #degrees
DEC = np.loadtxt(data_filename, usecols=(1,)) #degrees
ZPHOT = np.loadtxt(data_filename, usecols=(2,))


rRA = np.loadtxt(randomData_filename, usecols=(0,)) #degrees
rDEC = np.loadtxt(randomData_filename, usecols=(1,)) #degrees
rZPHOT = np.loadtxt(randomData_filename, usecols=(2,))


N = 64 #number of divisions from Jackknife

#because the catalog is circular, we change to circular coordinates, later we will change back to run CUTE

DECrad = np.deg2rad(DEC)
DEC = np.cos((np.pi/2. - DECrad))

'''TEST random points over 1/8 of the sky 

Nrd = 100

ramin = 0
ramax = 90
decmin = 0
decmax = 90

RA = np.random.uniform(ramin,ramax,Nrd)
DEC = np.random.uniform(decmin,decmax,Nrd)
ZPHOT = np.random.uniform(0.6,0.8,Nrd)

'''

#getting the vertices of rectangular catalog
RAmax = max(RA)
RAmin = min(RA) #because MICE cats gives a uncommon point at -270 that screws the intervals--I deleted that point in my cat (point appears at 1662794)
RAbinSize = (RAmax - RAmin) / 8

DECmax = max(DEC)
DECmin = min(DEC)
DECbinSize = (DECmax - DECmin) / 8



#intervals will run over all 8 ranges for both RA and DEC
RAintervalMax = RAmin
RAintervalMin = RAmin #intervals start at the same point
DECintervalMax = DECmin
DECintervalMin = DECmin #intervals start at the same point



#the loop below does N_points*64 ifs
i=0
r=0
d=0

n_cat=0


RA_sub = []
DEC_sub = []
Z_sub = [] #necessary for CFs but doesnt involve intervals

randRA_sub = []
randDEC_sub = []
randZ_sub = []


n=0
ncut = 0


while r < 8: #64 bins, 8 ranges for RA and 8 ranges for DEC
	
	#updating intervals RA
	RAintervalMin = RAintervalMax
	RAintervalMax += RAbinSize

	'''TEST
	print 'RAints'
	print RAintervalMin, RAintervalMax
	'''

	while d < 8:

		#updating intervals DEC
		DECintervalMin = DECintervalMax #min becomes the past max
		DECintervalMax += DECbinSize #max goes forward
		'''TEST
		print 'DECints'
		print DECintervalMin, DECintervalMax
		'''

		while i<len(RA):#running over data catalog

			if RA[i]>=RAintervalMin and RA[i]<=RAintervalMax and DEC[i]>=DECintervalMin and DEC[i]<=DECintervalMax: #contains small area to be removed

				n +=1#n is only used for tests

			else: #saves the values that arent into sub area array

				RA_sub.append(RA[i])
				DEC_sub.append(DEC[i])
				Z_sub.append(ZPHOT[i])

			i+=1
			
		'''TEST
		print n
		'''
		i=0#reseting counter

		'''
		print len(RA_copy)
		print len(DEC_copy)
		print len(Z_copy)
		'''

		#a weights column is necessary to run CUTE
		num = len(RA_sub)

		weights = np.full(shape=num,fill_value=1,dtype=np.float)


		jackarea = np.array([RA_sub, DEC_sub, Z_sub, weights]).transpose() # THIS IS WHERE DATA IS SAVED


		f = open("jacksubArea.txt", "w+")
		np.savetxt('jacksubArea.txt',  np.array([RA_sub, DEC_sub, Z_sub, weights]).transpose(), fmt='%e',delimiter='\t')
		f.close()


		#now doing the same for RANDOM CAT

		while i<len(rRA):#running over random catalog

			if rRA[i]>=RAintervalMin and rRA[i]<=RAintervalMax and rDEC[i]>=DECintervalMin and rDEC[i]<=DECintervalMax: #contains small area to be removed

				n +=1#n is only used for tests

			else: #saves the values that arent into sub area array

				randRA_sub.append(RA[i])
				randDEC_sub.append(DEC[i])
				randZ_sub.append(ZPHOT[i])

			i+=1

		i=0

		num = len(randRA_sub)

		rweights = np.full(shape=num,fill_value=1,dtype=np.float)

		Rjackarea = np.array([randRA_sub, randDEC_sub, randZ_sub, rweights]).transpose() # THIS IS WHERE RANDOM DATA IS SAVED


		n_cat+=1 #index of jack catalogs


		f = open("randomsubArea.txt", "w+")
		np.savetxt('randomsubArea.txt',  np.array([randRA_sub, randDEC_sub, randZ_sub, weights]).transpose(), fmt='%e',delimiter='\t') #it will keep overwritting
		f.close()

		print 'done{0}'.format(n_cat)
		print n
		print ncut
		ncut=0







		#----OPTION SAVE ALL SUB_CATS: after running over the catalog and appending all points within the intervals, the sub_areas can be saved below ---------------
		'''
		f = open("jackArea{0}.txt".format(n_cat), "w+")
		np.savetxt('jackArea{0}.txt'.format(n_cat),  np.array([RA_cop,DEC_cop,ZPHOT_cop,weights]).transpose(),  fmt='%e',delimiter='\t')
		f.close()


		f = open("jacksubArea.txt", "w+")
		np.savetxt('jacksubArea.txt',  np.array([RA_sub, DEC_sub, Z_sub, weights]).transpose(), fmt='%e',delimiter='\t')
		f.close()


		'''
		#-----------------------------------------------------------------------------------------------------------------------------------------------------------




		#---- Running CUTE  

		#obs.:always check that the number of each line that is being changed is correct
		#jacksubArea.txt and randomsubArea.txt are changed each ittertion, but name remains the same
		print 'starting CUTE run{0}.out\n'.format(n_cat)
		
		with open('/home/luke/Desktop/IC/Codes/CUTE-master/ACF_ini.ini', 'r') as file: #opens ini file to read
	    	ini_file = file.readlines()# read a list of lines into data

	    # changing name of output file on ini file
		ini_file[1] = 'output_filename= output_from_area{0}.out\n'.format(n_cat)

	    with open('/home/luke/Desktop/IC/Codes/CUTE-master/ACF_ini.ini', 'w') as file: #writes over ini file with chnages to lines
   			file.writelines(ini_file)

		subprocess.call(['nohup','condor_do','CUTE','ACF_ini.ini'])# ACTUALLY Running , could use shell=True makes executes command through the shell
		



		#reseting jackareas random and data
		randRA_sub = []
		randDEC_sub = []
		randZ_sub = []	

		RA_sub = []
		DEC_sub = []
		Z_sub = []

		d+=1


	#once the while that runs over DEC runs 8 times it has to be reset in order to run over the next RA column, as well as the counter d
	DECintervalMax = DECmin
	DECintervalMin = DECmin
	d = 0

	r+=1




###################################################################################################################################

#########################---------------COVARIANCE MATRIX (from CUTE outputs) ------------------###################################

#output format from CUTE: first column = separations, second = w(separation), the others are DD,RR,RD,DR
#in the case of 64 jackknife divisions and 30 angles, it will construct 30x30 matrix
# w is used to refer to a calculated ACF and theta to the angle for which it was calculated
#we consider the CUTE input was 'bins = 30', so 30 angle bins
#it can be written in a simpler way, I wrote this way so I could understand each step better






#getting the mean for each angule
w = []
thetas = []
w_mean = [] #it contains a mean for each angle, so 30 means
N = 64

#loop below does 30*N calculations (it could be done faster but it would be harder to read the code)
i=0
k=0
while i<30: #runs over 30 angules

	w_sum = 0#summatory empty

	while k<N: #runs over N subsamples

		w = np.loadtxt('output_from_area{0}.out'.format(k), usecols=(1,)) #opens output file of each 

		w_sum += w[i] #gets i w(theta) of k subsample

		k+=1

	w_mean.append(w_sum / N) 

	k=0
	i+=1




#getting all 30 w_thetas for all 64 subsamples
sample, angle = 64, 30
w_jack = [[0 for x in range(angle)] for y in range(sample)]  #creates 64 sets of 30 (matrix 64 x 30) to store all w_thetas [sample_index][angle_index] 




i=0
while i < sample:

	w_jack[i] = np.loadtxt('output_from_area{0}.out'.format(i), usecols=(1,))

	i+=1

#obs.: indexes go from 0 to 63 and 0 to 29








#now with the w_means[] and w_jack[][] we compute the values of the matrix -- loop does 30*30*64 iterations

angle1, angle2 = 30, 30
C_matrix = [[0 for x in range(angle1)] for y in range(angle2)]

Norm_factor = (N-1)/N

i=0 #angle 1
j=0 #angle 2
k=0 #sub samples

while i<30:

	while j<30:

		summation = 0

		while k<N:

			product = (w_jack[k][i] - w_mean[i])*(w_jack[k][j] - w_mean[j])

			summation += product

			k+=1

		C_matrix[i][j] =(Norm_factor*summation)

		j+=1
		k=0

	i+=1
	j=0






#logs = np.log10(C_matrix)  trying log scale

'''
#plotting inverse
inverse = np.linalg.inv(C_matrix)
plot= plt.matshow(inverse)

plt.xticks(range(len(thetas)), thetas, rotation=90)
plt.yticks(range(len(thetas)), thetas)

plt.xlabel('Separation [degrees]')
plt.ylabel('Separation [degrees]')
#plt.colorbar()
plt.show()
plt.savefig('ĩnverse2.png')
'''

#plotting normal matrix
plot= plt.matshow(C_matrix) #matrix has to be in form matrix , norm=LogNorm(vmin=0, vmax=1)

plt.xticks(range(len(thetas)), thetas, rotation=90)
plt.yticks(range(len(thetas)), thetas)

plt.xlabel('Separation [degrees]')
plt.ylabel('Separation [degrees]')
#plt.colorbar()
plt.show()
plt.savefig('normalCM.png')













#################################---- JACKKNIFE USING HEALPY ----#########################################(only testing)

'''
data_filename = '/home/anderson/CUTE/MICE/shell_0.dat' # columns RA DEC ZPHOT
random_filename= '/home/anderson/CUTE/MICE/random.dat'

RA = np.loadtxt(data_filename, usecols=(0,)) #degrees
DEC = np.loadtxt(data_filename, usecols=(1,)) #degrees
ZPHOT = np.loadtxt(data_filename, usecols=(2,))


phi=np.deg2rad(full_catalogRA)
theta=np.deg2rad(full_catalogDEC)


nside = 512


i=0
list_pixels = []
while i<len(theta):

	pixelID = hp.ang2pix(nside,theta[i],phi[i]) #running through list of RA DEC ZPHOT
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

	f = open("jackArea{0}.txt".format(k), "w+")
	numpy.savetxt('jackArea{0}.txt'.format(k),np.array([area_ra,area_dec]).transpose(),fmt='%e',delimiter='\t')
	f.close()
	#numpy savetxt

	k+=1


'''

