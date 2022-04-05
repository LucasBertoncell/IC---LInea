import numpy as np
import matplotlib.pyplot as plt



##----ESTIMATING VARIANCE OF ESTIMATORS BASED ON PAPER BY LANDY AND SZALAY----

#note that this is not standard variance, but the variance related to Poisson error

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

DD = []
DR = []
RD = []
RR = []
x = []
xi_of_x = []

file_CUTEoutput = '/home/luke/Desktop/IC/dados/acf_shell3_noPM.dat'#'/home/luke/Desktop/IC/dados/acf_shell3_noPM.dat'

with open(file_CUTEoutput) as CUTE: 

    for line in CUTE:
        column = line.split()
        if not line.startswith('#'): #skipping column labels

            a = float(column[0])
            b = float(column[1])
            c = float(column[2])
            d = float(column[3])
            e = float(column[4])
            f = float(column[5])

            x.append(a)
            xi_of_x.append(b)
            DD.append(c)
            DR.append(d)
            RD.append(e) # not used
            RR.append(f)


n = 123786 #data points
nr = 1000000 #random points at random.dat


theta = x

#---------------------definitions to get variance


np = sum(RR)/len(RR) # average number of unique pair counts (the average is done over multiple random data sets, I WILL BE ONLY COUNTING IT FOR 1 SET, NOT AN AVERAGE)


'''def nt(theta): # average number of unique triplet counts (the average is done over multiple random data sets, I WILL BE ONLY COUNTING IT FOR 1 SET, NOT AN AVERAGE)

	return '''
nt = np /12



Gp = np/(nr*(nr-1)/2) #geometric form factor, prob of finding any two points separeted by distance theta+- dtheta/2


Gt = nt/(nr*(nr-1)*(nr-2)/2)# prob of finding  two points separeted by distance theta+- dtheta/2 given a point


p = 2/(nr*(nr-1)*Gp) #just a simplification

def d(theta):

	wmean = sum(xi_of_x)
	return (1+xi_of_x[x.index(theta)]/(1+wmean)) #just a simplification, mean d, for each theta separation there will be a diferent d


t = (1/n)*(Gt/((Gp)**2)-1) #just a simplification




#----------------using these function to find variances of 4 estimators

PH_var = []
DP_var = []
HAM_var = []
LZ_var = []

i=0
while i < len(theta):

	PH_var.append((d(theta[i])**2)*(4*p + t))
	DP_var.append((d(theta[i])**2)*(p + t))
	HAM_var.append((d(theta[i])**2)*p)
	LZ_var.append((d(theta[i])**2)*p)

	i+=1


plt.title('Estimators Variances')
plt.plot(x, PH_var)
plt.plot(x, DP_var)
plt.plot(x, HAM_var)
plt.plot(x, LZ_var)
plt.legend(['Peebles & Hauser', 'Davis & Peebles', 'Hamilton', 'Landy-Szalay'], loc='upper left')
plt.xlabel('theta')
plt.ylabel('variance')

plt.show()


errorPH = []
errorDP = []
errorHAM = []
errorLZ = []


i=0
while i < len(theta):

	errorPH.append(m.sqrt(PH_var[i])*10000)
	errorDP.append(m.sqrt(DP_var[i]))
	errorHAM.append(m.sqrt(HAM_var[i]))
	errorLZ.append(m.sqrt(LZ_var[i]))

	i+=1

plt.errorbar(x, xi_of_x, yerr = errorPH, fmt='o')
plt.show()