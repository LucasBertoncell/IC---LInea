
import numpy as np
import matplotlib.pyplot as plt

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

file_CUTEoutput = '/home/luke/Desktop/IC/dados/acf_shell3_noPM.dat'

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





## RODANDO OS ESTIMáDORES

#using wc -l
Nrd = 20000000 
N = 123786
# 1426966   968327    263726   123786


##NORMALIZED PAIR COUNTS-----------------------------------------------------

DDnorm = []
i=0
while i < len(DD):
	DDnorm.append( DD[i]/( (N*(N - 1)/2)))
	i+=1
RRnorm = []
i=0
while i < len(RR):
	RRnorm.append( RR[i]/(  (Nrd*(Nrd - 1)/2) ) )
	i+=1
DRnorm = []
i=0
while i < len(DD):
	DRnorm.append(DR[i]/(Nrd*N))
	i+=1


#NORMALIZED ESTIMATORS--------------------------------------------------------------

PH = []
i=0
while i < len(DD):


	PH.append((DDnorm[i]/RRnorm[i]) - 1.)
	i+=1


DP = []
i=0
while i < len(DD):

	DP.append((DDnorm[i]/DRnorm[i]) - 1.)
	i+=1


HAM = []
i=0
while i < len(DD):

	HAM.append((DDnorm[i]*RRnorm[i])/(DRnorm[i]**2) - 1)
	i+=1


LZ = []
i=0
while i < len(DD):

	LZ.append( (DDnorm[i] - 2*DRnorm[i] + RRnorm[i]) / RRnorm[i])
	i+=1

##--------------------------------------------------------------------------------------using the non normalized pair counts



norm_fac = Nrd/N # what is the diference from term Nrd-1, N-1 ...



PH = []
i=0
while i < len(DD):


	PH.append((norm_fac**2)*(DD[i]/RR[i]) - 1.)
	i+=1


DP = []
i=0
while i < len(DD):

	DP.append((norm_fac**2)*(DD[i]/DR[i]) - 1.)
	i+=1


HAM = []
i=0
while i < len(DD):

	HAM.append((norm_fac**2)*((DD[i]*RR[i])/(DR[i]**2))-1)
	i+=1


LZ = []
LZ2 = []
LZ3 = []
fac1 = (Nrd*(Nrd - 1. )) / (N*(N - 1.))
fac2 = ((Nrd - 1.) / N)

fac3 = (Nrd*Nrd)/(N*N)
fac4 = Nrd / N

i=0

while i < len(DD):

	LZ.append(1 + (norm_fac**2)*(DD[i]/RR[i]) - 2*norm_fac*DR[i]/RR[i])
	LZ2.append( (fac1*DD[i] - fac2*DR[i] + RR[i])/RR[i])

	LZ3.append( (fac3*DD[i] - fac4*DR[i] + RR[i])/RR[i])

	i+=1

##----------------------------------------------------------------------------------VARIANCE


import math as m

mean = sum(PH)/len(PH)

i=0
dif_from_mean=[]

while i<len(PH):
	dif_from_mean.append((mean - PH[i])**2)
	i+=1

Var = m.sqrt(sum(dif_from_mean)/(len(dif_from_mean) -1))









## -----------------------------getting percentual difference



plt.plot(x, PH,'r--', label = "Peebles&Hauser")
plt.xlabel("$\Theta$ [deg]")
plt.ylabel("$w(\Theta)$")
plt.title("Peebles & Hauser estimator - MICE catalog 0")
plt.show()

plt.plot(x, LZsquared)
plt.xlabel("$\Theta$ [deg]")
plt.ylabel("$w(\Theta) $\Theta$^2 $")
plt.title("Landy-Szalay - MICE catalog 0")
plt.show()



# plot with various axes scales
plt.figure(1)

plt.subplot(231)
plt.plot(x, PHsquared)
plt.title('Peebles and Hauser')

plt.subplot(232)
plt.plot(x, DPsquared )
plt.title('Davis and Peebles')

plt.subplot(233)
plt.plot(x, HAMsquared )
plt.title('Hamilton')

plt.subplot(234)
plt.plot(x, LZ )
plt.title('Landy-Szalay')

plt.subplot(235)
plt.plot(x, LZsquared )
plt.title('Landy-Szalay * theta^2 ')

plt.subplot(236)
plt.plot(x, diff )
plt.title('Percentual Difference: LZ and HAM')

plt.show()






LZsquared = []
l=0
while l < len(DD):

	LZsquared.append(LZ[l]*x[l]*x[l])
	l+=1

PHsquared = []
l=0
while l < len(DD):

	PHsquared.append(PH[l]*x[l]*x[l])
	l+=1

HAMsquared = []
l=0
while l < len(DD):

	HAMsquared.append(HAM[l]*x[l]*x[l])
	l+=1

DPsquared = []
l=0
while l < len(DD):

	DPsquared.append(DP[l]*x[l]*x[l])
	l+=1




diffPH_LS = []
i = 0
while i < len(HAM):

	diffPH_LS.append(((LZ[i]  - PH[i])/PH[i])*100) #array
	i += 1


diffHAM_LS = []
i = 0
while i < len(HAM):

	diffHAM_LS.append(((LZ[i]  - HAM[i])/HAM[i])*100) #array
	i += 1


plt.figure(1, figsize= (15,15))

plt.subplot(421)
plt.title('Estimators Deviation')
plt.plot(x, PH)
plt.plot(x, DP)
plt.plot(x, HAM)
plt.plot(x, LZ)
plt.legend(['Peebles & Hauser', 'Davis & Peebles', 'Hamilton', 'Landy-Szalay'], loc='upper left')


plt.subplot(422)
plt.title('Estimators Deviation (times theta^2 for magnification) ')
plt.plot(x, PHsquared)
plt.plot(x, DPsquared)
plt.plot(x, HAMsquared)
plt.plot(x, LZsquared)
plt.legend(['Peebles & Hauser', 'Davis & Peebles', 'Hamilton', 'Landy-Szalay'], loc='upper left')


plt.subplot(423)
plt.title('Percentual Difference - Peebles & Hauser and Landy-Szalay ')
plt.plot(x, diffPH_LS, 'r')
plt.xlabel('theta')
plt.ylabel('percentage [%]')


plt.subplot(424)
plt.title('Percentual Difference - Hamilton and Landy-Szalay ')
plt.plot(x, diffHAM_LS, 'r')
plt.xlabel('theta')
plt.ylabel('percentage [%]')

plt.subplots_adjust( wspace=0.5, hspace=0.8)
plt.show()





#SAVING THE FILE 
#np.savetxt('estimators_shell0.dat', zip(PH, xi_of_x, DD, DR, RD, RR), fmt=" %5.9f %5.9f %5.9f %5.9f %5.9f %5.9f")



plt.plot(x, PH)
plt.plot(x, LZ)

plt.title('ACF: Peebles & Hauser e Landy-Szalay')
plt.legend(['Peebles & Hauser', 'Landy-Szalay'], loc='upper left')
plt.xlabel('theta^2 [rad]')
plt.ylabel('ACF')
plt.show()


plt.plot(x, diffPH_LS, 'r')
plt.title('Percentual Difference: Peebles & Hauser e Landy-Szalay')
plt.xlabel('theta [rad]')
plt.ylabel('porcentagem [%]')
plt.show()