import numpy as np
import matplotlib.pyplot as plt
from math import *
from scipy.integrate import quad
import sympy as sp

npontos = 200

ramin = 10 #issues if using 0?
ramax = 12
decmin = 10
decmax = 20
Area = (ramax - ramin)*(decmax - decmin)

degtorad = pi/180

ramin_rad = ramin*degtorad
ramax_rad = ramax*degtorad

decmin_rad = decmin*degtorad
decmax_rad = decmax*degtorad

RA = (ramax_rad-ramin_rad)*np.array(np.random.uniform(0,1,npontos))+ramin_rad;
DEC = np.arcsin(   np.array(np.random.uniform(0,1,npontos))*(np.sin(decmax_rad)-np.sin(decmin_rad))+np.sin(decmin_rad)    );

#RANDOM REDSHIFTS------------------------------------------------

zmin = 0.6
zmax = 0.8

Z = np.random.uniform(zmin,zmax,npontos)

#DISTANCIA COMOVEL-----------------------------------------------

H = 70
m = 0.3 
A = 0.7
c = 299792


Dist_entre_pares = []
Contado = np.ones(npontos) #1s pra ja contados
z = sp.Symbol('z')

def X(z):
	return (299792/70)*1/(sqrt(0.3*(z + 1)**3 + 0.7))

def COStheta(dec,dec_, ra,ra_):
	return sin(dec)*sin(dec_) + cos(dec)*cos(dec_)*cos(ra - ra_)

def r(coseno,z1,z2):
	return sqrt((X(z1))**2 + (X(z2))**2 - 2*X(z1)*X(z2)*coseno)

i = 0
j = 0

while (i < npontos):

	#whiles vao diminuindo, o ultimo while nao tem que contar (todos pares ja foram contados) ---> recurcao?
	while (j < npontos):

		if Contado[j] == 1 and i != j: #ponto nao foi contado && nao calculo dist do ponto com ele mesmo

			coseno = COStheta(DEC[i], DEC[j], RA[i], RA[j])

			distancia = r(coseno, Z[i], Z[j])

			Dist_entre_pares.append(distancia)

		j=j+1

	Contado[i] = 1 #esse ponto ja teve sua distancia calculada em relacao a todos outros pontos
	i=i+1
	j=0

plt.hist(Dist_entre_pares, bins = 15)
plt.xlabel(r'Distancia comovel entre 2 galaxias (Mpc)')
plt.ylabel(r'Ocorrencias')
plt.title('HIstograma: Distancia comovel entre todos pares de galaxias de uma amostra')
plt.show()