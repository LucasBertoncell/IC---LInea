"""
Lucas Bertoncello de Oliveira
182381


len()
int()
"""

import random 
""" 
"Mersenne Twister pseudorandom number generator as the core generator"
uniform()

"""
import numpy as np 
"""
array()
append()
histogram2d()

"""
import matplotlib.pyplot as plt 
import matplotlib as mpl
"""
scatter()
imshow()
axes()


"""
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import math as m


array_X = [] #Criando arrays vazios para x e y
array_Y = []

#-----------PARAMETROS------------------
n_Pontos = 1000

xmin = 10
xmax = 20
ymin = 10
ymax = 20

ramin = 10
ramax = 20
decmin = 0
decmax = 90

nbin = 10

raio = 5
#---------------------------------------
i = 0
x = 0
y = 0
while (i < 100):
	
	x = x+1

	while (y<=10):
		
		array_X.append(x) # acrescenta cada x aleatorio no final do array_X
		array_Y.append(y)
		y = y+1

	i = i+1


print("\n1 - Scatter retangular\n2 - Histogram2d\n3 - Histogram2d na superficie de uma esfera\n4 - random points no plano vs esfera\n5 - random points Masked\n6 - Sair")
loop = True
while loop:
	option=input("Opcao:") 

	if option == '1': # scatter------------------------------------------------------------------------------------------


i = 0
while(i<n_Pontos): #Gerando coordenadas (x,y) aleatorias n vezes
	x = random.uniform(xmin,xmax) # retorna um float aleatorio dentro do intervalo (min, max)
	y = random.uniform(ymin,ymax)

	"""
	para pontos unicos, nenhum x ou y repetidos
	if x not in array_X:
		array_X.append(x)
	if y not in array_Y:
		array_Y.append(y)
	"""
	
	array_X.append(x) # acrescenta cada x aleatorio no final do array_X
	array_Y.append(y) 
	i = i+1


plt.xlabel(r'$RA$',fontsize = 15)
plt.ylabel(r'$COS(DEC)$',fontsize = 15)

plt.axes().set_aspect('equal', 'box') #'box' ajusta o tamanho das axis para um quadrad,'equal' faz as escalas de x e y iguais

plt.scatter(array_X,array_Y, marker='.', edgecolor='none', facecolor='k', alpha=0.7)

plt.show()





	elif option == '2':# histogram2d------------------------------------------------------------------------------------



i = 0
while(i<n_Pontos): #Gerando coordenadas (x,y) aleatorias n vezes

	x = random.uniform(xmin,xmax) # retorna um float aleatorio dentro do intervalo (min, max)
	y = random.uniform(ymin,ymax)
	
	array_X.append(x) # acrescenta cada x aleatorio no final do array_X
	array_Y.append(y) 
	i = i+1


hist,xedges,yedges = np.histogram2d(array_X,array_Y,bins=nbin) #"Compute the bi-dimensional histogram of two data samples"
															   # as duas dimencoes, (X,Y) tem bins = nbin
															   # edges sao a beira de cada histograma

plt.xlabel(r'$X$',fontsize = 15)
plt.ylabel(r'$Y$',fontsize = 15)

extent = (np.amin(array_X), np.amax(array_X), np.amin(array_Y), np.amax(array_Y))# limites

plt.imshow(hist.T, interpolation='bessel', extent=extent, cmap='jet', norm=LogNorm())
# usa um array e retorna uma imagem MxN(dependendo do tipo de array)



plt.colorbar()
plt.show()






	elif option == '3': #gerar pontos numa area sobre a superficie de uma esfera-----------------------------------------
		pi = math.pi

		
i = 0
while(i<n_Pontos): # Gerando coordenadas (x,y) aleatorias n vezes

	phi = random.uniform(0,2*pi)  # gerando angulos, em rad, aleatorios
	
	costheta = random.uniform(-1,1)

	theta = math.acos( costheta )

	#theta = random.uniform(0,2*pi)
	
	x = raio * math.cos(theta) * math.sin(phi)# coordenadas, de esfericas para cartesianas
	y = raio * math.sin(theta) * math.sin(phi)

	# obs.: RAIO CONSTANTE (pontos na supericie da esfera)
	#raio variante : R = raio *(u ** (1 / 3)) , u = random.uniform(0,1)
	array_X.append(x)
	array_Y.append(y)

	i = i+1
		

		# posso tambem gerar ponto no plano cartesiano e traze-los para o esferico:

		i = 0
		while(i<n_Pontos): #Gerando coordenadas (x,y) aleatorias n vezes

			xmin = -5
			xmax = 5
			ymin = -5
			ymax = 5
			x = random.uniform(xmin,xmax) # retorna um float aleatorio dentro do intervalo (min, max)
			y = random.uniform(ymin,ymax)

			XsqMaisYsq = x**2 + y**2
			zsqard = (raio**2 - x**2 - y**2)

			if zsqard < 0:
				zsqard = zsqard*(-1)

			z = math.sqrt(zsqard)
			
			theta = math.atan2(z, math.sqrt(XsqMaisYsq)) #conversao para esfericas
			phi = math.atan2(y,x)

			X = raio * math.cos(theta) * math.sin(phi)# coordenadas, de esfericas para cartesianas
			Y = raio * math.sin(theta) * math.sin(phi)

			array_X.append(X)
			array_Y.append(Y) 
			i = i+1
		

			#metodo pronto:
			'''
			def appendSpherical_np(x,y,z):
				ptsnew = np.hstack((xyz, np.zeros(xyz.shape)))
				xy = xyz[:,0]**2 + xyz[:,1]**2
				ptsnew[:,3] = np.sqrt(xy + xyz[:,2]**2)
				ptsnew[:,4] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
				#ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
				ptsnew[:,5] = np.arctan2(xyz[:,1], xyz[:,0])
				return ptsnew
			z = 5;
			'''


hist,xedges,yedges = np.histogram2d(array_X,array_Y,bins=nbin)


plt.xlabel(r'$X$',fontsize = 10)
plt.ylabel(r'$Y$',fontsize = 10)

#plt.axes().set_aspect('equal', 'box')

#extent = (np.amin(array_X), np.amax(array_X), np.amin(array_Y), np.amax(array_Y))# limites

plt.hexbin(array_X,array_Y,bins="log", gridsize=500, cmap='magma')

plt.colorbar()
plt.show()


#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


	elif option == '4':

nsim = 100000

degtorad = m.pi/180

ramin_rad = ramin*degtorad
ramax_rad = ramax*degtorad

decmin_rad =  decmin*degtorad
decmax_rad = decmax*degtorad

#creating the random points
rasim_rad = (ramax_rad-ramin_rad)*np.array(np.random.uniform(0,1,nsim))+ramin_rad;
decsim_rad = np.arcsin(   np.array(np.random.uniform(0,1,nsim))*(np.sin(decmax_rad)-np.sin(decmin_rad))+np.sin(decmin_rad)    );

#Making the plot
#plot in a plane
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
#plt.subplots_adjust(wspace=0.5, hspace=0.5)

plt.show()

#plt.savefig('footPrint.png')


elif option == '5': # MASKED

nsim = 100000

degtorad = m.pi/180

ramin_rad = ramin*degtorad
ramax_rad = ramax*degtorad

decmin_rad =  decmin*degtorad
decmax_rad = decmax*degtorad

rasim_rad = (ramax_rad-ramin_rad)*np.array(np.random.uniform(0,1,nsim))+ramin_rad;
decsim_rad = np.arcsin(   np.array(np.random.uniform(0,1,nsim))*(np.sin(decmax_rad)-np.sin(decmin_rad))+np.sin(decmin_rad)    );



	elif option == '6':
		loop = False


""" Outras opcoes:

plt.scatter(x,y, marker='o',s = 50, edgecolor='none', facecolor='k', alpha=0.2)

fig, ax = plt.subplots()
ax = fig.add_subplot(111)
ax.set_aspect(1) # quadrado
--usado para arrumar a ration das axis pra um quadrado

x = numpy.random.poisson(lam=1.0, size=100)
y = numpy.random.poisson(lam=1.0, size=100)
--funcao poisson


x = np.array(array_X)
y = np.array(array_Y)
--tornando x e y em arrays

color:
	magma
	inferno
	cm.hot

plt.hist(array_X, 100)

		fig = plt.figure()
		ax = fig.add_subplot(111, title='Interpolated',aspect='equal', xlim=xedges[[0, -1]], ylim=yedges[[0, -1]])
		im = mpl.image.NonUniformImage(ax, interpolation='bilinear')

		xcenters = (xedges[:-1] + xedges[1:]) / 2
		ycenters = (yedges[:-1] + yedges[1:]) / 2

		im.set_data(xcenters, ycenters, hist.T)
		ax.images.append(im)
	

UNSEEN - nao mostrar pixels sem ponto (cinza)

plot de densidade --> healpy, 


ssh -i ~/Desktop/SSHkey.pub -p 53222 -XC lucas.bertoncello@devel2.linea.gov.br
scp -P 53222 mycode.tar.gz lucas.bertoncello@devel2.linea.gov.br
scp -P 53222    lucas.bertoncello@devel2.linea.gov.br:/mnt/scratch/users/fsobreira/devel/MICE/catalog/MICEv2_BLUE_mag23_04to06.fits .

??  //  ||  \\  
"""