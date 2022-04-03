#Este codigo gera um mapa de galaxias onde cada pixel na esfera celeste contem informacao sobre o numero de galaxias


import healpy as hp
from pylab import *

RA,DEC = loadtxt('catalog.dat',unpack=True) #lendo o catalago


plot(RA,DEC,'.')
show()

nside = 512  #resolucao do mapa. Este numero nao e arbitrario, tome cuidado para nao escolher um valor errado
npixel = hp.nside2npix(nside) #numero de pixel na esfera com resolucao nside
ID = arange(0,npixel,1) # cada pixel na esfera tem uma identificacao ID

map = zeros(npixel) # mapa inicial onde todos os pixels tem valor zero


phi = deg2rad(RA)
theta = deg2rad(90.-DEC)

ID_cat = hp.ang2pix(nside, theta, phi) #Indentificacao dos pixels que tem galaxias


#criando uma mascara para preencher os valores dos pixels que possuem galaxias
mask_cat =  np.in1d(ID_cat, ID)
ID_cat_mask = ID_cat[where(mask_cat)[0]]


ma_map,bin = histogram(ID_cat_mask,bins=ID) # contando quantas galaxias tem em cada pixel que possuem galaxias


ma_map = array(ma_map, dtype=float64)  #healpy aceita somente numeros floats
map[ID]=ma_map
map[map==0]=hp.UNSEEN # UNSEEN e uma funcaoo do healpy que faz os pixels com 0 galaxias ficarem invisiveis


#plot
map_ = hp.reorder(map,r2n=True)
hp.mollview(map_,title='galaxy map',nest=True)

show()


#  Density galaxy map

# intput :  MASK ID in ring format, nside from Mask and RA,DEC from catalog
# output : map with number os galaxies in each healpix number
#  Flavia Sobreira
#  24/02
 
from pylab import *
from healpy import *
import pyfits as pf
 
 
#Taking ID mask acoording to healpix resolution
print "taking mask ID creating the mask"
 
nside = 4096
ID = loadtxt('mask_id.in', unpack=True)
 
#creating map
print 'lendo o catalogo'
file = '../catalog/clean_mask_absMag_cuts/ALL_catalog_volumeLimited.fits'
hdulist = pf.open(file)
MICE = hdulist[1].data
RA = MICE.field('RA')
DEC = MICE.field('DEC')
 
phi = deg2rad(RA)
theta = deg2rad(90.-DEC)
 
#encontrando o ID pixels onde temos galaxias
ID_cat = ang2pix(nside, theta, phi) #ring
 
print "creating the map"
#vamos fazer um histograma nos ID das galaxias usando como bins os ID da mascara
# Notemos que todo ID e um numero inteiro, entao definimos a largura do nosso bin
#como  I_mask -0.5 e I_mask + 0.5
 
ID_ = np.hstack((array(ID-0.5),array(max(ID)+0.5)))  #aqui estou definindo os bins de ID de acordo com a mascara.
 
numberGal,bin = histogram(ID_cat,bins=sort(ID_))
numberGal = array(numberGal, dtype=float64)  #healpy aceita somente numeros floats
 
 
#salvando mapa no formato healpix
npix = nside2npix(nside)
ID_map = arange(0,npix,1)
map_gal = UNSEEN*ones(npix)
 
map_gal[ID] = numberGal
write_map("NumberGal.fits", map_gal)
 
mollview(map_gal)
show()
#savefig('numberGal_map.pdf'

















#Este programa le as coordenadas de 100.000 galaxias, plota em um plano cartersiado e depois no mapa da esfera celeste e mostra o histograma do numero de galaxias por pixel

import healpy as hp
from pylab import *
import matplotlib.pyplot as plt

RA, DEC = loadtxt('/home/carolina/Downloads/catalogo_galaxia.dat', unpack = True)
plot(RA, DEC, '.')
show()

phi = deg2rad(RA)
theta = deg2rad(90.-DEC)

nside = 512
npixel = hp.nside2npix(nside)

ID = arange(0, npixel, 1)

map = zeros(npixel)

ID_cat = hp.ang2pix(nside, theta, phi)

hp.mollview(ID_cat, title = "Mapa de galaxias da esfera celeste")
show()

plt.hist(ID_cat, 'auto')
plt.title("Numero de galaxias por pixel")
plt.show()



i=0
ID2 = list(range(len(ID_cat) + 2 ))
#cria uma nova array, com 2 elementos a mais que o catalogo
for i in (range(0, len(ID2)))
	ID2[i] = 0

#zerando todos os elementos da array
j=0

for i in range(1, len(ID2) - 1)
	ID2[i] = ID_cat[j]
	j = j + 1

#colocando os elementos do catalogo na nova array, mas deixando vazio o primeiro e ultimo elemento
ID2 = np.array(ID2, dtype=float) #todos os elementos sao float