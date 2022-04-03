
import numpy as np

from astropy.io import fits

import matplotlib.pyplot as plt





hdu_list = fits.open('/home/luke/Desktop/IC/dados/3491.fits', memmap=True)
hdu_list.info()

print(hdu_list[1].columns)

hdu_list.close()




CAT = fits.open('/home/luke/Desktop/IC/dados/3491.fits') # estou dando a path total do arquivo, isso significa que vcs nao precisam copiar este arquivo para o seu home e salvar memoria 
RA =  CAT[1].data.field('ra_gal').copy()  #lendo a coluna RA e salvando num array
DEC =  CAT[1].data.field('dec_gal').copy()  #DEC
ZTRUE =  CAT[1].data.field('z_desdm_mc').copy()  #redshift

#print len(RA),len(DEC),len(ZTRUE) # imprimindo o comprimento de cada array
np.savetxt('/home/luke/Desktop/IC/dados/MICE_RA_DEC_ZPHOT.cat',  np.array([RA,DEC,ZTRUE]).transpose(),  fmt='%e',delimiter='\t')
