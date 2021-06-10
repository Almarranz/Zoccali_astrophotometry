#!/usr/bin/env python
# coding: utf-8

# In[2]:


from astropy.io import fits
import glob
import random
import statistics
from astropy.io.fits import getheader
import os
from astropy import stats
from astropy.stats import sigma_clipped_stats
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.spatial import distance
from astropy.stats import sigma_clip
from astropy.wcs import WCS
from astropy.utils.data import get_pkg_data_filename
from astropy.table import QTable
import pandas as pd
import csv
from csv import writer
from csv import reader


# In[3]:


band='Ks'
exptime=10
#chip=1
folder='im_jitter_NOgains/'
#folder='im_sky_ESOReflex/'
results='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder+'/results/'
pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/pruebas/'
psf='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder
tmp='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder+'tmp/'
indir= '/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder
jitter='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/07.1_Reduce_aligned/054_'+band+'/dit_'+str(exptime)+'/'+folder
pruebas_GNS=indir+'pruebas'
results_GNS=indir+'results/'


# In[ ]:


#hay que hacer los pesos de las images para comprobar las areas de corte de las listas.
for k in range(1,5):
    hdu0 = fits.open(jitter+'mask_cube_chip'+str(k)+'_canvas.fits')[0]
    hdu2= fits.ImageHDU()
    new_hdul = fits.HDUList([hdu0, hdu2])
    wt,header=fits.getdata(jitter+'mask_cube_chip'+str(k)+'_canvas.fits',1,header=True)
    wt_one=np.zeros(shape=(wt.shape[1],wt.shape[2]))
    for i in range(wt.shape[0]):
        wt_one+=wt[i,:,:]
        new_hdul.writeto(tmp+'wt_chip'+str(k)+'.fits',overwrite=True)
        fits.update(tmp+'wt_chip'+str(k)+'.fits',wt_one,header,1,overwrite=True)



# In[48]:


#Add the offset to 'place' the coordiantes on the cube aligned
for chip in range(1,5):
    n=np.loadtxt(jitter+'xy_off_xy_alig_chip'+str(chip)+'.txt')
    print(len(n))
    for i in range(len(n)):
        #print(n[i])
        xc=250-n[i,0]+n[i,2]
        yc=250-n[i,1]+n[i,3]
        #print(xc,yc,i+1)
        coord_b=np.loadtxt(indir+'WHOLE_stars_im'+str(i+1)+'_chip'+str(chip)+'_.txt')
        coord=np.loadtxt(indir+'WHOLE_stars_im'+str(i+1)+'_chip'+str(chip)+'_.txt')
        coord[:,0]+=xc
        coord[:,1]+=yc
        np.savetxt(tmp+'cube_stars_im'+str(i+1)+'_chip'+str(chip)+'.txt',coord,fmt='%.5e')
    print('Done with chip %s'%(chip))
        


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




