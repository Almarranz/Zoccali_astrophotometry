#!/usr/bin/env python
# coding: utf-8

# In[ ]:
#generates a file with coordinates and mag and uncertaintes 

from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
from astropy.table import QTable
import json


# In[ ]:


band='H'
exptime=10
#chip=1
folder='im_jitter_NOgains/'
#folder='im_sky_ESOReflex/'
results='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder+'/results/'
pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/pruebas/'
indir = '/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder
psf='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder
tmp='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder+'tmp/'
name='NPL_054'


# In[ ]:


dic_x={}
dic_y={}
dic_mag={}
dic_dmag={}

for chip in range(1,5):
    with open(tmp+'x_commons_chip'+str(chip)+'.json') as fx:
        dic_x['x_chip%s'%(chip)] = [tuple(x) for x in json.load(fx)]
    with open(tmp+'y_commons_chip'+str(chip)+'.json') as fy:
        dic_y['y_chip%s'%(chip)] = [tuple(x) for x in json.load(fy)]
    dic_mag['mag_chip%s'%(chip)]=np.loadtxt(tmp+'mag_media_chip%s.txt'%(chip))   
    dic_dmag['dmag_chip%s'%(chip)]=np.loadtxt(tmp+'error_mag_med_chip%s.txt'%(chip))
for chip in range(1,5):
    x_mean=[np.mean((dic_x['x_chip%s'%(chip)][i])) for i in range(len(dic_x['x_chip%s'%(chip)]))]
    y_mean=[np.mean((dic_y['y_chip%s'%(chip)][i])) for i in range(len(dic_x['x_chip%s'%(chip)]))]
    dx=[0.106*np.std(dic_x['x_chip%s'%(chip)][i])/np.sqrt(len(dic_x['x_chip%s'%(chip)][i]))
        for i in range(len(dic_x['x_chip%s'%(chip)]))]
    dy=[0.106*np.std(dic_y['y_chip%s'%(chip)][i])/np.sqrt(len(dic_x['x_chip%s'%(chip)][i]))
        for i in range(len(dic_x['x_chip%s'%(chip)]))]
    
    
    f = fits.open(tmp+'wt_chip%s.fits'%(chip))
    w = WCS(f[1].header)
    print(w)
    coord_gal=[w.pixel_to_world(x_mean[i], y_mean[i]).galactic for i in range(len(x_mean))]
    coord=[w.pixel_to_world(x_mean[i], y_mean[i]) for i in range(len(x_mean))]
    t_ra = QTable([coord], names=["lines coord"])
    t_gal= QTable([coord_gal], names=["lines coord"])
    df_ra=t_ra.to_pandas()
    ra=df_ra.to_numpy()#ra is (ra,dec)
    df_gal=t_gal.to_pandas()
    gal=df_gal.to_numpy()
    total=np.c_[ra,x_mean,dx,y_mean,dy, dic_mag['mag_chip%s'%(chip)],dic_dmag['dmag_chip%s'%(chip)],gal[:,0],gal[:,1]]
    np.savetxt(tmp+name+'_chip%s.txt'%(chip),total,header='ra,dec,x_mean,dx,y_mean,dy,mag,dmag,l,b',fmt='%.6f')
    print('Done with chip %s'%(chip))

# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




