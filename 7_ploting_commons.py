#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.stats import gaussian_kde


# In[2]:


band='Ks'
exptime=10
#chip=1
folder='im_jitter_NOgains/'
#folder='im_sky_ESOReflex/'
results='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder+'/results/'
pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/pruebas/'
indir = '/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder
psf='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder
tmp='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder+'tmp/'
GNS='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/field_Bulge18/'
name='NPL_054'


# In[3]:


#RA,DRA,DEC,DDEC,raJ,draJ,decJ,ddecJ,raH,draH,decH,ddecH,raKs,draKs,decKs,ddecKs,J,dJ,H,dH,Ks,dKs=np.loadtxt(GNS+'catalogue_F18_c3.txt',unpack=True)
field18=np.loadtxt(GNS+'catalogue_F18_c3.txt')
GNS_gal = SkyCoord(ra=field18[:,0]*u.degree, dec=field18[:,2]*u.degree, frame='icrs').galactic
t = QTable([GNS_gal], names=['coord'])
df=t.to_pandas()
gal_GNS=df.to_numpy()
field18_gal1=np.c_[field18,gal_GNS[:,0],gal_GNS[:,1]]

if band=='H':
    col=18
elif band=='Ks':
    col=20
field18_gal=[field18_gal1[i] for i in range(len(field18_gal1)) if field18_gal1[i,col]!=99]
field18_gal=np.array(field18_gal)


# In[4]:


l_max=max(field18_gal[:,22])
b_max=max(field18_gal[:,23])
l_min=min(field18_gal[:,22])
b_min=min(field18_gal[:,23])
lim=[l_max,b_max,l_min,b_min]


# In[5]:


total=[]
for chip in range(1,5):
    lst_1=np.loadtxt(tmp+name+'_chip%s.txt'%(chip))
    lst_1_comm=[lst_1[i] for i in range(len(lst_1)) if (l_min<lst_1[i,8]<l_max and b_min<lst_1[i,9]<b_max)]
    total.extend(lst_1_comm)
np.savetxt(tmp+name+'_on_GNS.txt',total)


# In[6]:


RA,DRA,DEC,DDEC,raJ,draJ,decJ,ddecJ,raH,draH,decH,ddecH,raKs,draKs,decKs,ddecKs,J,dJ,H,dH,Ks,dKs,l,b=np.split(field18_gal,24,axis=1)
GNS_mag=[H,dH,Ks,dKs]
GNS_x=[raH,draH,raKs,draKs]
GNS_y=[decH,ddecH,decKs,ddecKs]


# In[7]:


total=np.array(total)
ra,dec,x_mean,dx,y_mean,dy,mag,dmag,l,b=np.split(total,10, axis=1)


# In[8]:


if exptime==2:
    l=max(mag)+0.2
elif exptime==10:
    l=max(mag)+0.2
if band=='H':
    m=0
    c='k'
elif band=='Ks':
    m=2
    c='b'
ms=6
y_l=0.015
alpha=0.1

fig,ax=plt.subplots(1,2,figsize=(20,10))
ax[0].plot(GNS_mag[m],GNS_mag[m+1],'x',scalex=True,alpha=alpha,color='red',markersize=ms)

#ax[0].legend(['GNS',])
ax[0].set_xlabel('[%s]'%(band),fontsize=20)
ax[0].set_ylabel('d[%s]'%(band),fontsize=20)
leg=ax[0].legend(['GNS #=%s'%(len(GNS_mag[0]))],fontsize=20,markerscale=3,shadow=True,loc=2)
for lh in leg.legendHandles: 
    lh._legmarker.set_alpha(1)
ax[0].grid()
#ax[0].set_xlim(10,l)
ax[0].set_ylim(0,0.20)
fig.set_alpha(0.2)
ax[1].plot(mag,dmag,'x',scalex=True,alpha=alpha,color='orange',markersize=ms)
ax[1].set_xlabel('[%s]'%(band),fontsize=20)
ax[1].grid()
ax[1].set_xlim(10,l)
ax[1].set_ylim(0,0.20)
ax[1].set_ylabel('d[%s]'%(band),fontsize=20)
leg=ax[1].legend(['The Brick (DIT=%s, #=%s)'%(exptime,len(mag))],fontsize=20,markerscale=3,shadow=True)
for lh in leg.legendHandles: 
    lh._legmarker.set_alpha(1)


ax[0].tick_params(axis='both', which='major', labelsize=20)
ax[1].tick_params(axis='both', which='major', labelsize=20)
#plt.suptitle('The Brick`s stars = %s (DIT=%s). GNS` stars =%s'%(len(brick_mag),exptime,len(GNS_x[0])),fontsize=20)
plt.savefig(results + 'Brick_GNS_mag_pos_band%s.png'%(band))
plt.show()


# In[9]:


if band=='H':
    col=18
elif band=='Ks':
    col=20
xy = np.vstack([field18_gal[:,col],field18_gal[:,col+1]])
z = gaussian_kde(xy)(xy)
xy1 = np.vstack([total[:,6],total[:,7]])
z1= gaussian_kde(xy1)(xy1)


if exptime==2:
    l=max(mag)+0.2
    map_c='viridis'
elif exptime==10:
    l=max(mag)+0.2
    map_c='inferno'
if band=='H':
    m=0
    c='k'
elif band=='Ks':
    m=2
    c='b'
ms=6
y_l=0.015

fig,ax=plt.subplots(1,2,figsize=(20,10))
ax[0].scatter(GNS_mag[m],GNS_mag[m+1],c=z, s=5,alpha=1,cmap=map_c)

#ax[0].legend(['GNS',])
ax[0].set_xlabel('[%s]'%(band),fontsize=20)
ax[0].set_ylabel('d[%s]'%(band),fontsize=20)
ax[0].legend(['GNS #=%s'%(len(GNS_x[0]))],fontsize=20,markerscale=3,shadow=True,loc=2)
ax[0].grid()
ax[0].set_xlim(10,l)
ax[0].set_ylim(0,0.20)

ax[1].scatter(mag,dmag,c=z1, s=5,alpha=1,cmap=map_c)
ax[1].set_xlabel('[%s]'%(band),fontsize=20)
ax[1].grid()
ax[1].set_xlim(10,l)
ax[1].set_ylim(0,0.20)
ax[1].set_ylabel('d[%s]'%(band),fontsize=20)
ax[1].legend(['The Brick (DIT=%s, #=%s)'%(exptime,len(mag))],fontsize=20,markerscale=3,shadow=True)
ax[0].tick_params(axis='both', which='major', labelsize=20)
ax[1].tick_params(axis='both', which='major', labelsize=20)
#plt.suptitle('The Brick`s stars = %s (DIT=%s). GNS` stars =%s'%(len(brick_mag),exptime,len(GNS_x[0])),fontsize=20)
plt.savefig(results + 'Brick_GNS_mag_pos_band%s.png'%(band))
plt.show()


# In[10]:



if exptime==2:
    l=20
elif exptime==10:
    l=19.5
if band=='H':
    m=0
    c='k'
elif band=='Ks':
    m=2
    c='b'
ms=6
y_l=0.015
alpha=0.1
fig,ax=plt.subplots(1,2,figsize=(20,10))
ax[0].plot(GNS_mag[m],GNS_x[m+1],'x',scalex=True,alpha=alpha,color='red',markersize=ms)
#ax[0].legend(['GNS',])
ax[0].set_xlabel('[%s]'%(band),fontsize=20)
ax[0].set_ylabel('X Uncertainty (arcsec)',fontsize=20)
ax[0].grid()
ax[0].plot(mag,total[:,3],'.',scalex=True,alpha=alpha,color=c,markersize=ms)
leg=ax[0].legend(['GNS','The Brick'],fontsize=20,markerscale=3,shadow=True,loc=2)
for lh in leg.legendHandles: 
    lh._legmarker.set_alpha(1)

ax[0].set_xlim(10,max(mag))
ax[0].set_ylim(0,max(total[:,3]/2))


ax[1].plot(GNS_mag[m],GNS_y[m+1],'x',scalex=True,alpha=alpha,color='red',markersize=ms)
ax[1].set_xlabel('[%s]'%(band),fontsize=20)
ax[1].grid()
ax[1].plot(mag,dy,'.',scalex=True,alpha=alpha,color=c,markersize=ms)
ax[1].set_xlim(10,max(mag))
ax[1].set_ylim(0,max(total[:,5]/2))
ax[1].set_ylabel('Y Uncertainty (arcsec)',fontsize=20)
leg=ax[1].legend(['GNS','The Brick'],fontsize=20,markerscale=3,shadow=True)
for lh in leg.legendHandles: 
    lh._legmarker.set_alpha(1)

#plt.suptitle('The Brick`s stars = %s (DIT=%s). GNS` stars =%s'%(len(brick_mag),exptime,len(GNS_x[0])),fontsize=20)
plt.savefig(results + 'Brick_GNS_mag_pos_band%s.png'%(band))
plt.show()

ax[0].tick_params(axis='both', which='major', labelsize=20)
ax[1].tick_params(axis='both', which='major', labelsize=20)


# In[11]:


if exptime==2:
    l=20
    map_c='viridis'
elif exptime==10:
    map_c='inferno'
    l=20
if band=='Ks':
    m=2
    c='b'
    idi=[20,13,15]
elif band=='H':
    m=0
    idi=[18,9,11]
    c='k'

ms=6
y_l=max(dx)*0.50

alpha=0.1

#RA,DRA,DEC,DDEC,raJ,draJ,decJ,ddecJ,raH,draH,decH,ddecH,raKs,draKs,decKs,ddecKs,J,dJ,H,dH,Ks,dKs,l,b=np.split(field18_gal,24,axis=1)
#ra,dec,x_mean,dx,y_mean,dy,mag,dmag,l,b=np.split(total,10, axis=1)
xy = np.vstack([field18_gal[:,idi[0]],field18_gal[:,idi[1]]])
z = gaussian_kde(xy)(xy)
xy1 = np.vstack([total[:,6],total[:,3]])
z1= gaussian_kde(xy1)(xy1)

xy2 = np.vstack([field18_gal[:,idi[0]],field18_gal[:,idi[2]]])
z2 = gaussian_kde(xy2)(xy2)
xy3 = np.vstack([total[:,6],total[:,5]])
z3= gaussian_kde(xy3)(xy3)


fig,ax=plt.subplots(2,2,figsize=(20,20))
ax[0,0].scatter(GNS_mag[m],GNS_x[m+1],c=z, s=5,alpha=1,cmap=map_c)
#ax[0].legend(['GNS',])
ax[0,0].set_xlabel('[%s]'%(band),fontsize=20)
ax[0,0].set_ylabel('dra (arcsec)',fontsize=20)
#fig.colorbar(ax[0,0].scatter(GNS_mag[m],GNS_x[m+1],c=z, s=5,alpha=1,cmap=map_c), ax=ax[0,0])
ax[0,0].grid()

ax[0,0].set_xlim(10,l)
ax[0,0].set_ylim(0,y_l)

ax[0,1].scatter(GNS_mag[m],GNS_y[m+1],c=z2, s=5,alpha=1,cmap=map_c)
#ax[0].legend(['GNS',])
ax[0,1].set_xlabel('[%s]'%(band),fontsize=20)
ax[0,1].set_ylabel('ddec (arcsec)',fontsize=20)
#fig.colorbar(ax[0,1].scatter(GNS_mag[m],GNS_y[m+1],c=z2, s=5,alpha=1,cmap=map_c), ax=ax[0,1])
ax[0,1].grid()

ax[0,1].set_xlim(10,l)
ax[0,1].set_ylim(0,y_l)


ax[0,1].legend(['GNS'],fontsize=20,markerscale=3,shadow=True,loc=2)
ax[0,0].legend(['GNS #=%s'%(len(GNS_x[0]))],fontsize=20,markerscale=3,shadow=True,loc=2)
############################################################

ax[1,0].scatter(mag,dx,c=z1, s=5,alpha=1,cmap=map_c)
#ax[0].legend(['GNS',])
ax[1,0].set_xlabel('[%s]'%(band),fontsize=20)
ax[1,0].set_ylabel('X Uncertainty (arcsec)',fontsize=20)
#fig.colorbar(ax[1,0].scatter(brick_mag,std_x,c=z1, s=5,alpha=1,cmap=map_c), ax=ax[1,0])
ax[1,0].grid()

ax[1,0].set_xlim(10,l)
ax[1,0].set_ylim(0,y_l)

ax[1,1].scatter(mag,dy,c=z3, s=5,alpha=1,cmap=map_c)
#ax[0].legend(['GNS',])
ax[1,1].set_xlabel('[%s]'%(band),fontsize=20)
ax[1,1].set_ylabel('Y Uncertainty (arcsec)',fontsize=20)
#fig.colorbar(ax[1,1].scatter(brick_mag,std_y,c=z3, s=5,alpha=1,cmap=map_c), ax=ax[1,1])
ax[1,1].grid()

ax[1,1].set_xlim(10,l)
ax[1,1].set_ylim(0,y_l)
ax[1,0].legend(['The Brick (DIT=%s, #=%s)'%(exptime,len(mag))],fontsize=20,markerscale=3,shadow=True,loc=2)
ax[1,1].legend(['The Brick'],fontsize=20,markerscale=3,shadow=True,loc=2)
#plt.suptitle('The Brick`s stars = %s (DIT=%s). GNS` stars =%s'%(len(brick_mag),exptime,len(GNS_x[0])),fontsize=20)
ax[0,0].tick_params(axis='both', which='major', labelsize=20)
ax[1,0].tick_params(axis='both', which='major', labelsize=20)
ax[0,1].tick_params(axis='both', which='major', labelsize=20)
ax[1,1].tick_params(axis='both', which='major', labelsize=20)

#plt.savefig(results + 'Brick_GNS_mag_pos_band%s_separate.png'%(band))


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




