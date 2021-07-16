#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 09:14:33 2021

@author: amartinez
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
from scipy.stats import gaussian_kde
import pint
#%%
band='H'
chip=3
#band='Ks'
exptime=10
#Uncoment one of these
#folder='im_dark/'
#folder='im_jitter_gains/'
pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/pruebas/'
folder='im_jitter_NOgains/'
py_pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/py_pruebas/'
images='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/07.1_Reduce_aligned/054_'+band+'/dit_'+str(exptime)+'/'+folder
tmp='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder+'tmp/'
sirius='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/SIRIUS/'
GNS='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/field12/'
scripts='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/scripts/'

ureg = pint.UnitRegistry()# to give units to values
# 'x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK,H-Ks'
GNS=np.loadtxt(tmp+'GNS_commons_w_Zoc_c%s.txt'%(chip))



#%%

dist=8*ureg.kpc
s=3
unc_xy=0.006
unc_v=1


### 'a ,d , m, dm, f, df,x,y,dx,dy,x_displacement,y_displacement.
stars=np.loadtxt(tmp+'Zoc_c%s_commons_w_GNS.txt'%(chip))

dx_dis=np.sqrt((GNS[:,1]*0.106)**2+stars[:,8]**2)
dy_dis=np.sqrt((GNS[:,3]*0.106)**2+stars[:,9]**2)
dxy=np.sqrt(dx_dis**2+dy_dis**2)
stars=np.c_[stars,dxy]

### 'a ,d , m, dm, f, df,x,y,dx,dy,x_displacement,y_displacement,dxy



vx=stars[:,10]/4
vy=stars[:,11]/4
dvel=np.sqrt(((vx*dx_dis/4)**2+(vy*dy_dis/4)**2)/(vx**2+vy**2))

vel=np.sqrt(stars[:,10]**2+stars[:,11]**2)
stars=np.c_[stars,vel,dvel]
vel_clip=sigma_clip(vel,sigma=s)
stars=stars[vel_clip.mask==False]
#'a ,d , m, dm, f, df,x,y,dx,dy,x_displacement,y_displacement,dxy,vel,dvel

low_xy=np.where(stars[:,12]<unc_xy)
stars=stars[low_xy]

low_v=np.where(stars[:,14]<unc_v)
stars=stars[low_v]
vel_x=stars[:,10]*0.106
vel_y=stars[:,11]*0.106

vel_x=vel_x/(4*365*24*3600)
vel_y=vel_y/(4*365*24*3600)

vel_x=vel_x*ureg.arcsec
vel_y=vel_y*ureg.arcsec
dist=dist.to('km')
vel_x=vel_x.to('rad')*dist
vel_y=vel_y.to('rad')*dist

fig,ax= plt.subplots(1,2,figsize=(20,10))
ls=[np.array(vel_x),np.array(vel_y)]

nam=['$v_{x}$(km/s)','$v_{y}$(km/s)']
#ax[0]=plt.suptitle('Sigma Threshold at reconstruct = %s, CHIP  %s'%(th,ch),fontsize=20)
#plt.clf()

for h in range(len(ls)):
    
    his=ax[h].hist(ls[h], bins=10,alpha=0.7, rwidth=1,color='blue',edgecolor='black',linewidth=2)
    ax[h].axvline(np.mean(ls[h]), color='r', linestyle='dashed', linewidth=3)
    ax[h].grid(axis='both', alpha=0.75)
    ax[h].legend(['Chip=%s, #%s, mean= %.2f, std=%.2f'%(chip,len(ls[h]),np.mean(ls[h]),np.std(ls[h]))],fontsize=15,markerscale=0,shadow=True,loc=3,handlelength=-0.0)
    ax[h].set_xlabel(nam[h],fontsize=20)
    ax[h].set_ylabel('# stars',fontsize=20)
    ax[h].tick_params(axis='x', labelsize=20)
    ax[h].tick_params(axis='y', labelsize=20)
    ax[h].text(his[1][0],max(his[0]/2),'%s$\sigma$ clipped velocity'%(s),color='k',fontsize=20,zorder=3,weight='bold') 
    # ax[h].text(his[1][0],max(his[0]/2-his[0]/10),r'$\sigma_{\vec {v}}$(ars/yr)<%s'%(unc/4),color='k',fontsize=20,zorder=3,weight='bold') 
    if unc_v<unc_xy:
        ax[h].text(his[1][0],max(his[0]/2-his[0]/10),r'$\sigma_{\vec {v}}$<%s(ars/yr)'%(unc_v),color='k',fontsize=20,zorder=3,weight='bold') 
    else:
        ax[h].text(his[1][0],max(his[0]/2-his[0]/10),r'$\sigma_{\vec {xy}}$<%s"'%(unc_xy),color='k',fontsize=20,zorder=3,weight='bold') 


#fig.text(0.5, 0, 'Difference in position for common stars to GNS, band %s.'%(band),fontsize=20, ha='center'1.4

#%%
if unc_v>unc_xy:
    #### Uncertainty in POSITION
    ### 'a ,d , m, dm, f, df,x,y,dx,dy,x_displacement,y_displacement.
    
    stars=np.loadtxt(tmp+'Zoc_c%s_commons_w_GNS.txt'%(chip))
    ### 'x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK,H-Ks'
    GNS=np.loadtxt(tmp+'GNS_commons_w_Zoc_c%s.txt'%(chip))
    
    dx_dis=np.sqrt((GNS[:,1]*0.106)**2+stars[:,8]**2)
    dy_dis=np.sqrt((GNS[:,3]*0.106)**2+stars[:,9]**2)
    dxy=np.sqrt(dx_dis**2+dy_dis**2)
    
    menor=0
    for i in range(len(dxy)):
        if dxy[i]<unc_xy:
            menor+=1
    # stars=np.c_[stars,dxy]
    ### 'a ,d , m, dm, f, df,x,y,dx,dy,x_displacement,y_displacement.
    fig, ax=plt.subplots(1,1,figsize=(10,10))
    ax.scatter(stars[:,2],dxy,color='k',alpha=0.3)
    ax.set_xlabel('[H]',fontsize=20)
    ax.set_ylabel(r'$\sigma_{\vec {xy}}$(arcsec)',fontsize=20)
    # ax.set_ylim(0,0.03)
    # ax.axhline(l_min, color='r', linestyle='dashed', linewidth=3)
    # ax.axhline(l_max, color='r', linestyle='dashed', linewidth=3)
    ax.axhline(unc_xy, color='g', linestyle='dashed', linewidth=3,zorder=3)
    ax.text(min(stars[:,2]),unc_xy+unc_xy/10,'#stars= %s'%(menor), color='g',fontsize=20,weight='bold',zorder=3)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    

else:
    #### Uncertainty in VELOCITY
    ### 'a ,d , m, dm, f, df,x,y,dx,dy,x_displacement,y_displacement.
    
    stars=np.loadtxt(tmp+'Zoc_c%s_commons_w_GNS.txt'%(chip))
    ### 'x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK,H-Ks'
    GNS=np.loadtxt(tmp+'GNS_commons_w_Zoc_c%s.txt'%(chip))
    
    dx_dis=np.sqrt((GNS[:,1]*0.106)**2+stars[:,8]**2)
    dy_dis=np.sqrt((GNS[:,3]*0.106)**2+stars[:,9]**2)
    dxy=np.sqrt(dx_dis**2+dy_dis**2)
    
    vel=np.sqrt(stars[:,10]**2+stars[:,11]**2)
    vx=stars[:,10]/4
    vy=stars[:,11]/4
    
    dvel=np.sqrt(((vx*dx_dis/4)**2+(vy*dy_dis/4)**2)/(vx**2+vy**2))
    menor=0
    for i in range(len(dvel)):
        if dvel[i]<unc_v/4:
            menor+=1
    # stars=np.c_[stars,dxy]
    ### 'a ,d , m, dm, f, df,x,y,dx,dy,x_displacement,y_displacement.
    fig, ax=plt.subplots(1,1,figsize=(10,10))
    ax.scatter(stars[:,2],dvel,color='k',alpha=0.3)
    ax.set_xlabel('[H]',fontsize=20)
    ax.set_ylabel(r'$\sigma_{\vec {v}}$(arcsec/yr)',fontsize=20)
    # ax.set_ylim(0,0.03)
    # ax.axhline(l_min, color='r', linestyle='dashed', linewidth=3)
    # ax.axhline(l_max, color='r', linestyle='dashed', linewidth=3)
    ax.axhline(unc_v/4, color='g', linestyle='dashed', linewidth=3,zorder=3)
    ax.text(min(stars[:,2]),unc_v/4+unc_v/40,'#stars= %s'%(menor), color='g',fontsize=20,weight='bold',zorder=3)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    
    
    
    #%% 
############# CMD
# # 'x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK'
# chip=3
# GNS=np.loadtxt(tmp+'GNS_commons_w_Zoc_c%s.txt'%(chip))
# fig, ax = plt.subplots(1,1,figsize=(10,10))
# good=np.where(GNS[:,12]<90)
# GNS=GNS[good]
# print(len(GNS))
# menor=0
# mayor=0
# for i in range(len(GNS)):
#     if GNS[i][10]-GNS[i][12]<1.3:
#         menor+=1
#     elif GNS[i][10]-GNS[i][12]>=1.3:
#         mayor+=1
# ax.scatter(GNS[:,10]-GNS[:,12],GNS[:,12],color='k',alpha=0.3)
# ax.tick_params(axis='x', labelsize=20)
# ax.tick_params(axis='y', labelsize=20)
# ax.set_xlabel('[H]-[Ks]',fontsize=20)
# ax.set_ylabel('[Ks]',fontsize=20)
# ax.axvline(1.3, color='r', linestyle='dashed', linewidth=3)
# ax.text(0,10,'stars = %s'%(menor),color='k',fontsize=20,zorder=3,weight='bold')
# ax.text(1.5,10,'stars = %s'%(mayor),color='k',fontsize=20,zorder=3,weight='bold')      
# ax.invert_yaxis()







