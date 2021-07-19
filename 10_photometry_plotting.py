#!/usr/bin/env python
# coding: utf-8

# In[2]:


import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from scipy.spatial import distance
from astropy.stats import sigma_clipped_stats
from astropy.stats import sigma_clip


# In[3]:


band='H'
exptime=10
#Uncoment one of these
#folder='im_dark/'
#folder='im_jitter_gains/'
pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/pruebas/'
folder='im_jitter_NOgains/'
py_pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/py_pruebas/'
images='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/07.1_Reduce_aligned/054_'+band+'/dit_'+str(exptime)+'/'+folder
tmp='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder+'tmp/'
results='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder+'/results/'

# In[6]:


#HAWKI dh vs H
fig,ax=plt.subplots(2,2,figsize=(20,20))
a=1
for i in range(0,2):
    for j in range(0,2):
        ra ,dec , m, dm, f, df,x,y,dx,dy = np.loadtxt(tmp+'stars_calibrated_%s_chip%s_sirius.txt'%(band,a),unpack=True)
        xy = np.vstack([m,dm])
        z = gaussian_kde(xy)(xy)
        
        ax[i,j].scatter(m,dm,c=z, s=5,alpha=1,cmap='inferno')
        ax[i,j].set_ylim(0,0.1)
        ax[i,j].grid()
        ax[i,j].set_xlabel('[%s]'%(band),fontsize=20)
        ax[i,j].set_ylabel('[d%s]'%(band),fontsize=20)
        ax[i,j].tick_params(direction='out', length=6, width=2, colors='k',
                       grid_color='k', grid_alpha=0.5)
        ax[i,j].tick_params(axis='x', labelsize=20)
        ax[i,j].tick_params(axis='y', labelsize=20)
        ax[i,j].legend(['Chip%s #%s'%(a,len(m))],fontsize=20,markerscale=0,shadow=True,loc=2,numpoints=1,handlelength=0,handletextpad=0)
        a=a+1

plt.savefig(results + 'dmag_vs_mag_align_with_SIRIUS.png')
# In[ ]:

sig=2
if band =='H':
    y1=27
    y0=26
    x0=11.5
    x1=15.5
if band =='Ks':
     y1=26
     y0=25
     x0=11
     x1=15 
# In[ ]:





# In[5]:


#HAWKI ZP with SIRUIS reference, common stars
#polts the stars used for ZP and crosses out the ones elminated for the sigmaclipping algoritm when calculating the mean
fig, ax=plt.subplots(4,1,figsize=(10,20))
for i in range(0,4):
    
    zp,mag_ref = np.loadtxt(tmp+'ZP_ref_sirius_%s_chip%s.txt'%(band,i+1),unpack=True)
    s=sigma_clipped_stats(zp,sigma=sig,maxiters=10)
    mask_sig=sigma_clip(zp,sigma=sig,maxiters=10)
    nope=np.where(mask_sig.mask==True)
    ax[i].scatter(mag_ref,zp,color='k',alpha=0.5)
    ax[i].legend(['Chip%s #%s'%(i+1,len(zp))],fontsize=20,markerscale=0,shadow=True,loc=2,handlelength=0,handletextpad=0)
    ax[i].scatter(mag_ref[nope],zp[nope],color='red',marker='x',s=100,alpha=0.7)
    ax[i].axhline(s[1],color='green',ls='--',lw=5)
    #ax[i].set_ylim(y0,y1)
    ax[i].set_xlim(x0,x1)
    ax[i].grid()
    ax[i].tick_params(axis='x', labelsize=20)
    ax[i].tick_params(axis='y', labelsize=20)
    
    ax[i].text(x1-2,min(zp),'ZP =%.3f $\pm$ %.3f'%(s[0],s[2]/np.sqrt(len(zp)-1)),fontsize='xx-large',color='green',zorder=3,weight='bold')
    fig.text(0.5, 0.08, '[%s]'%(band),fontsize=30, ha='center')
    fig.text(0.5, 0.06, 'mean ZP(2$\sigma$). Croosed out stars are out of 2$\sigma$',fontsize=20, ha='center')
    
    fig.text(-0, 0.5, 'ZP', va='center', rotation='vertical',fontsize=30)
    
plt.savefig(results + 'ZP_align_with_SIRIUS.png')
# In[5]:

'''
# Hawki calibrated stars
chip=1
ra ,dec , m, dm, f, df,x,y=np.loadtxt(tmp+'stars_calibrated_%s_chip%s_sirius.txt'%(band,chip),unpack=True)
haw_cal=np.loadtxt(tmp+'stars_calibrated_%s_chip%s_sirius.txt'%(band,chip))

x_si,y_si,m_si,dm_si=np.loadtxt(tmp+'ALL_SIRUS_on_%s_chip%s.txt'%(band,chip),unpack=True)
ALL_si=np.loadtxt(tmp+'ALL_SIRUS_on_%s_chip%s.txt'%(band,chip))

x_ref,y_ref,m_ref,dm_ref=np.loadtxt(tmp+'ref_sirius_%s_chip%s.txt'%(band,chip),unpack=True)
ALL_si_ref=np.loadtxt(tmp+'ref_sirius_%s_chip%s.txt'%(band,chip))

'''
# In[6]:

'''
distancia=1
diff=[]
for i in range(len(ra)): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
            dist=distance.cdist(haw_cal[i:i+1,6:8],ALL_si[:,0:2], 'euclidean')
            d=np.where(dist<distancia)
            if len(d[1])>0:
                diff.append((haw_cal[i],ALL_si[d[1][np.argmin(dist[d])]]))
print('comunes listas %s y %s ----> '%('HA','SIRIUS'),len(diff))
media=[]
mags=[]
for i in range(len(diff)):
    #print(diff[i][0][2]-diff[i][1][2])
    media.append(diff[i][0][2]-diff[i][1][2])
    mags.append((diff[i][0][2],diff[i][1][2]))
media=np.mean(np.array(media))
mags=(np.array(mags))
print(media)
'''

# In[7]:
#plots the difference in mage between ref stars in SIRIUS and Zocalli. indicates the std by bins of 1mg withs and indicates
#the number of stars in that bin and the number of stars in that bin rejected by the clipping algorithim ,

fig,ax = plt.subplots(4,1,figsize=(10,20))
for k in range(0,4):
    chip=k+1
    ra ,dec , m, dm, f, df,x,y,dx,dy=np.loadtxt(tmp+'stars_calibrated_%s_chip%s_sirius.txt'%(band,chip),unpack=True)
    haw_cal=np.loadtxt(tmp+'stars_calibrated_%s_chip%s_sirius.txt'%(band,chip))

    x_si,y_si,m_si,dm_si=np.loadtxt(tmp+'VALID_SIRUS_on_%s_chip%s.txt'%(band,chip),unpack=True)
    ALL_si=np.loadtxt(tmp+'VALID_SIRUS_on_%s_chip%s.txt'%(band,chip))

    x_ref,y_ref,m_ref,dm_ref=np.loadtxt(tmp+'ref_sirius_%s_chip%s.txt'%(band,chip),unpack=True)
    ALL_si_ref=np.loadtxt(tmp+'ref_sirius_%s_chip%s.txt'%(band,chip))
    
    
    zp,mag_ref = np.loadtxt(tmp+'ZP_ref_sirius_%s_chip%s.txt'%(band,chip),unpack=True)
    s=sigma_clipped_stats(zp,sigma=sig,maxiters=10)
    mask_sig=sigma_clip(zp,sigma=sig,maxiters=10)
    nope=np.where(mask_sig.mask==True)
    
    distancia=1
    diff=[]
    for i in range(len(ra)): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
            dist=distance.cdist(haw_cal[i:i+1,6:8],ALL_si_ref[:,0:2], 'euclidean')
            d=np.where(dist<distancia)
            if len(d[1])>0:
                diff.append((haw_cal[i],ALL_si_ref[d[1][np.argmin(dist[d])]]))
    print('comunes listas %s y %s ----> '%('HA','SIRIUS'),len(diff))
    resta=[]
    mags=[]
    nbins=(np.ceil(max(m_ref))-np.floor(min(m_ref)))
    mags_bins=np.floor(min(m_ref))+np.arange(nbins)
    print(mags_bins)
    for i in range(len(diff)):
        #print(diff[i][0][2]-diff[i][1][2])
        resta.append(diff[i][0][2]-diff[i][1][2])
        mags.append((diff[i][0][2],diff[i][1][2]))
    diff=np.array(diff)
    resta=np.array(resta)
    average=np.mean(np.array(resta))
    mags=(np.array(mags))
    mmag =np.zeros(shape=(int(nbins)))
    sig_mag =np.zeros(shape=(int(nbins)))
    for j in range(int(nbins)):
        thisbin=np.where((mag_ref>mags_bins[j])&(mag_ref<=mags_bins[j]+1))
        vals = resta[thisbin]
        nope_thisbin=np.where((mag_ref[nope]>mags_bins[j])&(mag_ref[nope]<mags_bins[j]+1))
        if len(vals)>1:
            bin_rej=sigma_clip(vals, sigma=sig, maxiters=5,masked=True)
            sig_bin=sigma_clipped_stats(vals,sigma=sig,maxiters=5)
            mmag[j]=sig_bin[0]
            sig_mag[j]=sig_bin[2]
        ax[k].errorbar(mags_bins[j]+0.5,mmag[j],sig_mag[j],color='red', elinewidth=3,capsize=10,capthick=2,barsabove=True,zorder=3)
        ax[k].text(mags_bins[j]+0.4,-0.5,'%.3f'%(sig_mag[j]),color='red',fontsize=14)  
        ax[k].text(mags_bins[j]+0.4,-0.7,'%s'%(len(vals)),color='blue',fontsize=14)
        ax[k].text(mags_bins[j]+0.4,-0.9,'%s'%(len(nope_thisbin[0])),color='orange',fontsize=14)
        
    ax[k].scatter(mags[:,1],mags[:,0]-mags[:,1],color='k',alpha=0.2)
    ax[k].legend(['Chip%s #%s'%(chip,len(mags))],fontsize=20,markerscale=0,shadow=True,loc=1, handlelength=-1)
    #for lh in leg.legendHandles: 
     #   lh.set_visible(False)
    ax[k].axhline(average,color='g',ls='--',lw=2)
    ax[k].text(min(mags[:,1]),average+0.5,'mean offset (%s$\sigma$) =%.3f, std= %.3f'%(sig,average,np.std(resta)),color='green',fontsize=14,zorder=3,weight='bold') 
    ax[k].grid()
    ax[k].tick_params(axis='x', labelsize=20)
    ax[k].tick_params(axis='y', labelsize=20)
    ax[k].set_ylim(-1,1)
    #ax[k].set_xlim(x0,x1)
    #ax[k].text(12,0.5,'ZP =%.3f $\pm$ %.3f'%(s[0],s[2]/np.sqrt(len(diff)-1)),fontsize='xx-large',color='g') 
    fig.text(0.5, 0.08, '$[%s]_{Zoc}$'%(band),fontsize=30, ha='center')
    fig.text(-0, 0.5, '$[%s]_{Zoc}-[%s]_{SIRref}$'%(band,band), va='center', rotation='vertical',fontsize=30)
    fig.text(0.5, 0.06, 'Red is std of stars in bins of 1mag width. Blue is #stars in that bin and orage is #stars out of 2$\sigma$',fontsize=12, ha='center')


#%%
# In[7]:
#This plot is not ready yet

sig=3
fig,ax = plt.subplots(4,1,figsize=(10,20))
for k in range(0,4):
    chip=k+1
    ra ,dec , m, dm, f, df,x,y,dx,dy=np.loadtxt(tmp+'stars_calibrated_%s_chip%s_sirius.txt'%(band,chip),unpack=True)
    haw_cal=np.loadtxt(tmp+'stars_calibrated_%s_chip%s_sirius.txt'%(band,chip))

    x_si,y_si,m_si,dm_si=np.loadtxt(tmp+'VALID_SIRUS_on_%s_chip%s.txt'%(band,chip),unpack=True)
    ALL_si=np.loadtxt(tmp+'VALID_SIRUS_on_%s_chip%s.txt'%(band,chip))

    x_ref,y_ref,m_ref,dm_ref=np.loadtxt(tmp+'ref_sirius_%s_chip%s.txt'%(band,chip),unpack=True)
    ALL_si_ref=np.loadtxt(tmp+'ref_sirius_%s_chip%s.txt'%(band,chip))
    
    
   
    
    distancia=1
    diff=[]
    for i in range(len(ra)): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
            dist=distance.cdist(haw_cal[i:i+1,6:8],ALL_si[:,0:2], 'euclidean')
            d=np.where(dist<distancia)
            if len(d[1])>0:
                diff.append((haw_cal[i],ALL_si[d[1][np.argmin(dist[d])]]))
    print('comunes listas %s y %s ----> '%('HA','SIRIUS'),len(diff))
    resta=[]
    mags=[]
    nbins=(np.ceil(max(m_si))-np.floor(min(m_si)))
    mags_bins=np.floor(min(m_si))+np.arange(nbins)
    print(mags_bins)
    for i in range(len(diff)):
        #print(diff[i][0][2]-diff[i][1][2])
        resta.append(diff[i][0][2]-diff[i][1][2])
        mags.append((diff[i][0][2],diff[i][1][2]))
    diff=np.array(diff)
    resta=np.array(resta)
    average=np.mean(np.array(resta))
    mags=(np.array(mags))
    
     
    s=sigma_clipped_stats(mags[:,0],sigma=sig,maxiters=10)
    # mask_sig=sigma_clip(mags[:,0],sigma=sig,maxiters=10)
    mask_sig=sigma_clip(resta,sigma=sig,maxiters=10)
    nope=np.where(mask_sig.mask==True)
    
    
    mmag =np.zeros(shape=(int(nbins)))
    sig_mag =np.zeros(shape=(int(nbins)))
    for j in range(int(nbins)):
        thisbin=np.where((mags[:,0]>mags_bins[j])&(mags[:,0]<=mags_bins[j]+1))
        vals = resta[thisbin]
        nope_thisbin=np.where((mags[:,0][nope]>mags_bins[j])&(mags[:,0][nope]<mags_bins[j]+1))
        if len(vals)>1:
            bin_rej=sigma_clip(vals, sigma=sig, maxiters=5,masked=True)
            sig_bin=sigma_clipped_stats(vals,sigma=sig,maxiters=5)
            mmag[j]=sig_bin[0]
            sig_mag[j]=sig_bin[2]
            ax[k].errorbar(mags_bins[j]+0.5,mmag[j],sig_mag[j],color='red', elinewidth=3,capsize=10,capthick=2,barsabove=True,zorder=3)
            ax[k].text(mags_bins[j]+0.4,-0.5-0.5,'%.3f'%(sig_mag[j]),color='red',fontsize=14,zorder=3)  
            ax[k].text(mags_bins[j]+0.4,-0.7-0.5,'%s'%(len(vals)),color='blue',fontsize=14)
            ax[k].text(mags_bins[j]+0.4,-0.9-0.5,'%s'%(len(nope_thisbin[0])),color='orange',fontsize=14)
        
    ax[k].scatter(mags[:,0],mags[:,0]-mags[:,1],color='k',alpha=0.2)
    ax[k].scatter(mags[:,0][nope],resta[nope],color='orange',marker='x',s=100,alpha=0.7)
    ax[k].legend(['Chip%s #%s'%(chip,len(mags))],fontsize=20,markerscale=0,shadow=True,loc=1, handlelength=-1)
    #for lh in leg.legendHandles: 
     #   lh.set_visible(False)
    ax[k].axhline(average,color='g',ls='--',lw=2)
    ax[k].text(min(mags[:,0]),0.75,'mean offset (%s$\sigma$) =%.3f, std= %.3f'%(sig,average,np.std(resta)),color='green',fontsize=14,zorder=3,weight='bold') 
    ax[k].grid()
    ax[k].tick_params(axis='x', labelsize=20)
    ax[k].tick_params(axis='y', labelsize=20)
    ax[k].set_ylim(-1.5,1.5)
    #ax[k].set_xlim(x0,x1)
    #ax[k].text(12,0.5,'ZP =%.3f $\pm$ %.3f'%(s[0],s[2]/np.sqrt(len(diff)-1)),fontsize='xx-large',color='g') 
    fig.text(0.5, 0.08, '$[%s]_{Zoc}$'%(band),fontsize=30, ha='center')
    fig.text(-0, 0.5, '$[H]_{Zoc}-[H]_{SIRref}$', va='center', rotation='vertical',fontsize=30)
    fig.text(0.5, 0.06, 'Red is std of stars in bins of 1mag width. Blue is #stars in that bin and orage is #stars out of 2$\sigma$',fontsize=12, ha='center')
   
# In[ ]:
'''
chip=4
ra ,dec , m, dm, f, df,x,y,dx,dy=np.loadtxt(tmp+'stars_calibrated_%s_chip%s_sirius.txt'%(band,chip),unpack=True)
haw_cal=np.loadtxt(tmp+'stars_calibrated_%s_chip%s_sirius.txt'%(band,chip))

x_si,y_si,m_si,dm_si=np.loadtxt(tmp+'VALID_SIRUS_on_%s_chip%s.txt'%(band,chip),unpack=True)
ALL_si=np.loadtxt(tmp+'VALID_SIRUS_on_%s_chip%s.txt'%(band,chip))

x_ref,y_ref,m_ref,dm_ref=np.loadtxt(tmp+'ref_sirius_%s_chip%s.txt'%(band,chip),unpack=True)
ALL_si_ref=np.loadtxt(tmp+'ref_sirius_%s_chip%s.txt'%(band,chip))

zp,mag_ref = np.loadtxt(tmp+'ZP_ref_sirius_%s_chip%s.txt'%(band,chip),unpack=True)
s=sigma_clipped_stats(zp,sigma=3.0,maxiters=10)

distancia=0.5
diff=[]
for i in range(len(ra)): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
            dist=distance.cdist(haw_cal[i:i+1,6:8],ALL_si[:,0:2], 'euclidean')
            d=np.where(dist<distancia)
            if len(d[1])>0:
                diff.append((haw_cal[i],ALL_si[d[1][np.argmin(dist[d])]]))
print('comunes listas %s y %s ----> '%('HA','SIRIUS'),len(diff))
media=[]
mags=[]
for i in range(len(diff)):
    #print(diff[i][0][2]-diff[i][1][2])
    media.append(diff[i][0][2]-diff[i][1][2])
    mags.append((diff[i][0][2],diff[i][1][2]))
media=np.mean(np.array(media))
mags=(np.array(mags))
print(media)

x_shift=[(diff[i][0][6]-diff[i][1][0]) for i in range(len(diff))]
y_shift=[(diff[i][0][7]-diff[i][1][1]) for i in range(len(diff))]

x_shift=np.array(x_shift)*0.106
y_shift=np.array(y_shift)*0.106
fig,ax= plt.subplots(1,2,figsize=(20,10))
ls=[x_shift,y_shift]
nam=['x diff (arcsec)','y diff (arcsec)']
#ax[0]=plt.suptitle('Sigma Threshold at reconstruct = %s, CHIP  %s'%(th,ch),fontsize=20)
#plt.clf()
for h in range(len(ls)):
    ax[h].hist(ls[h], bins='auto',alpha=0.7, rwidth=0.85,color='g')
    ax[h].axvline(np.mean(ls[h]), color='r', linestyle='dashed', linewidth=3)
    ax[h].grid(axis='both', alpha=0.75)
    ax[h].legend(['Chip%s: mean= %.4f, std=%.3f'%(chip,np.mean(ls[h]),np.std(ls[h]))],fontsize=20,markerscale=0,shadow=True,loc=3,handlelength=0)
    ax[h].set_xlabel(nam[h],fontsize=20)
    ax[h].set_ylabel('# stars',fontsize=20)
    ax[h].tick_params(axis='x', labelsize=20)
    ax[h].tick_params(axis='y', labelsize=20)
fig.text(0.5, 0, 'Difference in position for common stars to SIRIUS, band %s'%(band),fontsize=20, ha='center')
'''
'''
ax[1].hist(x_shift, bins=10,alpha=0.7, rwidth=0.85,color='g')
ax[1].axvline(np.mean(x_shift), color='r', linestyle='dashed', linewidth=3)
ax[1].grid(axis='both', alpha=0.75)
ax[1].legend(['Chip%s: mean= %.4f, std=%.3f'%(1,np.mean(y_shift),np.std(y_shift))],fontsize=20,markerscale=0,shadow=True,loc=3,handlelength=0)
ax[1].set_xlabel('diff y (arcsec)',fontsize=20)
ax[1].set_ylabel('# stars',fontsize=20)
'''

