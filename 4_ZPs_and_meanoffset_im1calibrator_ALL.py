#!/usr/bin/env python
# coding: utf-8

# In[1]:



import glob
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import distance
from astropy.stats import sigma_clip



# In[2]:


#esto hace que se represnten las figuras en una nueva ventana.
#%matplotlib
#from matplotlib import interactive # para que se represeten la figuras
#interactive(True)


# In[3]:


band='H'
exptime=10
#chip=4
folder='im_jitter_NOgains/'
#folder='im_sky_ESOReflex/'
results='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder+'/results/'
pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/pruebas/'
indir = '/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder
psf='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder
tmp='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder+'tmp/'


# In[4]:

for chip in range(1,5):
    distancia=0.5
    
    
    # In[5]:
    
    
    #sub = sorted(glob.glob(indir+ 'stars_*chip'+str(chip)+'*.txt'), key=os.path.getmtime)
    sub = sorted(glob.glob(tmp+ 'cube_list_E_im*_chip'+str(chip)+'*aligned.txt'), key=os.path.getmtime)
    for i in sub:
        print(i)
    dic_stars={}
    s=0
    for i in sub:
        s+=1
        dic_stars['stars_im'+str(s)]=np.loadtxt(i)
    len(dic_stars)
    
    
    # In[6]:
    
    
    dic_listas={}
    for j in range(2,len(dic_stars)+1):
        diff=[]
        im_a=1
        im_b=j
        for i in range(dic_stars['stars_im'+str(im_a)].shape[0]): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
                    dist=distance.cdist(dic_stars['stars_im'+str(im_a)][i:i+1,0:2],dic_stars['stars_im'+str(im_b)][:,0:2], 'euclidean')
                    d=np.where(dist<=distancia)
                    if len(d[1])>0:
                        diff.append((dic_stars['stars_im'+str(im_a)][i],dic_stars['stars_im'+str(im_b)][d[1][np.argmin(dist[d])]]))
        dic_listas['stars_1'+str(j)]=diff
        print('comunes listas %s y %s ----> '%(im_a,im_b),len(diff))
    
    
    # In[7]:
    
    
    dic_mag1={} #magnitudes de las estrellas comunes de la primera lista con cada una de las otras
    if band=='H':
        ZP=26.32
    else:
        ZP=25.63
    for j in range(2,len(dic_stars)+1):
        mag1=[]
        for i in range(len(dic_listas['stars_1'+str(j)])):
            mag=ZP-2.5*np.log10(dic_listas['stars_1'+str(j)][i][0][2]/exptime)
            mag1.append(mag)
        dic_mag1['mag1_lis_1'+str(j)]=mag1
    
    
    # In[8]:
    
    
    dic_zp={}# zp cada una de las listas usando como calibrador la magintud de las estrellas comunes de la lista 1
    for j in range(2,len(dic_stars)+1):
        zp_im=[]
        for i in range(len(dic_listas['stars_1'+str(j)])):
            zp=dic_mag1['mag1_lis_1'+str(j)][i]+2.5*np.log10(dic_listas['stars_1'+str(j)][i][1][2]/exptime)
            zp_im.append(zp)
        print('zp_im sin clip %s'%(len(zp_im)))
        print(np.mean(zp_im),np.std(zp_im))
        zp_im=sigma_clip(zp_im,sigma=2,masked=False)
        dic_zp['zp_'+str(j)]=zp_im
        print('zp_im con clip %s'%(len(zp_im)))
        print(np.mean(zp_im),np.std(zp_im))
    
    
    
    # In[9]:
    
    
    zp_im=[]#zp de cada de las images y su error rms(Root-mean-square deviation)
    zp_e=[]
    for i in range(2,len(dic_stars)+1):
        zp_im.append(np.mean(dic_zp['zp_'+str(i)]))
        zp_e.append(np.std(dic_zp['zp_'+str(i)])/np.sqrt(len(dic_zp['zp_'+str(i)])))
        #print(np.mean(dic_zp['zp_'+str(i)]),np.std(dic_zp['zp_'+str(i)])/np.sqrt(len(dic_zp['zp_'+str(i)])))
        print(np.mean(dic_zp['zp_'+str(i)]),np.std(dic_zp['zp_'+str(i)]))
    if folder=='im_sky_ESOReflex/':
        np.savetxt(results+"zp_%s_ESO_chip%s.txt"%(band,chip), zp_im, fmt="%.3f")
    else:
        np.savetxt(results+"zp_%s_chip%s.txt"%(band,chip), zp_im, fmt="%.3f")
    
    
    # In[10]:
    
    
    dic_mag={}
    for j in range(2,len(dic_stars)+1):
        mag=[]
        for i in range(len(dic_listas['stars_1'+str(j)])):
            mag.append(zp_im[j-2]-2.5*np.log10(dic_listas['stars_1'+str(j)][i][1][2]/exptime))
        print(zp_im[j-2])
        dic_mag['mag_'+str(j)]=mag
    
    
    # In[11]:
    
    
    zp_im
    
    
    # In[12]:
    
    
    dic_diff={}
    for j in range(2,len(dic_stars)+1):
        b=[]
        for i in range(len(dic_mag['mag_'+str(j)])):
            a=dic_mag['mag_'+str(j)][i]-dic_mag1['mag1_lis_1'+str(j)][i]
            b.append(a)
        dic_diff['dmag_'+str(j)]=b
    
    
    # In[13]:
    
    
    if band=='H':
        x0=12
        x1=14
        x2=16
        x3=18
    else:
        x0=11
        x1=13
        x2=15
        x3=17
    for r in range(len(dic_stars)-1):
            fig, ax = plt.subplots(1,figsize=(10,10))
            ax.plot(dic_mag['mag_'+str(r+2)][:],dic_diff['dmag_'+str(r+2)][:],'x',scalex=True,alpha=0.4,color='k')
            ax.set_xlabel('[%s] im%s'%(band,r+2),fontsize=10)
            ax.set_ylabel('d[%s] '%(band),fontsize=10)
            
            #ax[r].ylabel('d[%s]'%(band),fontsize=20)
            #ax.set_ylim(-1,1)
            ax.set_title('Diff mag im1 and im%s, chip%s '%(r+2,chip),fontsize=10)
            magim1=np.array(dic_mag['mag_'+str(r+2)][:])
            rg0=np.where((magim1>x0)&(magim1<=x1))
            rg1=np.where((magim1>x1)&(magim1<=x2))
            rg2=np.where((magim1>x2)&(magim1<=x3))
            diffmag=np.array(dic_diff['dmag_'+str(r+2)][:])
            mean=np.mean(diffmag)
            mean0=np.mean(diffmag[rg0])
            std0=np.std(diffmag[rg0])
            mean1=np.mean(diffmag[rg1])
            std1=np.std(diffmag[rg1])
            mean2=np.mean(diffmag[rg2])
            std2=np.std(diffmag[rg2])
            ax.axhline(mean, color='g', linestyle='dashed', linewidth=3)
            ax.text(min(dic_mag['mag_'+str(r+2)][:])+1,1,'mean offset =%.2e'%(mean),color='g')
            ax.errorbar(x1-1,mean1,std0,color='red', elinewidth=3,capsize=10,capthick=2,barsabove=True,zorder=3)
            ax.errorbar(x2-1,mean2,std1,color='red',elinewidth=3,capsize=10,capthick=2,barsabove=True,zorder=3)
            ax.errorbar(x3-1,mean2,std2,color='red',elinewidth=3,capsize=10,capthick=2,barsabove=True,zorder=3)
            ax.text(x1-1.5,0.75,'%s<mg≤%s'%(x0,x1),color='r',zorder=3)
            ax.text(x2-1.5,+0.75,'%s<mg≤%s'%(x1,x2),color='r',zorder=3)
            ax.text(x3-1.5,0.75,'%s<mg≤%s'%(x2,x3),color='r',zorder=3)
            
            ax.text(x1-1.2,-0.75,'%.3f'%(std0),color='r',zorder=3)
            ax.text(x2-1.2,-0.75,'%.3f'%(std1),color='r',zorder=3)
            ax.text(x3-1.2,-0.75,'%.3f'%(std2),color='r',zorder=3)
            if folder =='im_sky_ESOReflex/':
                ax.text(x0,0.5,'ESO sky',fontsize=20,color='k')
                plt.savefig(pruebas + 'diff%s_im_%s%s_ESO_chip%s.png'%(im_a,band,r+2,chip))
            else:
                print('#####')
                plt.savefig(tmp + 'diff%s_im_%s%s_chip%s.png'%(im_a,band,r+2,chip))
    
    
    # In[14]:
    
    
    fig, ax = plt.subplots(1,figsize=(10,10))
    color=['k','g','r','b']
    for r in range(0,4):
        ax.plot(dic_mag['mag_'+str(r+2)][:],dic_diff['dmag_'+str(r+2)][:],'x',scalex=True,alpha=0.2,color=color[r])
        ax.set_xlabel('[%s] im%s'%(band,r+2),fontsize=10)
        #ax[r].ylabel('d[%s]'%(band),fontsize=20)
        ax.set_ylim(-1,1)
    
