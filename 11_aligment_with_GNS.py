#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
from astropy.table import QTable
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.spatial import distance
from astropy.wcs import WCS
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
from astropy.io import fits


# In[2]:


# como comprobacion has hecho este mismo scripts usando IDL con el metodo de los initial offset. El resultado es el mismo.
# el script de IDL se llama aligment_with_GNS.pro
band='Ks'
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
GNS='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/field_Bulge18/'
scripts='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/scripts/'

stream = open(scripts+'polywarp.py')
read_file = stream.read()
exec(read_file)


# In[3]:


distancia=1
if band =='H':
    r_gns=8
    d_gns=10
    brillo=18
if band=='Ks':
    r_gns=12
    d_gns=14
    brillo=20


# In[4]:


RA,DRA,DEC,DDEC,raJ,draJ,decJ,ddecJ,raH,draH,decH,ddecH,raKs,draKs,decKs,ddecKs,J,dJ,H,dH,Ks,dKs=np.loadtxt(GNS+'catalogue_F18_c3.txt',unpack=True)
field18=np.loadtxt(GNS+'catalogue_F18_c3.txt')
GNSH_gal = SkyCoord(ra=field18[:,0]*u.degree, dec=field18[:,2]*u.degree, frame='icrs').galactic
t = QTable([GNSH_gal], names=['coord'])
df=t.to_pandas()
gal_GNS=df.to_numpy()
field18=np.c_[field18,gal_GNS[:,0],gal_GNS[:,1]]


# In[5]:


l_max=max(field18[:,22])
l_min=min(field18[:,22])
b_max=max(field18[:,23])
b_min=min(field18[:,23])


# In[ ]:





# In[6]:


total=[]
for chip in range(1,5):
    #ra ,dec , m, dm, f, df,x,y = np.loadtxt(tmp+'stars_calibrated_%s_chip%s_sirius.txt'%(band,chip),unpack=True)
    lst_1=np.loadtxt(tmp+'stars_calibrated_%s_chip%s_sirius.txt'%(band,chip),unpack=False)
    brick_gal = SkyCoord(ra=lst_1[:,0]*u.degree, dec=lst_1[:,1]*u.degree, frame='icrs').galactic
    t = QTable([brick_gal], names=['coord'])
    df=t.to_pandas()
    gal_brick=df.to_numpy()
    lst_1=np.c_[lst_1,gal_brick[:,0],gal_brick[:,1]]
    lst_1_comm=[lst_1[i] for i in range(len(lst_1)) if (l_min<lst_1[i,8]<l_max and b_min<lst_1[i,9]<b_max)]
    total.extend(lst_1_comm)
    print('Chip %s galaticeado'%(chip))
np.savetxt(tmp+'calibrated_%s_on_GNS.txt'%(band),total)


# In[7]:


brick=np.loadtxt(tmp+'calibrated_%s_on_GNS.txt'%(band))


# In[ ]:





# In[8]:


#fig,ax=plt.subplots(2,1,figsize=(20,10))
#ax[0].scatter(field18[:,-2],field18[:,-1])
#ax[1].scatter(brick[:,6],brick[:,7])
#np.savetxt(pruebas+'birck_p.txt',brick)


# In[9]:


brick=np.loadtxt(tmp+'calibrated_%s_on_GNS.txt'%(band))
field18=np.loadtxt(GNS+'catalogue_F18_c3.txt')
print(len(field18))
valid=np.where(field18[:,r_gns]>0)
field18=field18[valid]
print(len(field18))
x_gns=[]
y_gns=[]
f = fits.open(GNS+'field18_c3.fits')
w = WCS(f[0].header)
print(w)
x_gns,y_gns=w.world_to_pixel_values(field18[:,r_gns],field18[:,d_gns])
field18=np.c_[field18,x_gns,y_gns]
#field18= RA,DRA,DEC,DDEC,raJ,draJ,decJ,ddecJ,raH,draH,decH,ddecH,raKs,draKs,decKs,ddecKs,J,dJ,H,dH,Ks,dKs,x,y
x_brick,y_brick=w.world_to_pixel_values(brick[:,0],brick[:,1])
brick[:,6]=x_brick
brick[:,7]=y_brick



for loop in range(10):
    diff=[]
    for i in range(len(field18)): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
        dist=distance.cdist(field18[i:i+1,-2:],brick[:,6:8], 'euclidean')
        d=np.where(dist<distancia)
        if len(d[1])>0:
            diff.append((field18[i],brick[d[1][np.argmin(dist[d])]]))
    print('comunes listas %s y %s ----> '%('brick','GNS'),len(diff))
    diff=np.array(diff)
    x1=[]
    x2=[]
    y1=[]
    y2=[]
    for j in range(len(diff)):
        x1.append(diff[j][0][-2])
        y1.append(diff[j][0][-1])
        x2.append(diff[j][1][6])
        y2.append(diff[j][1][7])

    x1=np.array(x1)
    y1=np.array(y1)
    x2=np.array(x2)
    y2=np.array(y2)

    Kx=[]
    Ky=[]
    degree=1
    Kx,Ky=polywarp(x1,y1,x2,y2,degree=degree)
    #print(Kx[0,0])
    xi=np.zeros(len(brick))
    yi=np.zeros(len(brick))
    x=[]
    y=[]
    x=brick[:,6]
    y=brick[:,7]
    x=np.array(x)
    y=np.array(y)
    for k in range(degree+1):
        for m in range(degree+1):
            xi=xi+Kx[k,m]*x**k*y**m
            yi=yi+Ky[k,m]*x**k*y**m
    brick[:,6]=xi
    brick[:,7]=yi


# In[10]:


#now secod degree alignment
for loop in range(10):
    diff=[]
    indx=[]
    for i in range(len(field18)): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
        dist=distance.cdist(field18[i:i+1,-2:],brick[:,6:8], 'euclidean')
        d=np.where(dist<distancia)
        if len(d[1])>0:
            indx.append(i)
            diff.append((field18[i],brick[d[1][np.argmin(dist[d])]]))
    print('comunes listas %s y %s ----> '%('brick','GNS'),len(diff))
    diff=np.array(diff)
    x1=[]
    x2=[]
    y1=[]
    y2=[]
    for j in range(len(diff)):
        x1.append(diff[j][0][-2])
        y1.append(diff[j][0][-1])
        x2.append(diff[j][1][6])
        y2.append(diff[j][1][7])

    x1=np.array(x1)
    y1=np.array(y1)
    x2=np.array(x2)
    y2=np.array(y2)

    Kx=[]
    Ky=[]
    degree=2
    Kx,Ky=polywarp(x1,y1,x2,y2,degree=degree)
    #print(Kx[0,0])
    xi=np.zeros(len(brick))
    yi=np.zeros(len(brick))
    x=[]
    y=[]
    x=brick[:,6]
    y=brick[:,7]
    x=np.array(x)
    y=np.array(y)
    for k in range(degree+1):
        for m in range(degree+1):
            xi=xi+Kx[k,m]*x**k*y**m
            yi=yi+Ky[k,m]*x**k*y**m
    brick[:,6]=xi
    brick[:,7]=yi


# In[11]:


ra_brick,dec_brick=w.pixel_to_world_values(brick[:,6],brick[:,7])
brick[:,0]=ra_brick
brick[:,1]=dec_brick
np.savetxt(tmp+'brick_%s_alig_with_gns.txt'%(band),brick)


# In[12]:


'''fig,ax=plt.subplots(1,2,figsize=(40,10))
ax[0].scatter(field18[:,-2],field18[:,-1])
ax[1].scatter(brick[:,6],brick[:,7])
#np.savetxt(pruebas+'birck_p.txt',brick)
'''


# In[13]:


#RA,DRA,DEC,DDEC,raJ,draJ,decJ,ddecJ,raH,draH,decH,ddecH,raKs,draKs,decKs,ddecKs,J,dJ,H,dH,Ks,dKs=np.loadtxt(GNS+'catalogue_F18_c3.txt',unpack=True)
#ra ,dec , m, dm, f, df,x,y,l,b=calibrated_%s_on_GNS.txt


# In[14]:


resta=np.array([diff[k][0][brillo]-diff[k][1][2] for k in range(len(diff))])
mags=np.array([(diff[k][0][brillo],diff[k][1][2])for k in range(len(diff))])
nbins=np.ceil(max(mags[:,0]))-np.floor(min(mags[:,1]))
mmag =np.zeros(shape=(int(nbins)))
mags_bins=np.floor((min(mags[:,1])+np.arange(nbins)))
sig_mag =np.zeros(shape=(int(nbins)))
resta_clip=sigma_clipped_stats(resta,sigma=2,maxiters=5)
mask_sig=sigma_clip(resta,sigma=2,maxiters=10)
nope=np.where(mask_sig.mask==True)
# In[15]:


mags


# In[16]:


fig,ax=plt.subplots(1,figsize=(20,10))
for j in range(int(nbins)):
        thisbin=np.where((mags[:,0]>mags_bins[j])&(mags[:,0]<=mags_bins[j]+1))
        vals = resta[thisbin]
        #nope_thisbin=np.where((mag_ref[nope]>mags_bins[j])&(mag_ref[nope]<mags_bins[j]+1))
        if len(vals)>1:
            bin_rej=sigma_clip(vals, sigma=2, maxiters=5,masked=True)
            sig_bin=sigma_clipped_stats(vals,sigma=2,maxiters=5)
            mmag[j]=sig_bin[0]
            sig_mag[j]=sig_bin[2]
        ax.errorbar(mags_bins[j]+0.5,mmag[j],sig_mag[j],color='red', elinewidth=3,capsize=10,capthick=2,barsabove=True,zorder=3)
        ax.text(mags_bins[j]+0.4,min(resta)+0.01,'%.3f'%(sig_mag[j]),color='red',fontsize=14)  
        ax.text(mags_bins[j]+0.4,-min(resta),'%.1f%%'%(len(vals)/len(diff)*100),color='blue',fontsize=14)
        #ax.text(mags_bins[j]+0.4,-0.9,'%s'%(len(nope_thisbin[0])),color='orange',fontsize=14)
#ax.text((min(mags[:,0])),max(resta)-0.5,'mean = %.3f, std=%.3f'%(np.mean(resta),np.std(resta)),weight='bold',color ='g',fontsize=20,zorder=3)  
ax.text((min(mags[:,0])),max(resta)-0.5,'2$\sigma$ mean = %.3f, std=%.3f'%(resta_clip[0],resta_clip[2]),weight='bold',color ='g',fontsize=20,zorder=3) 
ax.scatter(mags[:,0],resta,color='k',alpha=0.3)
ax.scatter(mags[:,0][nope],resta[nope],color='red',marker='x',s=100,alpha=0.3)
ax.legend(['#%s'%(len(diff))],handlelength=-1,markerscale=0,shadow=True,fontsize=20)
#ax.axhline(np.mean(resta),color='g',ls='--',lw=3,zorder=3)
ax.axhline(resta_clip[0],color='g',ls='--',lw=3,zorder=3)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.set_ylabel('$[%s]_{GNS}-[%s]_{Zoc}$'%(band,band),size=20)
ax.set_xlabel('$[%s]_{Zoc}$'%(band),size=20)
ax.grid()
ax.set_title('Red is the std of stars in bins 1mag width. Blue is the percetange of stars in that bin. Crossed out stars are out of 2$\sigma$',fontsize=15, ha='center')

# In[17]:


print(diff[0][0])
print(diff[0][1])


# In[18]:


# Ahora la diffrencia en la posicion
#RA,DRA,DEC,DDEC,raJ,draJ,decJ,ddecJ,raH,draH,decH,ddecH,raKs,draKs,decKs,ddecKs,J,dJ,H,dH,Ks,dKs=np.loadtxt(GNS+'catalogue_F18_c3.txt',unpack=True)
#ra ,dec , m, dm, f, df,x,y,l,b=calibrated_%s_on_GNS.txt

x_shift=[(diff[i][0][-2]-diff[i][1][6]) for i in range(len(diff))]
y_shift=[(diff[i][0][-1]-diff[i][1][7]) for i in range(len(diff))]

x_shift=np.array(x_shift)*0.106
y_shift=np.array(y_shift)*0.106
fig,ax= plt.subplots(1,2,figsize=(20,10))
ls=[x_shift,y_shift]
nam=['x diff (arcsec)','y diff (arcsec)']
#ax[0]=plt.suptitle('Sigma Threshold at reconstruct = %s, CHIP  %s'%(th,ch),fontsize=20)
#plt.clf()
for h in range(len(ls)):
    ax[h].hist(ls[h], bins=10,alpha=0.7, rwidth=0.85,color='g')
    ax[h].axvline(np.mean(ls[h]), color='r', linestyle='dashed', linewidth=3)
    ax[h].grid(axis='both', alpha=0.75)
    ax[h].legend(['Chip%s: mean= %.4f, std=%.4f'%(1,np.mean(ls[h]),np.std(ls[h]))],fontsize=20,markerscale=0,shadow=True,loc=3,handlelength=0)
    ax[h].set_xlabel(nam[h],fontsize=20)
    ax[h].set_ylabel('# stars',fontsize=20)
    ax[h].tick_params(axis='x', labelsize=20)
    ax[h].tick_params(axis='y', labelsize=20)
fig.text(0.5, 0, 'Difference in position for common stars to GNS, band %s'%(band),fontsize=20, ha='center')


# In[ ]:





# In[19]:


# La fotometria de H no sale muy bien.Vamos a comparar la fotometria con Sirius
diff_1=[]
for chip in range(1,5):
    #x_valid,y_valid,m_valid,dm_valid=np.loadtxt(tmp + 'VALID_SIRUS_on_' + band + '_chip' + str(chip) + '.txt',unpack=True)
    si_valid=np.loadtxt(tmp + 'VALID_SIRUS_on_' + band + '_chip' + str(chip) + '.txt')
    #ra ,dec , m, dm, f, df,x,y=np.loadtxt(tmp + 'stars_calibrated_' + band + '_chip' + str(chip) + '_sirius.txt',unpack=True)
    brick_si=np.loadtxt(tmp + 'stars_calibrated_' + band + '_chip' + str(chip) + '_sirius.txt')
    for i in range(len(si_valid)): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
        dist=distance.cdist(si_valid[i:i+1,0:2],brick_si[:,6:8], 'euclidean')
        d=np.where(dist<distancia)
        if len(d[1])>0:
            diff_1.append((si_valid[i],brick_si[d[1][np.argmin(dist[d])]]))
    print('comunes listas %s y %s ----> '%('SI','Brick'),len(diff_1))


# In[20]:


len(si_valid)


# In[21]:


resta_1=np.array([diff_1[k][0][2]-diff_1[k][1][2] for k in range(len(diff_1))])
mags_1=np.array([(diff_1[k][0][2],diff_1[k][1][2])for k in range(len(diff_1))])
nbins_1=np.ceil(max(mags_1[:,1]))-np.floor(min(mags_1[:,1]))
mmag_1 =np.zeros(shape=(int(nbins_1)))
mags_bins_1=np.floor((min(mags_1[:,1])+np.arange(nbins_1)))
sig_mag_1 =np.zeros(shape=(int(nbins_1)))
resta_clip_1=sigma_clipped_stats(resta_1,sigma=2,maxiters=5)
mask_sig_1=sigma_clip(resta_1,sigma=2,maxiters=10)
nope_1=np.where(mask_sig_1.mask==True)

# In[22]:


fig,ax=plt.subplots(1,figsize=(20,10))
for j in range(int(nbins_1)):
        thisbin_1=np.where((mags_1[:,0]>mags_bins_1[j])&(mags_1[:,0]<=mags_bins_1[j]+1))
        vals_1 = resta_1[thisbin_1]
        #nope_thisbin=np.where((mag_ref[nope]>mags_bins[j])&(mag_ref[nope]<mags_bins[j]+1))
        if len(vals_1)>1:
            bin_rej_1=sigma_clip(vals_1, sigma=2, maxiters=5,masked=True)
            sig_bin_1=sigma_clipped_stats(vals_1,sigma=2,maxiters=10)
            mmag_1[j]=sig_bin_1[0]
            sig_mag_1[j]=sig_bin_1[2]
        ax.errorbar(mags_bins_1[j]+0.5,mmag_1[j],sig_mag_1[j],color='red', elinewidth=3,capsize=10,capthick=2,barsabove=True,zorder=3)
        ax.text(mags_bins_1[j]+0.4,min(resta_1)-min(resta_1)*0.33,'%.3f'%(sig_mag_1[j]),color='red',fontsize=14,zorder=3)  
        ax.text(mags_bins_1[j]+0.4,min(resta_1)-min(resta_1)*0.4,'%.1f%%'%(len(vals_1)/len(diff_1)*100),color='blue',fontsize=14,zorder=3)
        #ax.text(mags_bins[j]+0.4,-0.9,'%s'%(len(nope_thisbin[0])),color='orange',fontsize=14)
ax.text((min(mags_1[:,0])+.6),max(resta_1)-max(resta_1)*0.5,'mean = %.3f, std=%.3f'%(np.mean(resta_1),np.std(resta_1)),weight='bold',color ='g',fontsize=20,zorder=3)    
ax.scatter(mags_1[:,0],resta_1,color='k',alpha=0.3)
ax.legend(['#%s'%(len(diff_1))],handlelength=-1,markerscale=0,shadow=True,fontsize=20)
ax.axhline(np.mean(resta_1),color='g',ls='--',lw=3,zorder=3)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.set_ylabel('$[%s]_{SI}-[%s]_{Zoc}$'%(band,band),size=20)
ax.set_xlabel('$[%s]_{Zoc}$'%(band),size=20)
ax.grid()
ax.set_title('Difference in magnitud for SIRIUS and Zoc. Red is the std of stars in bins 1mag width. Blue is the percetange of stars in that bin',fontsize=17, ha='center')

# In[ ]:





# In[ ]:
# for some reason the offset of the 2sigma clipped data is bigger
fig,ax=plt.subplots(1,figsize=(20,10))
for j in range(int(nbins_1)):
        thisbin_1=np.where((mags_1[:,0]>mags_bins_1[j])&(mags_1[:,0]<=mags_bins_1[j]+1))
        vals_1 = resta_1[thisbin_1]
        #nope_thisbin=np.where((mag_ref[nope]>mags_bins[j])&(mag_ref[nope]<mags_bins[j]+1))
        if len(vals_1)>1:
            bin_rej_1=sigma_clip(vals_1, sigma=2, maxiters=5,masked=True)
            sig_bin_1=sigma_clipped_stats(vals_1,sigma=2,maxiters=10)
            mmag_1[j]=sig_bin_1[0]
            sig_mag_1[j]=sig_bin_1[2]
        ax.errorbar(mags_bins_1[j]+0.5,mmag_1[j],sig_mag_1[j],color='red', elinewidth=3,capsize=10,capthick=2,barsabove=True,zorder=3)
        ax.text(mags_bins_1[j]+0.4,min(resta_1)-min(resta_1)*0.33,'%.3f'%(sig_mag_1[j]),color='red',fontsize=14,zorder=3)  
        ax.text(mags_bins_1[j]+0.4,min(resta_1)-min(resta_1)*0.4,'%.1f%%'%(len(vals_1)/len(diff_1)*100),color='blue',fontsize=14,zorder=3)
        #ax.text(mags_bins[j]+0.4,-0.9,'%s'%(len(nope_thisbin[0])),color='orange',fontsize=14)
ax.text((min(mags_1[:,0])+.6),max(resta_1)-max(resta_1)*0.5,'mean = %.3f, std=%.3f'%(resta_clip_1[0],resta_clip_1[2]),weight='bold',color ='g',fontsize=20,zorder=3)    
ax.scatter(mags_1[:,0],resta_1,color='k',alpha=0.3)
ax.scatter(mags_1[:,0][nope_1],resta_1[nope_1],color='red',marker='x',s=100,alpha=0.3)
ax.legend(['#%s'%(len(diff_1))],handlelength=-1,markerscale=0,shadow=True,fontsize=20)
ax.axhline(np.mean(resta_1),color='g',ls='--',lw=3,zorder=3)
ax.tick_params(axis='x', labelsize=20)
ax.tick_params(axis='y', labelsize=20)
ax.set_ylabel('$[%s]_{SI}-[%s]_{Zoc}$'%(band,band),size=20)
ax.set_xlabel('$[%s]_{Zoc}$'%(band),size=20)
ax.grid()
ax.set_title('Difference in magnitud for SIRIUS and Zoc. Red is the std of stars in bins 1mag width. Blue is the percetange of stars in that bin',fontsize=17, ha='center')



# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




