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
band='H'
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

distancia=1
#l_mag=0.01


#slect apropiate columns from GNS
if band =='H':
    r_gns=8
    d_gns=10
    brillo=18
if band=='Ks':
    r_gns=12
    d_gns=14
    brillo=20


# In[3]:


#galactic coordinates for GNS
RA,DRA,DEC,DDEC,raJ,draJ,decJ,ddecJ,raH,draH,decH,ddecH,raKs,draKs,decKs,ddecKs,J,dJ,H,dH,Ks,dKs=np.loadtxt(GNS+'catalogue_F18_c3.txt',unpack=True)
field18=np.loadtxt(GNS+'catalogue_F18_c3.txt')
GNSH_gal = SkyCoord(ra=field18[:,0]*u.degree, dec=field18[:,2]*u.degree, frame='icrs').galactic
t = QTable([GNSH_gal], names=['coord'])
df=t.to_pandas()
gal_GNS=df.to_numpy()
field18=np.c_[field18,gal_GNS[:,0],gal_GNS[:,1]]
#field18=RA,DRA,DEC,DDEC,raJ,draJ,decJ,ddecJ,raH,draH,decH,ddecH,raKs,draKs,decKs,ddecKs,J,dJ,H,dH,Ks,dKs,l,b
l_max=max(field18[:,22])
l_min=min(field18[:,22])
b_max=max(field18[:,23])
b_min=min(field18[:,23])

#galactic coordinate for the Zoccaly to set them on field18
fig,ax=plt.subplots(1,1,figsize=(20,10))
dic_chips={}
for chip in range(1,5):
    #ra ,dec , m, dm, f, df,x,y,dx,dy= np.loadtxt(tmp+'stars_calibrated_%s_chip%s_sirius.txt'%(band,chip),unpack=True)
    lst_1=np.loadtxt(tmp+'stars_calibrated_%s_chip%s_sirius.txt'%(band,chip),unpack=False)
    print(len(lst_1),chip)
    brick_gal = SkyCoord(ra=lst_1[:,0]*u.degree, dec=lst_1[:,1]*u.degree, frame='icrs').galactic
    t = QTable([brick_gal], names=['coord'])
    df=t.to_pandas()
    gal_brick=df.to_numpy()
    lst_1=np.c_[lst_1,gal_brick[:,0],gal_brick[:,1]]
    #ra ,dec , m, dm, f, df,x,y,dx,dy,l,b
    lst_1_comm=[lst_1[i] for i in range(len(lst_1)) if (l_min<lst_1[i,10]<l_max and b_min<lst_1[i,11]<b_max)]
    dic_chips['chip_%s'%(chip)]=lst_1_comm
    dic_chips['chip_%s'%(chip)]=np.array( dic_chips['chip_%s'%(chip)])
    print('Chip %s galaticeado.'%(chip),len(dic_chips['chip_%s'%(chip)]))
    f = fits.open(GNS+'field18_c3.fits')
    w = WCS(f[0].header)
    #print(w)
    #extract coordinates X,Y within the frame of the image field18
    x_gns,y_gns=w.world_to_pixel_values(dic_chips['chip_%s'%(chip)][:,0],dic_chips['chip_%s'%(chip)][:,1])
    dic_chips['chip_%s'%(chip)][:,6]=x_gns
    dic_chips['chip_%s'%(chip)][:,7]=y_gns
    #field18= RA,DRA,DEC,DDEC,raJ,draJ,decJ,ddecJ,raH,draH,decH,ddecH,raKs,draKs,decKs,ddecKs,J,dJ,H,dH,Ks,dKs,x,y
    ax.scatter(dic_chips['chip_%s'%(chip)][:,6],(dic_chips['chip_%s'%(chip)][:,7]),alpha=0.2)
    np.savetxt(pruebas+'chunk_chip%s.txt'%(chip+1),dic_chips['chip_%s'%(chip)])


# In[4]:


#eliminates stars overlapping borders of chips(stars closer than 0.5'')
diff=[]
distancia=5
index=[]
c=1
for b in range(c,5):
    for j in range(c+1,5):
        diff=[]
        index=[]
        for i in range(len(dic_chips['chip_%s'%(c)])): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
            #dist=distance.cdist(dic_chips['chip_2'][i:i+1,6:8],dic_chips['chip_1'][:,6:8], 'euclidean')
            dist=distance.cdist(dic_chips['chip_%s'%(c)][i:i+1,6:8],dic_chips['chip_%s'%(j)][:,6:8], 'euclidean')
            d=np.where(dist<distancia)
            if len(d[1])==1:
                index.append(int(d[1][[np.argmin(dist[d])]]))
                diff.append((dic_chips['chip_%s'%(c)][i],dic_chips['chip_%s'%(j)][d[1][np.argmin(dist[d])]]))
            elif len(d[1])>1:
                index.append(int(d[1][[np.argmin(dist[d])]]))
                index.append(int(d[1][[np.argmax(dist[d])]]))
                diff.append((dic_chips['chip_%s'%(c)][i],dic_chips['chip_%s'%(j)][d[1][np.argmin(dist[d])]]))
        dic_chips['chip_%s'%(j)]=np.delete(dic_chips['chip_%s'%(j)],index,axis=0)
        print('Comunes y eliminadas de chip %s y chip %s ---->'%(b,j),len(diff),len(index))
    c+=1


# In[5]:


#make plot to check eliminated stars
fig,ax=plt.subplots(1,1,figsize=(20,10))
for chip in range(1,5):
    ax.scatter(dic_chips['chip_%s'%(chip)][:,6],dic_chips['chip_%s'%(chip)][:,7],alpha=0.2)
    


# In[6]:


#load field18 and add XY coordinates using astrometric information from field18_c3.fits (this images is from ESO fase 3)
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


# In[7]:


#bric=[]
#brick=np.concatenate((dic_chips['chip_1'],dic_chips['chip_2'],dic_chips['chip_3'],dic_chips['chip_4']),axis=0)
#fig,ax=plt.subplots(1,1,figsize=(20,10))
#ax.scatter(brick[:,6],brick[:,7],alpha=0.2)


# In[8]:


brick_all=[]
brick_all=np.concatenate((dic_chips['chip_1'],dic_chips['chip_2'],dic_chips['chip_3'],dic_chips['chip_4']),axis=0)
per=10# choose values with uncertainty in position smoller than 0.7 * mean(dx)
#valid=np.where((brick_all[:,8]<np.mean(brick_all[:,8])*per)&(brick_all[:,9]<np.mean(brick_all[:,9])*per))
valid=np.where((brick_all[:,8]<0.013)&(brick_all[:,9]<0.013))#0.013pixels cooresponding to a distant of 2mas
brick=brick_all[valid]
print(len(brick),len(brick_all))
##########################################################
lon,lat=w.pixel_to_world_values(brick_all[:,6],brick_all[:,7])
brick_all[:,0]=lon
brick_all[:,1]=lat
np.savetxt(pruebas+'zoc_alig_GNS_%s.txt'%(band),brick_all)
##########################################################
#now we are looping with a degree 1 ,2,...,
distancia=1
for degree in range(1,3):
    for loop in range(5):
        print(degree)
        brick=brick_all[valid]
        print(len(brick))
        diff=[]
        for i in range(len(field18)): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
            dist=distance.cdist(field18[i:i+1,-2:],brick[:,6:8], 'euclidean')
            d=np.where(dist<distancia)
            if len(d[1])>0:
                diff.append((field18[i],brick[d[1][np.argmin(dist[d])]]))
        print('comunes listas %s y %s ----> '%('brick','GNS'),len(diff))
        #diff=[diff[j] for j in range(len(diff)) if abs(diff[j][0][brillo]-diff[j][1][2])<l_mag]
        #print(len(diff))
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
        #degree=1
        Kx,Ky=polywarp(x1,y1,x2,y2,degree=degree)
        #print(Kx[0,0])
        xi=np.zeros(len(brick_all))
        yi=np.zeros(len(brick_all))
        x=[]
        y=[]
        x=brick_all[:,6]
        y=brick_all[:,7]
        x=np.array(x)
        y=np.array(y)
        for k in range(degree+1):
            for m in range(degree+1):
                xi=xi+Kx[k,m]*x**k*y**m
                yi=yi+Ky[k,m]*x**k*y**m
        brick_all[:,6]=xi
        brick_all[:,7]=yi
##########################################################
distancia=.4
diff=[]
for i in range(len(field18)): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
    dist=distance.cdist(field18[i:i+1,-2:],brick_all[:,6:8], 'euclidean')
    d=np.where(dist<distancia)
    if len(d[1])>0:
        diff.append((field18[i],brick_all[d[1][np.argmin(dist[d])]]))
diff=np.array(diff)
l_mag=0.6

x_shift=[(diff[t][0][-2]-diff[t][1][6]) for t in range(len(diff)) if abs(diff[t][0][brillo]-diff[t][1][2])<l_mag]
y_shift=[(diff[t][0][-1]-diff[t][1][7]) for t in range(len(diff)) if abs(diff[t][0][brillo]-diff[t][1][2])<l_mag]
print(len(x_shift))

x_shift=np.array(x_shift)*0.106
y_shift=np.array(y_shift)*0.106
fig,ax= plt.subplots(1,2,figsize=(20,10))
ls=[x_shift,y_shift]
nam=['x diff (arcsec)','y diff (arcsec)']
#ax[0]=plt.suptitle('Sigma Threshold at reconstruct = %s, CHIP  %s'%(th,ch),fontsize=20)
#plt.clf()
for h in range(len(ls)):
    sig_h=sigma_clipped_stats(ls[h],sigma=2,maxiters=20,cenfunc='mean')
    ax[h].hist(ls[h], bins=10,alpha=0.7, rwidth=0.85,color='g')
    #ax[h].axvline(np.mean(ls[h]), color='r', linestyle='dashed', linewidth=3)
    ax[h].axvline(sig_h[0], color='r', linestyle='dashed', linewidth=3)
    ax[h].grid(axis='both', alpha=0.75)
    #ax[h].legend(['Chip%s: mean= %.4f, std=%.4f'%(1,np.mean(ls[h]),np.std(ls[h]))],fontsize=20,markerscale=0,shadow=True,loc=3,handlelength=0)
    ax[h].legend([' dis=%s, #%s, mean= %.4f, std=%.4f'%(distancia,len(x_shift),sig_h[0],sig_h[2])],fontsize=20,markerscale=0,shadow=True,loc=3,handlelength=-0.0)
    #ax[h].legend([' #stars %s'%(len(x_shift))],fontsize=20,markerscale=0,shadow=True,loc=2,handlelength=0)
    ax[h].set_xlabel(nam[h],fontsize=20)
    ax[h].set_ylabel('# stars',fontsize=20)
    ax[h].tick_params(axis='x', labelsize=20)
    ax[h].tick_params(axis='y', labelsize=20)
fig.text(0.5, 0, 'Difference in position for common stars to GNS, band %s. Degr=%s'%(band,degree),fontsize=20, ha='center')

# In[10]:


distancia=0.3
diff=[]
for i in range(len(field18)): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
    dist=distance.cdist(field18[i:i+1,-2:],brick_all[:,6:8], 'euclidean')
    d=np.where(dist<distancia)
    if len(d[1])>0:
        diff.append((field18[i],brick_all[d[1][np.argmin(dist[d])]]))
diff=np.array(diff)
l_mag=1

x_shift=[(diff[t][0][-2]-diff[t][1][6]) for t in range(len(diff)) if abs(diff[t][0][brillo]-diff[t][1][2])<l_mag]
y_shift=[(diff[t][0][-1]-diff[t][1][7]) for t in range(len(diff)) if abs(diff[t][0][brillo]-diff[t][1][2])<l_mag]
print(len(x_shift))

x_shift=np.array(x_shift)*0.106
y_shift=np.array(y_shift)*0.106
fig,ax= plt.subplots(1,2,figsize=(20,10))
ls=[x_shift,y_shift]
nam=['x diff (arcsec)','y diff (arcsec)']
#ax[0]=plt.suptitle('Sigma Threshold at reconstruct = %s, CHIP  %s'%(th,ch),fontsize=20)
#plt.clf()
for h in range(len(ls)):
    sig_h=sigma_clipped_stats(ls[h],sigma=2,maxiters=20,cenfunc='mean')
    ax[h].hist(ls[h], bins=10,alpha=0.7, rwidth=0.85,color='g')
    #ax[h].axvline(np.mean(ls[h]), color='r', linestyle='dashed', linewidth=3)
    ax[h].axvline(sig_h[0], color='r', linestyle='dashed', linewidth=3)
    ax[h].grid(axis='both', alpha=0.75)
    #ax[h].legend(['Chip%s: mean= %.4f, std=%.4f'%(1,np.mean(ls[h]),np.std(ls[h]))],fontsize=20,markerscale=0,shadow=True,loc=3,handlelength=0)
    ax[h].legend([' dis=%s,#%s, mean= %.4f, std=%.4f'%(distancia,len(x_shift),sig_h[0],sig_h[2])],fontsize=20,markerscale=0,shadow=True,loc=3,handlelength=-0.0)
    #ax[h].legend([' #stars %s'%(len(x_shift))],fontsize=20,markerscale=0,shadow=True,loc=2,handlelength=0)
    ax[h].set_xlabel(nam[h],fontsize=20)
    ax[h].set_ylabel('# stars',fontsize=20)
    ax[h].tick_params(axis='x', labelsize=20)
    ax[h].tick_params(axis='y', labelsize=20)
fig.text(0.5, 0, 'Difference in position for common stars to GNS, band %s. Degr=%s'%(band,degree),fontsize=20, ha='center')
