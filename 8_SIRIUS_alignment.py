#!/usr/bin/env python
# coding: utf-8

# In[1]:


from astropy.io import fits
import numpy as np
from scipy.spatial import distance
from astropy.wcs import WCS


# In[2]:


#Works well with chips 1 and 4. Does not work with chip 2 and 3. Don`t know why yet.

band='H'
exptime=10
folder='im_jitter_NOgains/'
red='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/06_Reduce/054_H/dit_10/im_jitter_NOgains/'
pruebas='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/pruebas/'
sirius='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/SIRIUS/'
tmp='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder+'tmp/'
results='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder+'/results/'
name='NPL_054'
scripts='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/scripts/'


# In[3]:


chip=1
magerr_si=0.02
a_si, dec_si, J_si,dJ_si,H_si,dH_si,K_si,dK_si =np.loadtxt(sirius+'the_brick.txt',unpack=True)
si=np.loadtxt(sirius+'the_brick.txt')

ra,dec,x_mean,dx,y_mean,dy,mag,dmag,l,b =np.loadtxt(results+name+'_chip%s.txt'%(chip),unpack=True)
ch=np.loadtxt(results+name+'_chip%s.txt'%(chip))

six=[]
siy=[]
f = fits.open(tmp+'wt_chip%s.fits'%(chip))
w = WCS(f[1].header)
six,siy=w.world_to_pixel_values(si[:,0],si[:,1])
si=np.c_[si,six,siy]

print(w)


# In[4]:


stream = open(scripts+'polywarp.py')
read_file = stream.read()
exec(read_file)


# In[5]:


if band =='H': 
    m_si = H_si
    dm_si = dH_si
a_si, dec_si, J_si,dJ_si,H_si,dH_si,K_si,dK_si =np.loadtxt(sirius+'the_brick.txt',unpack=True)
count=0
for i in range(len(m_si)):
    if m_si[i]==0.0:
        #print(i,(m_si[i].dtype))
        count+=1
print(count)


# In[6]:


a_si, dec_si, J_si,dJ_si,H_si,dH_si,K_si,dK_si =np.loadtxt(sirius+'the_brick.txt',unpack=True)

if band =='H': 
    m_si = H_si
    dm_si = dH_si

if band =='Ks': 
    
    m_si = K_si
    dm_si = dK_si
    
valid = np.where(m_si >0.0)
a_si = a_si[valid]
dec_si = dec_si[valid]
m_si = m_si[valid]
dm_si = dm_si[valid]


#good = np.where((m_si >11.0) & (dm_si < magerr_si))
good = np.where((m_si >8.00)&  (m_si <13.00)& (dm_si < magerr_si*10000) )
good_dxy=np.where((dx < 0.002)&(dy < 0.002))
len(si[good])


# In[7]:


#ds=2.944444e-5*3
ds=12 # compare the lists (SIRIUS with Zocallis data), them with the common candidtaes within a distance <ds apply a polywarp, the apply the solotion to the whole list
      # and then loop again, making the ds 1 pixel smaller with each step.
for l in range(1,10):
    print(ds)
    x_good=ch[good_dxy][:,2]
    y_good=ch[good_dxy][:,4]
    xy_good=np.array([x_good,y_good]).T

    
    distancia = ds #2.944444e-5=1pixel in degrees #it is the same compare within the coordinates xy,with the distance in pixel, that did it with the RA,DEC with the distance in degrees!!!!!
    diff=[]
    dic_listas={}
    for i in range(len(si[good])): #distance of common point beetween chips. 
                dist=distance.cdist(si[good][i:i+1,[8,9]],xy_good, 'euclidean')
                d=np.where(dist<distancia)
                if len(d[1])>0:
                    diff.append((si[good][i],xy_good[d[1][np.argmin(dist[d])]]))
    dic_listas['Sirius_chip'+str(chip)]=diff
    diff=np.array(diff)
    print('comunes listas %s y %s ----> '%('Sirius',chip),len(diff))

    x1=[]
    y1=[]
    x2=[]
    y2=[]
    for i in range(len(dic_listas['Sirius_chip'+str(chip)])):
        x1.append(dic_listas['Sirius_chip'+str(chip)][i][0][8])
        y1.append(dic_listas['Sirius_chip'+str(chip)][i][0][9])
        x2.append(dic_listas['Sirius_chip'+str(chip)][i][1][0])
        y2.append(dic_listas['Sirius_chip'+str(chip)][i][1][1])

    x1=np.array(x1)
    y1=np.array(y1)
    x2=np.array(x2)
    y2=np.array(y2)

    Kx=[]
    Ky=[]
    deg=1
    Kx,Ky=polywarp(x1,y1,x2,y2,degree=deg)

    print(Kx)
    print(Ky)

    xi=np.zeros(len(ch))
    yi=np.zeros(len(ch))
    x=[]
    y=[]
    x=ch[:,2]
    y=ch[:,4]
    x=np.array(x)
    y=np.array(y)
    #for k in range(degree+1):
     #   for m in range(degree+1):
      #      xi=xi+Kx[k,m]*(x**k)*(y**m)
       #     yi=yi+Ky[k,m]*(x**k)*(y**m)
    for k in range(deg+1):
        for m in range(deg+1):
            xi=xi+Kx[k,m]*x**k*y**m
            yi=yi+Ky[k,m]*x**k*y**m


    ch[:,2]=xi
    ch[:,4]=yi
    
    ds=ds-1

np.savetxt(pruebas+'xy1_alig_chip%s.txt'%(chip),ch)


# In[8]:


# this second degree alignment doesnt make it bette. Quite the oposite. Maybe it is not well done??
#For chips 1 and 4 the first degree alignmet is good enought though 
'''
ds=5
for l in range(1,10):
    print(ds)
    x_good=ch[good_dxy][:,2]
    y_good=ch[good_dxy][:,4]
    xy_good=np.array([x_good,y_good]).T

    
    distancia = ds #2.944444e-5=1pixel in degrees #it is the same compare within the coordinates xy,with the distance in pixel, that did it with the RA,DEC with the distance in degrees!!!!!
    diff=[]
    dic_listas={}
    for i in range(len(si[good])): #distance of common point beetween chips. 
                dist=distance.cdist(si[good][i:i+1,[8,9]],xy_good, 'euclidean')
                d=np.where(dist<distancia)
                if len(d[1])>0:
                    diff.append((si[good][i],xy_good[d[1][np.argmin(dist[d])]]))
    dic_listas['Sirius_chip'+str(chip)]=diff
    diff=np.array(diff)
    print('comunes listas %s y %s ----> '%('Sirius',1),len(diff))

    x1=[]
    y1=[]
    x2=[]
    y2=[]
    for i in range(len(dic_listas['Sirius_chip'+str(chip)])):
        x1.append(dic_listas['Sirius_chip'+str(chip)][i][0][8])
        y1.append(dic_listas['Sirius_chip'+str(chip)][i][0][9])
        x2.append(dic_listas['Sirius_chip'+str(chip)][i][1][0])
        y2.append(dic_listas['Sirius_chip'+str(chip)][i][1][1])

    x1=np.array(x1)
    y1=np.array(y1)
    x2=np.array(x2)
    y2=np.array(y2)

    Kx=[]
    Ky=[]
    deg=2
    Kx,Ky=polywarp(x1,y1,x2,y2,degree=deg)

    print(Kx)
    print(Ky)

    xi=np.zeros(len(ch))
    yi=np.zeros(len(ch))
    x=[]
    y=[]
    x=ch[:,2]
    y=ch[:,4]
    x=np.array(x)
    y=np.array(y)
    #for k in range(degree+1):
     #   for m in range(degree+1):
      #      xi=xi+Kx[k,m]*(x**k)*(y**m)
       #     yi=yi+Ky[k,m]*(x**k)*(y**m)
    for k in range(deg+1):
        for m in range(deg+1):
            xi=xi+Kx[k,m]*x**k*y**m
            yi=yi+Ky[k,m]*x**k*y**m


    ch[:,2]=xi
    ch[:,4]=yi
    
    #ds=ds-1

np.savetxt(pruebas+'xy2_alig_chip%s.txt'%(chip),ch)
'''


# In[ ]:





# In[ ]:




