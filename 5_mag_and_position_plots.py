#!/usr/bin/env python
# coding: utf-8

# In[1]:
#This scripts plots the dmag vs mag and dx vs mag
#Generate a josn file with arrays of differentes legnts of the diferentes positions 
#and mag. There are different legnth couse not all the stars have the same numeber
#of mesaurements
#also generated mean fluxes a df

import glob
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import distance
from scipy.stats import gaussian_kde
import json


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


# In[3]:


for chip in range(1,5):
    if band=='H':
        ind=['A','B','C','D','E','F','G','H']#,'I'] # solo las listas con mas de 3 pointings
    elif band=='Ks':
        ind=['B','D','E','F','H']
    mag_med=[]
    error_mag_med=[]
    x_all=[]
    y_all=[]
    f_med=[]
    df_med=[]
    if folder=='im_sky_ESOReflex/':
        zp_list=np.loadtxt(results+"zp_%s_ESO_chip%s.txt"%(band),chip)
        print('Usando ESO')
    else:
        zp_list=np.loadtxt(results+"zp_%s_chip%s.txt"%(band,chip))
        print('Usando my sky, chip%s'%(chip))
    if band=='H':
        zp1=26.32
    else:
        zp1=25.63
    zp_im=np.insert(zp_list,0,zp1)#'np.insert()' could be useful to add a whole column of galactic coordinates
    len(zp_im)
    for ind in ind:
        name='cube_list_'+ind+'_im*_chip%saligned.txt'%(chip)
        zp_ind = np.loadtxt(tmp+'zp_list_%s_chip%s.txt'%(ind,chip)).astype(int)

        sub = sorted(glob.glob(tmp+ name), key=os.path.getmtime)
        #for i in sub:
            #print(i)
        dic_stars={}
        distancia=1
        s=0
        for i in sub:
            s+=1
            dic_stars['stars_im'+str(s)]=np.loadtxt(i)
        print(ind,len(dic_stars))
        diff=[]
        diff_Ax=[]
        diff_Ay=[]
        im_a=1
        im_b=2
        list_D=[]
        for i in range(dic_stars['stars_im'+str(im_a)].shape[0]): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
                    dist=distance.cdist(dic_stars['stars_im'+str(im_a)][i:i+1,0:2],dic_stars['stars_im'+str(im_b)][:,0:2], 'euclidean')
                    d=np.where(dist<distancia)
                    if len(d[1])>0:
                        diff.append((dic_stars['stars_im'+str(im_a)][i],dic_stars['stars_im'+str(im_b)][d[1][np.argmin(dist[d])]]))
                        diff_Ax.append((dic_stars['stars_im'+str(im_a)][i][0],dic_stars['stars_im'+str(im_b)][d[1][np.argmin(dist[d])]][0]))
                        diff_Ay.append((dic_stars['stars_im'+str(im_a)][i][1],dic_stars['stars_im'+str(im_b)][d[1][np.argmin(dist[d])]][1]))

        for i in range(len(diff)):
            mean_x=np.mean([diff[i][0][0],diff[i][1][0]])
            mean_y=np.mean([diff[i][0][1],diff[i][1][1]])
            list_D.append((mean_x,mean_y,diff[i][0][2],diff[i][1][2]))
        list_D=np.array(list_D)
        print('comunes a listas de 1 a %s'%(im_b),len(list_D))
        for j in range(3,len(zp_ind)+1):#loop comparando las comunes de 1,2 con las demas(123,1234,...,12...8)
            aux=[]
            aux_x=[]
            aux_y=[]
            for i in range(list_D.shape[0]): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
                        dist=distance.cdist(list_D[i:i+1,0:2],dic_stars['stars_im'+str(j)][:,0:2], 'euclidean')
                        d=np.where(dist<distancia)
                        if len(d[1])>0:
                            aux.append((list_D[i],dic_stars['stars_im'+str(j)][d[1][np.argmin(dist[d])]]))
                            aux_x.append((diff_Ax[i]+(dic_stars['stars_im'+str(j)][d[1][np.argmin(dist[d])]][0],)))
                            aux_y.append((diff_Ay[i]+(dic_stars['stars_im'+str(j)][d[1][np.argmin(dist[d])]][1],)))
            list_D=np.zeros(shape=(len(aux),j+2))
            diff_Ax=aux_x
            diff_Ay=aux_y
            for i in range(len(aux)):
                mean_x=np.mean([aux[i][0][0],aux[i][1][0]])
                mean_y=np.mean([aux[i][0][1],aux[i][1][1]])
                list_D[i][0:2]=mean_x,mean_y
                for k in range(2,j+1):
                    list_D[i][k]=aux[i][0][k]
                list_D[i][j+1]=aux[i][1][2]
            list_D=np.array(list_D)
            #np.savetxt(tmp+'x_band'+band+'_dit'+str(exptime)+'_chip'+str(chip)+'.txt',aux_x,fmt='%.5f')
            #np.savetxt(tmp+'y_band'+band+'_dit'+str(exptime)+'_chip'+str(chip)+'.txt',aux_y,fmt='%.5f')
            print('comunes a listas de 1 a %s'%(j),len(list_D))
        x_all.append(aux_x)
        y_all.append(aux_y)

        if band=='H':
            ZP=26.32
            c='k'
        else:
            ZP=25.63
            c='b'
        fluxes=[list_D[i][2:] for i in range(len(list_D))]
        for i in range(len(fluxes)):
            suma=[]
            suma_f=[]
            for j in range(len(fluxes[0])):
                #print(zp_im[zp_ind[j]-1])
                suma.append(zp_im[zp_ind[j]-1]-2.5*np.log10(fluxes[i][j]/exptime))
                suma_f.append(fluxes[i][j])
            mag_med.append((np.mean(suma)))
            error_mag_med.append((np.std(suma)/np.sqrt(len(fluxes[0]))))  
            f_med.append(np.mean(suma_f))
            df_med.append((np.std(suma_f)/np.sqrt(len(fluxes[0]))))
    fig,ax=plt.subplots(1,figsize=(10,10))
    #plt.gca().invert_yaxis()
    #plt.plot(mag,e_mag,'x',scalex=True,alpha=0.5,color='k')
    ax.plot(mag_med,error_mag_med,'.',scalex=True,alpha=1,color=c,markersize=5)
    ax.set_ylim(0,0.2)
    #plt.plot(x_peq,26.3 - 2.5 *np.log10(m_peq/1.26),'x',color='red')

    #plt.axhline(np.mean(e_mag), color='r', linestyle='dashed', linewidth=2)

    ax.set_xlabel('[%s]'%(band),fontsize=20)
    ax.set_ylabel('[d%s]'%(band),fontsize=20)
    ax.set_title('Chip=%s. DIT=%s, stars=%s'%(chip,exptime,len(mag_med)),fontsize=10)
    
    dx=[]
    dy=[]
    x_all=np.array(x_all)
    y_all=np.array(y_all)
    x_col=[]
    y_col=[]
    for i in range(len(x_all)):
        for j in range(len(x_all[i])):
            dx.append(np.std(x_all[i][j])/np.sqrt(len(x_all[i][j]))*0.106)
            dy.append(np.std(y_all[i][j])/np.sqrt(len(x_all[i][j]))*0.106)
            x_col.append((x_all[i][j]))
            y_col.append((y_all[i][j]))
    with open(tmp+'x_commons_chip'+str(chip)+'.json', 'w') as fx:
        json.dump(x_col,fx)
    with open(tmp+'y_commons_chip'+str(chip)+'.json', 'w') as fy:
        json.dump(y_col,fy)       
    fig,ax=plt.subplots(1,2,figsize=(20,10))
    np.savetxt(tmp+'mag_media_chip%s.txt'%(chip),mag_med)
    np.savetxt(tmp+'error_mag_med_chip%s.txt'%(chip),error_mag_med)
    np.savetxt(tmp+'f_chip%s.txt'%(chip),f_med,fmt='%.6e')
    np.savetxt(tmp+'df_med_chip%s.txt'%(chip),df_med,fmt='%.6e')
    
    ax[0].plot(mag_med,dx,'.',scalex=True,alpha=0.3,color=c,markersize=1)
    ax[0].set_xlabel('[%s]'%(band),fontsize=20)
    ax[0].grid()
    ax[0].set_ylabel('X std/sqrt(N) (arcsec)',fontsize=20)
    #ax[0].set_ylim(0,0.06)
    #ax[1].set_ylim(0,0.06)
    ax[1].set_xlabel('[%s]'%(band),fontsize=20)
    ax[1].set_ylabel('Y std/sqrt(N) (arcsec)',fontsize=20)
    ax[1].grid()
    #ax[0].axhline(np.mean(std_x), color='r',linestyle='dashed', linewidth=3)
    ax[1].plot(mag_med,dy,'.',scalex=True,alpha=0.3,color=c,markersize=1)
    #ax[1].axhline(np.mean(std_y), color='r',linestyle='dashed', linewidth=3)
    plt.suptitle('Total stars = %s.Chip=%s'%(len(dx),chip),fontsize=20)
    #if folder =='im_sky_ESOReflex/':
     #   ax[0].text(12,0.005,'ESO sky',fontsize=20)
      #  plt.savefig(results+'mag_vs_xy_ESO.png',overwrite=True)
    #else:

    if exptime==2:
        map_c='viridis'
    else:
         map_c='inferno'
    xy = np.vstack([mag_med,dx])
    xy1 = np.vstack([mag_med,dy])
    z = gaussian_kde(xy)(xy)
    z1= gaussian_kde(xy1)(xy1)

    fig, ax = plt.subplots(1,2,figsize=(15,10))
    #ax.invert_yaxis()
    #ax.set_ylim(16,10)
    ax[0].scatter(mag_med, dx, c=z, s=5,alpha=1,cmap=map_c)
    ax[0].set_xlabel('[%s]'%(band),fontsize=20)
    ax[0].grid()
    ax[0].set_ylim(0,0.02)
    ax[1].set_ylim(0,0.02)
    ax[0].set_ylabel('X std/sqrt(N) (arcsec)',fontsize=20)
    ax[1].scatter(mag_med, dy, c=z1, s=5,alpha=1,cmap=map_c)
    ax[1].set_xlabel('[%s]'%(band),fontsize=20)
    ax[1].set_ylabel('Y std/sqrt(N) (arcsec)',fontsize=20)
    ax[1].grid()
    plt.suptitle('Total stars = %s.Chip=%s'%(len(dx),chip),fontsize=20)
    #fig.colorbar(plot)
    #if folder =='im_sky_ESOReflex/':
     #   ax[0].text(12,0.005,'ESO sky',fontsize=20)
      #  plt.savefig(results+'mag_vs_xy_ESO2.png',overwrite=True)
    #else:
    #plt.savefig(results+'mag_vs_xy2_IDL.png',overwrite=True)
    #plt.show()


# In[ ]:

print(ind)



# In[ ]:




