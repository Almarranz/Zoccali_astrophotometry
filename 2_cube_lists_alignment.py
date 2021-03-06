#!/usr/bin/env python
# coding: utf-8

# In[1]:



import numpy as np
from scipy.spatial import distance



# In[2]:


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
scripts='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/scripts/'


# In[3]:


stream = open(scripts+'polywarp.py')
read_file = stream.read()
exec(read_file)


# In[4]:


m=0
#m=0
for chip in range(1,5):
    n=np.loadtxt(jitter+'xy_off_xy_alig_chip'+str(chip)+'.txt')
    lim_x=[401,2150]
    lim_y=[401,2150]
    dic_stars={}
    dic_whole={}
    for i in range(1,len(n)+1):
        lista=np.loadtxt(tmp+'cube_stars_im'+str(i)+'_chip'+str(chip)+'.txt')
        #print('lista len:%s Chip:%s im:%s'%(len(lista),chip,i))
        dic_whole['cube_stars_im'+str(i)+'_chip'+str(chip)]=lista
        list_E=[lista[r] for r in range(len(lista)) if lim_x[0]-m<lista[r,0]<lim_x[1]+m and lim_y[0]+m<lista[r,1]<=lim_y[1]-m ]
        list_E=np.array(list_E)
        dic_stars['listE_im'+str(i)]=list_E
    dic_listas={}
    for l in range(1,len(dic_stars)+1):
    #for l in range(1,4):    
        distancia=2
    #for j in range(2,3):
        diff=[]
        im_a=1
        im_b=l
        for i in range(dic_stars['listE_im'+str(im_a)].shape[0]): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la m??s perque??a
                    dist=distance.cdist(dic_stars['listE_im'+str(im_a)][i:i+1,0:2],dic_stars['listE_im'+str(im_b)][:,0:2], 'euclidean')
                    d=np.where(dist<distancia)
                    if len(d[1])>0:
                        diff.append((dic_stars['listE_im'+str(im_a)][i],dic_stars['listE_im'+str(im_b)][d[1][np.argmin(dist[d])]]))
        dic_listas['stars_1'+str(l)]=diff
        print('comunes listas %s y %s ----> '%(im_a,im_b),len(diff))
        x1=[]
        y1=[]
        x2=[]
        y2=[]
        for i in range(len(dic_listas['stars_1'+str(l)])):
            x1.append(dic_listas['stars_1'+str(l)][i][0][0])
            y1.append(dic_listas['stars_1'+str(l)][i][0][1])
            x2.append(dic_listas['stars_1'+str(l)][i][1][0])
            y2.append(dic_listas['stars_1'+str(l)][i][1][1])

        x1=np.array(x1)
        y1=np.array(y1)
        x2=np.array(x2)
        y2=np.array(y2)

        Kx=[]
        Ky=[]
        degree=1
        Kx,Ky=polywarp(x1,y1,x2,y2,degree=degree)
        #print(Kx[0,0])
        xi=np.zeros(len(dic_whole['cube_stars_im'+str(l)+'_chip'+str(chip)]))
        yi=np.zeros(len(dic_whole['cube_stars_im'+str(l)+'_chip'+str(chip)]))
        x=[]
        y=[]
        x=dic_whole['cube_stars_im'+str(l)+'_chip'+str(chip)][:,0]
        y=dic_whole['cube_stars_im'+str(l)+'_chip'+str(chip)][:,1]
        x=np.array(x)
        y=np.array(y)
        #for k in range(degree+1):
         #   for m in range(degree+1):
          #      xi=xi+Kx[k,m]*(x**k)*(y**m)
           #     yi=yi+Ky[k,m]*(x**k)*(y**m)
        for k in range(degree+1):
            for m in range(degree+1):
                xi=xi+Kx[k,m]*x**k*y**m
                yi=yi+Ky[k,m]*x**k*y**m
        im2=dic_whole['cube_stars_im'+str(l)+'_chip'+str(chip)]
        im2[:,0]=xi
        im2[:,1]=yi
        '''
        diff=[]
        im_a=1
        for i in range(dic_stars['listE_im'+str(im_a)].shape[0]): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la m??s perque??a
                    dist=distance.cdist(dic_stars['listE_im'+str(im_a)][i:i+1,0:2],im2[:,0:2], 'euclidean')
                    d=np.where(dist<distancia)
                    if len(d[1])>0:
                        #diff.append((dic_stars['stars_im'+str(im_a)][i],im2[d[1][np.argmin(dist[d])]]))
                        diff.append((im2[d[1][np.argmin(dist[d])]]))
        dic_listas['stars_1'+str(l)]=diff
        '''
        np.savetxt(tmp+'cube_im'+str(l)+'_chip'+str(chip)+'_aligned_py.txt',im2,fmt='%.5f')
        print('Saved aligned list of im%s'%(l))
    print('Done with chip %s'%(chip))
        #print('comunes listas %s y %s ----> '%(im_a,l),len(diff))

