#!/usr/bin/env python
# coding: utf-8

# In[1]:



import numpy as np



# In[2]:


band='H'
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


# In[9]:


m=0
#m=0
lim_x=[401,2150]
lim_y=[401,2150]
#cut the corresponding area of the images (acording to the wt)
#the commented lines (with '''....''') add galactic coordenates at the end of each list(it is time comsuming)
for chip in range(1,5):
    n=np.loadtxt(jitter+'xy_off_xy_alig_chip'+str(chip)+'.txt')
    '''
    fn = get_pkg_data_filename(tmp+'wt_chip1.fits')#, package='astropy.wcs.tests')
    f = fits.open(fn)
    w = WCS(f[1].header)
    '''
    dic_zp_ind={}
    list_ind=['list_A','list_B','list_C','list_D','list_E','list_F','list_G','list_H','list_I']
    for li in list_ind:
        dic_zp_ind[li]=[]
    for i in range(1,len(n)+1):
        '''
        coord=np.loadtxt(tmp+'cube_im'+str(i)+'_chip'+str(chip)+'_aligned_py.txt')
        stars_g=[w.pixel_to_world(coord[i][0], coord[i][1]).galactic for i in range(len(coord))]
        t=[]
        t = QTable([stars_g], names=['coord'])
        df=t.to_pandas()
        arr=df.to_numpy()
        l=arr[:,0]
        b=arr[:,1]
        '''
        lista=np.loadtxt(tmp+'cube_im'+str(i)+'_chip'+str(chip)+'_aligned_py.txt')
        '''
        lista=np.c_[lista,l,b]
        '''
        lista=np.loadtxt(tmp+'cube_im'+str(i)+'_chip'+str(chip)+'_aligned_py.txt')
        #print('lista len:%s Chip:%s im:%s'%(len(lista),chip,i))
        list_A=[lista[r] for r in range(len(lista)) if lista[r,0]<lim_x[0]+m and lista[r,1]<lim_y[0]-m]
        list_B=[lista[r] for r in range(len(lista)) if lista[r,0]<lim_x[0]+m and lim_y[0]+m<lista[r,1]<=lim_y[1]-m]
        list_C=[lista[r] for r in range(len(lista)) if lista[r,0]<lim_x[0]+m and lista[r,1]>lim_y[1]+m]
        
        list_D=[lista[r] for r in range(len(lista)) if lim_x[0]-m<lista[r,0]<lim_x[1]+m and lista[r,1]<lim_y[0]-m ]
        list_E=[lista[r] for r in range(len(lista)) if lim_x[0]-m<lista[r,0]<lim_x[1]+m and lim_y[0]+m<lista[r,1]<=lim_y[1]-m ]
        list_F=[lista[r] for r in range(len(lista)) if lim_x[0]-m<lista[r,0]<lim_x[1]+m and lista[r,1]>lim_y[1]+m]
        
        list_G=[lista[r] for r in range(len(lista)) if lista[r,0]>lim_x[1]-m and lista[r,1]<lim_y[0]-m ]
        list_H=[lista[r] for r in range(len(lista)) if lista[r,0]>lim_x[1]-m and lim_y[0]-m<lista[r,1]<lim_y[1]+m ]
        list_I=[lista[r] for r in range(len(lista)) if lista[r,0]>lim_x[1]-m and lista[r,1]>lim_y[1]+m]
        
        dic_lists={}
        list_good={}
        dic_lists={'list_A':list_A,'list_B':list_B,'list_C':list_C,'list_D':list_D,'list_E':list_E
              ,'list_F':list_F,'list_G':list_G,'list_H':list_H,'list_I':list_I}
        vocab = list( dic_lists.keys())
        
        list_good=sorted(dic_lists, key=lambda s: len(dic_lists[s]), reverse=True)
        print('######### im%s ######### '%(i))
        for u in range(4):
            print(list_good[u])
            np.savetxt(tmp+'cube_'+list_good[u]+'_im'+str(i)+'_chip'+str(chip)+'aligned.txt',dic_lists[list_good[u]], fmt='%.5f')
            dic_zp_ind[list_good[u]].append(i)
        #print('A %s,B %s,C %s,D %s,E %s,F %s,G %s,H %s,I %s'%(len(list_A),len(list_B),len(list_C),len(list_D),len(list_E),len(list_F),len(list_G),len(list_H),len(list_I)))
    for key in dic_zp_ind:
        np.savetxt(tmp+'zp_'+key+'_chip%s.txt'%(chip),dic_zp_ind[key],fmt='%i')
    print('Chip %s done'%(chip))  
    


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




