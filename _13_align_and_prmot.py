#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  6 09:45:23 2021

@author: amartinez
"""

# In[1]:



import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import distance
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
from scipy.stats import gaussian_kde
import sys
import pint

# In[2]:


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
tmp='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder+'tmp_bs/'
sirius='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/SIRIUS/'
GNS_ori='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/field12/'
scripts='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/scripts/'

#More than 1 ref star??? 
# yes: more=1(two stars) or more=2(three stars)
# no: more=0
more=0

stream = open(scripts+'polywarp.py')
read_file = stream.read()
exec(read_file)

distancia=1
chip=3# in this case only chip 2 and chip 3 have common elements with the GNS on the brick
unc=10

#############################################
#Varibles for foreground or inplace stars
#Uses foreground stars for alignment GNS_campo =1
#Doesnt use foreground stars for alignment GNS_campo =0
#Uses all stars GNS_campo=2
GNS_campo=0
#Uses foreground stars for proper motions Zoc_campo =1
#Doesnt use foreground stars for proper motions Zoc_campo =0
#Uses all stars Vel_campo=2
Vel_campo=0
#############################################
# Sigma cliping by velocityies and maximun limit in xy uncertainty fof proper motion calculation
# it would make a graph of unc_xy or unc_v for the smaller value. If the smoller valueis <1 it would discard those stars with 
# bigger uncertainties  
s=10
unc_xy= 1 #uncertainty limit for the position vector
unc_v=2 # uncertainty limit for the velocity vector
unc=10 #uncertainty limit for position on GNS stars
unc_z=1 #uncertainty limit for position on GNS stars
# Chose one lists for aligment
lst=3
if lst ==3:
    field12=np.loadtxt(GNS_ori+'field12_on_brick_reduced.txt')
elif lst==2 :
    field12=np.loadtxt(GNS_ori+'field12_on_brick_accu.txt')
elif lst==1:
    field12=np.loadtxt(GNS_ori+'field12_on_brick.txt')
lst_save=[lst]   
np.savetxt(tmp+'lst_chip%s.txt'%(chip),lst_save,fmt='%.0f')

#eliminates stas with H-Ks<1.3
h_ks=field12[:,10]-field12[:,12]
field12=np.c_[field12,h_ks]
if GNS_campo==1:
    in_place=np.where((h_ks<1.3))
elif GNS_campo==0:
    in_place=np.where((h_ks>1.3))
elif GNS_campo==2:
    in_place=np.where((h_ks))
field12=field12[in_place]
field12[:,0]*=0.5
field12[:,2]*=0.5
field12[:,1]*=0.5
field12[:,3]*=0.5
np.savetxt(GNS_ori+'field12_no_foreground.txt',field12,header='x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK, H-Ks,')

# x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK,H_Ks=np.loadtxt(GNS+'field12_no_foreground.txt',unpack=True)
# sys.exit("STOP")
# Here we are to select GNS stars by their uncertainty.
dxy_gns=np.sqrt((field12[:,1]*0.106)**2+(field12[:,3]*0.106)**2)
# low_g=np.where(dxy_gns<unc)
low_g=np.where((field12[:,1]*0.106<unc) & (field12[:,3]*0.106<unc) )
field12=field12[low_g]
x_gns=field12[:,0]
y_gns=field12[:,2]

gns_xy=np.array([x_gns,y_gns]).T

for chip in range(chip,chip+1):
    #a ,d , m, dm, f, df,x,y,dx,dy=np.loadtxt(tmp+'stars_calibrated_'+band+'_chip'+strn(chip)+'_sirius.txt',unpack=True)
    brick=np.loadtxt(tmp+'stars_calibrated_'+band+'_chip'+str(chip)+'_sirius.txt')
    # Here we are to select ZOC stars by their uncertainty.
    # dxy_zoc=np.sqrt(brick[:,8]**2+brick[:,9]**2)
    # low_z=np.where(dxy_zoc<unc)
    low_z=np.where((brick[:,8]<unc_z)&(brick[:,9]<unc_z))
    brick=brick[low_z]
    #We have to add the coordinates offset between the two lists
    #I this case I have choose this bu coomparing in Aladin
    if band=='H' and chip==3:
        ###########################################################################################################################
        #### Important note: Check if the list from where you are taking the GNS stars have the coordinates rescaled by 0.5 factor. 
        ###########################################################################################################################
        if more==0: #only 1 reference star





            xm_ref,ym_ref= 1070.18   ,1768.7# xm_ref is GNS
            xm_ref=xm_ref*0.5
            ym_ref=ym_ref*0.5
            xm,    ym    =1.174991699200000085e+03, 6.671052856000000020e+02

            xoff = xm_ref - xm
            yoff = ym_ref - ym
        if more==1: # 2 reference stars
            xm_ref1,ym_ref1=  1070.18,   1768.7# xm_ref is GNS
            xm_ref1=xm_ref1*0.5
            ym_ref1=ym_ref1*0.5
            xm1,    ym1    =1.174991699200000085e+03, 6.671052856000000020e+02
            xoff1 = xm_ref1 - xm1
            yoff1 = ym_ref1 - ym1
            
            xm_ref2,ym_ref2=   1206.31, 2003.85# xm_ref is GNS
            xm_ref2=xm_ref2*0.5
            ym_ref2=ym_ref2*0.5
            xm2,    ym2    =1.244916381799999954e+03,7.833552246000000423e+02

            xoff2 = xm_ref2 - xm2
            yoff2 = ym_ref2 - ym2
        
            xoff=np.mean([xoff1,xoff2])
            yoff=np.mean([yoff1,yoff2])
        if more==2:
 
            xm_ref1,ym_ref1=  1070.18,   1768.7# xm_ref is GNS
            xm_ref1=xm_ref1*0.5
            ym_ref1=ym_ref1*0.5
            xm1,    ym1    =1.174991699200000085e+03, 6.671052856000000020e+02
            xoff1 = xm_ref1 - xm1
            yoff1 = ym_ref1 - ym1
            
            xm_ref2,ym_ref2=   1206.31, 2003.85# xm_ref is GNS
            xm_ref2=xm_ref2*0.5
            ym_ref2=ym_ref2*0.5
            xm2,    ym2    =1.244916381799999954e+03,7.833552246000000423e+02

            xoff2 = xm_ref2 - xm2
            yoff2 = ym_ref2 - ym2
      



            xm_ref3,ym_ref3=   1564.25,1768.78# xm_ref is GNS
            xm_ref3=xm_ref3*0.5
            ym_ref3=ym_ref3*0.5
            xm3,    ym3    =1.421817382799999905e+03 , 6.629620972000000165e+02
            xoff3 = xm_ref3 - xm3
            yoff3 = ym_ref3 - ym3
        
            xoff=np.mean([xoff1,xoff2,xoff3])
            yoff=np.mean([yoff1,yoff2,yoff3])
        ##########################################################
    elif band=='H' and chip==2:
        ##########################################################
        if more ==0:
            xm_ref,ym_ref=694.834*0.5 , 569.783*0.5 # xm_ref is GNS
            xm,ym        = 9.832876586999999518e+02,2.260803955100000167e+03
    
            xoff = xm_ref - xm
            yoff = ym_ref - ym
        ##########################################################
        if more==1:
           
            xm_ref1,ym_ref1=694.834, 569.783 # xm_ref is GNS
            xm1,ym1        = 9.832876586999999518e+02,2.260803955100000167e+03
            xm_ref1=xm_ref1*0.5
            ym_ref1=ym_ref1*0.5
            xoff1 = xm_ref1 - xm1
            yoff1 = ym_ref1 - ym1
               


            xm_ref2,ym_ref2=   1444.05 , 747.72# xm_ref is GNS
            xm_ref2=xm_ref2*0.5
            ym_ref2=ym_ref2*0.5
            xm2,    ym2    =1.358597778300000073e+03,2.344092285199999878e+03

            xoff2 = xm_ref2 - xm2
            yoff2 = ym_ref2 - ym2
        
            xoff=np.mean([xoff1,xoff2])
            yoff=np.mean([yoff1,yoff2])
            
            
    
    print(' #'*20,'\n','xoff=%s yoff=%s'%(xoff,yoff),'\n'+' #'*20)    
    brick[:,6]+=xoff
    brick[:,7]+=yoff
    
    brick=np.array(brick)
    # In[4]:
    
    
    
    per=10# choose values with uncertainty in position smoller than 0.7 * mean(dx)
    #valid=np.where((brick_all[:,8]<np.mean(brick_all[:,8])*per)&(brick_all[:,9]<np.mean(brick_all[:,9])*per))
    #valid=np.where((brick[:,8]<0.013)&(brick[:,9]<0.013))#0.013pixels cooresponding to a distant of 2mas
    #brick=brick_all[valid]
    #print(len(brick),len(brick_all))
    ##########################################################
    
    ##########################################################
    #now we are looping with a degree 1 ,2,...,
   
    ciclo=10
    for degree in range(1,4):#Using degree 3 polynomial doesnt seem to improve things
        for loop in range(1,ciclo+1):
            print('Degree %s,iteration %s'%(degree,loop))
            diff=[]
            for i in range(len(gns_xy)): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
                dist=distance.cdist(gns_xy[i:i+1,0:2],brick[:,6:8], 'euclidean')
                d=np.where(dist<distancia)
                if len(d[1])>0:
                    diff.append((field12[i],brick[d[1][np.argmin(dist[d])]]))
            print('comunes listas %s y %s ----> '%('brick','GNS'),len(diff))
            #diff=[diff[j] for j in range(len(diff)) if abs(diff[j][0][brillo]-diff[j][1][2])<l_mag]
            #print(len(diff))
            diff=np.array(diff)
            x1=[]
            x2=[]
            y1=[]
            y2=[]
            for j in range(len(diff)):
                x1.append(diff[j][0][0])
                y1.append(diff[j][0][2])
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
        ciclo+=5
    ##########################################################
#     gns_txt=[diff[i][0][0:] for i in range(len(diff))]
#     zoc_txt=[diff[i][1][0:] for i in range(len(diff))]
#     x_dis=[diff[i][1][6]-diff[i][0][0] for i in range(len(diff))]
#     y_dis=[diff[i][1][7]-diff[i][0][2] for i in range(len(diff))]
#     displa=np.array([x_dis,y_dis]).T
#     zoc_txt=np.c_[zoc_txt,displa]
#     np.savetxt(tmp+'GNS_commons_w_Zoc_c%s.txt'%(chip),gns_txt,header='x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK, H-Ks')
#     np.savetxt(tmp+'Zoc_c%s_commons_w_GNS.txt'%(chip),zoc_txt,header='a ,d , m, dm, f, df,x,y,dx,dy,x_dis,y_dis. X and Y are the correspondig coorinates wit GNS, They are not the original ones!!!!')
#     np.savetxt(tmp+'dis_xy_chip%s.txt'%(chip),displa,header='Displacement in pixels.')


############ Histogram of common with GNS after aligment of Zocallis #################
diff=[]
if lst ==3:
    field12=np.loadtxt(GNS_ori+'field12_on_brick_reduced.txt')
elif lst==2 :
    field12=np.loadtxt(GNS_ori+'field12_on_brick_accu.txt')
elif lst==1:
    field12=np.loadtxt(GNS_ori+'field12_on_brick.txt')
field12[:,0]*=0.5
field12[:,2]*=0.5
field12[:,1]*=0.5
field12[:,3]*=0.5

h_ks=field12[:,10]-field12[:,12]
field12=np.c_[field12,h_ks]

x_gns=field12[:,0]
y_gns=field12[:,2]
gns_xy=np.array([x_gns,y_gns]).T

for i in range(len(gns_xy)): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
    dist=distance.cdist(gns_xy[i:i+1,0:2],brick[:,6:8], 'euclidean')
    d=np.where(dist<distancia)
    if len(d[1])>0:
        diff.append((field12[i],brick[d[1][np.argmin(dist[d])]]))
diff=np.array(diff)
l_mag=1

x_shift=[(diff[t][1][6]-diff[t][0][0]) for t in range(len(diff))]# if abs(diff[t][0][brillo]-diff[t][1][2])<l_mag]
y_shift=[(diff[t][1][7]-diff[t][0][2]) for t in range(len(diff))]# if abs(diff[t][0][brillo]-diff[t][1][2])<l_mag]
print(len(x_shift))

x_shift=np.array(x_shift)*0.106
y_shift=np.array(y_shift)*0.106
fig,ax= plt.subplots(1,2,figsize=(20,10))
ls=[x_shift,y_shift]
nam=['x diff (arcsec)','y diff (arcsec)']
#ax[0]=plt.suptitle('Sigma Threshold at reconstruct = %s, CHIP  %s'%(th,ch),fontsize=20)
#plt.clf()
for h in range(len(ls)):
    sig_h=sigma_clipped_stats(ls[h],sigma=5,maxiters=20,cenfunc='mean')
    ax[h].hist(ls[h], bins=10,alpha=0.7, rwidth=1,color='g',edgecolor='black',linewidth=2)
    #ax[h].axvline(np.mean(ls[h]), color='r', linestyle='dashed', linewidth=3)
    ax[h].axvline(sig_h[0], color='r', linestyle='dashed', linewidth=3)
    ax[h].grid(axis='both', alpha=0.75)
    #ax[h].legend(['Chip%s: mean= %.4f, std=%.4f'%(1,np.mean(ls[h]),np.std(ls[h]))],fontsize=20,markerscale=0,shadow=True,loc=3,handlelength=0)
    ax[h].legend(['Chip=%s,dmax=%spix, #%s, mean= %.4f, std=%.4f'%(chip,distancia,len(x_shift),sig_h[0],sig_h[2])],fontsize=15,markerscale=0,shadow=True,loc=3,handlelength=-0.0)
    #ax[h].legend([' #stars %s'%(len(x_shift))],fontsize=20,markerscale=0,shadow=True,loc=2,handlelength=0)
    ax[h].set_xlabel(nam[h],fontsize=20)
    ax[h].set_ylabel('# stars',fontsize=20)
    ax[h].tick_params(axis='x', labelsize=20)
    ax[h].tick_params(axis='y', labelsize=20)
fig.text(0.5, 0, 'Difference in position for common stars to GNS, band %s. Degr=%s'%(band,degree),fontsize=20, ha='center')
##########################################################

#%%
######### Savig common elements of both lits , plus a column for displacement on Zoc´s list
gns_txt=[diff[i][0][0:] for i in range(len(diff))]
zoc_txt=[diff[i][1][0:] for i in range(len(diff))]
x_dis=[diff[i][1][6]-diff[i][0][0] for i in range(len(diff))]
y_dis=[diff[i][1][7]-diff[i][0][2] for i in range(len(diff))]
displa=np.array([x_dis,y_dis]).T
zoc_txt=np.c_[zoc_txt,displa]
np.savetxt(tmp+'GNS_commons_w_Zoc_c%s.txt'%(chip),gns_txt,header='x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK, H-Ks')
np.savetxt(tmp+'Zoc_c%s_commons_w_GNS.txt'%(chip),zoc_txt,header='a ,d , m, dm, f, df,x,y,dx,dy,x_dis,y_dis. X and Y are the correspondig coorinates wit GNS, They are not the original ones!!!!')
np.savetxt(tmp+'dis_xy_chip%s.txt'%(chip),displa,header='Displacement in pixels.')

#%%

###### from here on proper motions ###########

ureg = pint.UnitRegistry()# to give units to values
# 'x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK,H-Ks'
GNS=np.loadtxt(tmp+'GNS_commons_w_Zoc_c%s.txt'%(chip))

# Whether we what to use foreground, background or all stars for the pm analis.
if Vel_campo==1:
    valid=np.where(GNS[:,14]<1.3)
elif Vel_campo==0:
    valid=np.where(GNS[:,14]>1.3)
elif Vel_campo==2:
    valid=np.where(GNS[:,14])


GNS=GNS[valid]
### 'a ,d , m, dm, f, df,x,y,dx,dy,x_displacement,y_displacement.
stars=np.loadtxt(tmp+'Zoc_c%s_commons_w_GNS.txt'%(chip))
stars=stars[valid]


low_g=np.where((GNS[:,1]*0.106<unc) & (GNS[:,3]*0.106<unc) )
GNS=GNS[low_g]
stars=stars[low_g]

dist=8*ureg.kpc




dx_dis=np.sqrt((GNS[:,1]*0.106)**2+stars[:,8]**2)
dy_dis=np.sqrt((GNS[:,3]*0.106)**2+stars[:,9]**2)
dxy=np.sqrt(dx_dis**2+dy_dis**2)
stars=np.c_[stars,dxy]

### 'a ,d , m, dm, f, df,x,y,dx,dy,x_displacement,y_displacement,dxy



vx=stars[:,10]/4*0.106
vy=stars[:,11]/4*0.106
dvel=np.sqrt(((vx*dx_dis/4)**2+(vy*dy_dis/4)**2)/(vx**2+vy**2))

vel=np.sqrt((stars[:,10]*0.106)**2+(stars[:,11]*0.106)**2)
stars=np.c_[stars,vel,dvel]
# sys.exit("STOP")
vel_clip=sigma_clip(vel,sigma=s)
stars=stars[vel_clip.mask==False]
#Aslo clipping GNS stars
GNS=GNS[vel_clip.mask==False]
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
ls_tosave=np.array(ls).T
#%%
##############################################################################
#save to try MCMC fitting
# save list with velocities in arsec/year

drap_x=np.sqrt(stars[:,8]**2+(GNS[:,1]*0.106)**2)/4
drap_y=np.sqrt(stars[:,9]**2+(GNS[:,3]*0.106)**2)/4
rap=np.c_[stars[:,10]*0.106/4,stars[:,11]*0.106/4,drap_x,drap_y]*1000
np.savetxt('/Users/amartinez/Desktop/PhD/python/Gaussian_fit/'+'arcsec_vx_vy_chip%s.txt'%(chip),rap,header='vx,vy,dvx,dvy (in arcsec/yr)')

#trnasfor dx and dy form GNS and Zoc to km/s and compute dvx and dvy separetly
ls_tosave=np.array(ls).T
GNS=GNS[low_xy]

dx1_tosave=(GNS[:,1]*0.106)
dx2_tosave=stars[:,8]

dx1_tosave=dx1_tosave*ureg.arcsec
dx2_tosave=dx2_tosave*ureg.arcsec

dx1_tosave=dx1_tosave.to('rad')*dist
dx2_tosave=dx2_tosave.to('rad')*dist
dvx_tosave=np.sqrt((dx1_tosave/(4*365*24*3600))**2+(dx2_tosave/(4*365*24*3600))**2)


dy1_tosave=(GNS[:,3]*0.106)*ureg.arcsec
dy2_tosave=stars[:,9]*ureg.arcsec

dy1_tosave=dy1_tosave.to('rad')*dist
dy2_tosave=dy2_tosave.to('rad')*dist
dvy_tosave=np.sqrt((dy1_tosave/(4*365*24*3600))**2+(dy2_tosave/(4*365*24*3600))**2)

# ls_tosave=np.c_[ls_tosave,dvx_tosave,dvy_tosave]
np.savetxt('/Users/amartinez/Desktop/PhD/python/Gaussian_fit/'+'vx_vy_chip%s.txt'%(chip),ls_tosave,header='vx,vy,dvx,dvy (in km/s)')


##############################################################################
#%%
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
    if s<8:
        ax[h].text(his[1][0],max(his[0]/2),'%s$\sigma$ clipped velocity'%(s),color='k',fontsize=20,zorder=3,weight='bold') 
    # ax[h].text(his[1][0],max(his[0]/2-his[0]/10),r'$\sigma_{\vec {v}}$(ars/yr)<%s'%(unc/4),color='k',fontsize=20,zorder=3,weight='bold') 
    if unc_v<1 or unc_xy<1:
        if unc_v<unc_xy:
            ax[h].text(his[1][0],max(his[0]/2+1-his[0]/10),r'$\sigma_{\vec {v}}$<%s(ars/yr)'%(unc_v),color='k',fontsize=20,zorder=3,weight='bold') 
        else:
            ax[h].text(his[1][0],max(his[0]/2+1-his[0]/10),r'$\sigma_{\vec {xy}}$<%s"'%(unc_xy),color='k',fontsize=20,zorder=3,weight='bold') 
    if unc<1:
         ax[h].text(his[1][0],max(his[0]/2-his[0]/10),r'GNS($\sigma_{x},\sigma_{y}$)<%s"'%(unc),color='k',fontsize=20,zorder=3,weight='bold')
    if unc_z<1:
         ax[h].text(his[1][0],max(his[0]/2-his[0]/10*1.5),r'Zoc($\sigma_{x},\sigma_{y}$)<%s"'%(unc_z),color='k',fontsize=20,zorder=3,weight='bold')

if GNS_campo==1:
    if Vel_campo==0:
        fig.suptitle('Foreground stars in alignment and NO Foreground stars in proper motions',fontsize=20)
    elif Vel_campo==1:
        fig.suptitle('Foreground stars in alignment and Foreground stars in proper motions',fontsize=20)
    elif Vel_campo==2:
        fig.suptitle('Foreground stars in alignment ALL stars in proper motions',fontsize=20)
elif GNS_campo==0:
    if Vel_campo==0:
        fig.suptitle('NO Foreground stars in alignment and NO Foreground stars in proper motions',fontsize=20)
    elif Vel_campo==1:
        fig.suptitle('NO Foreground stars in alignment and Foreground stars in proper motions',fontsize=20)
    elif Vel_campo==2:
        fig.suptitle('NO Foreground stars in alignment ALL stars in proper motions',fontsize=20)
elif GNS_campo==2:
    if Vel_campo==0:
        fig.suptitle('ALL stars in alignment and NO Foreground stars in proper motions',fontsize=20)
    elif Vel_campo==1:
        fig.suptitle('ALL stars in alignment and Foreground stars in proper motions',fontsize=20)
    elif Vel_campo==2:
        fig.suptitle('ALL stars in alignment ALL stars in proper motions',fontsize=20)


#fig.text(0.5, 0, 'Difference in position for common stars to GNS, band %s.'%(band),fontsize=20, ha='center'1.4

#%% 
######### Uncertainty in position or velocities plots ###############
# if unc_v>unc_xy:
#     #### Uncertainty in POSITION
#     ### 'a ,d , m, dm, f, df,x,y,dx,dy,x_displacement,y_displacement.
#     stars_all=np.loadtxt(tmp+'Zoc_c%s_commons_w_GNS.txt'%(chip))
#     # ### 'x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK,H-Ks'
#     GNS_all=np.loadtxt(tmp+'GNS_commons_w_Zoc_c%s.txt'%(chip))
    
#     dx_dis_all=np.sqrt((GNS_all[:,1]*0.106)**2+stars_all[:,8]**2)
#     dy_dis_all=np.sqrt((GNS_all[:,3]*0.106)**2+stars_all[:,9]**2)
#     dxy_all=np.sqrt(dx_dis_all**2+dy_dis_all**2)
    
#     menor=0
#     for i in range(len(dxy)):
#         if dxy[i]<unc_xy:
#             menor+=1
#     # stars=np.c_[stars,dxy]
#     ### 'a ,d , m, dm, f, df,x,y,dx,dy,x_displacement,y_displacement.
#     fig, ax=plt.subplots(1,1,figsize=(10,10))
#     ax.scatter(stars_all[:,2],dxy_all,color='k',alpha=0.3)
#     ax.set_xlabel('[H]',fontsize=20)
#     ax.set_ylabel(r'$\sigma_{\vec {xy}}$(arcsec)',fontsize=20)
#     # ax.set_ylim(0,0.03)
#     # ax.axhline(l_min, color='r', linestyle='dashed', linewidth=3)
#     # ax.axhline(l_max, color='r', linestyle='dashed', linewidth=3)
#     if unc_xy <1:
#         ax.axhline(unc_xy, color='g', linestyle='dashed', linewidth=3,zorder=3)
#         ax.text(min(stars[:,2]),unc_xy+unc_xy/10,'#stars= %s'%(len(stars[:,2])), color='g',fontsize=20,weight='bold',zorder=3)
#     ax.tick_params(axis='x', labelsize=20)
#     ax.tick_params(axis='y', labelsize=20)
    

# else:
#     #### Uncertainty in VELOCITY
  
#     ### 'a ,d , m, dm, f, df,x,y,dx,dy,x_displacement,y_displacement.
#     stars_all=np.loadtxt(tmp+'Zoc_c%s_commons_w_GNS.txt'%(chip))
#     # ### 'x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK,H-Ks'
#     GNS_all=np.loadtxt(tmp+'GNS_commons_w_Zoc_c%s.txt'%(chip))
    
#     dx_dis_all=np.sqrt((GNS_all[:,1]*0.106)**2+stars_all[:,8]**2)
#     dy_dis_all=np.sqrt((GNS_all[:,3]*0.106)**2+stars_all[:,9]**2)
#     dxy_all=np.sqrt(dx_dis_all**2+dy_dis_all**2)
    
#     vel_all=np.sqrt((stars_all[:,10]*0.106)**2+(stars_all[:,11]*0.106)**2)
#     vx_all=stars_all[:,10]/4*0.106
#     vy_all=stars_all[:,11]/4*0.106
    
#     dvel_all=np.sqrt(((vx_all*dx_dis_all/4)**2+(vy_all*dy_dis_all/4)**2)/(vx_all**2+vy_all**2))
#     menor=0
#     for i in range(len(dvel)):
#         if dvel[i]<unc_v/4:
#             menor+=1
#     # stars=np.c_[stars,dxy]
#     ### 'a ,d , m, dm, f, df,x,y,dx,dy,x_displacement,y_displacement.
#     fig, ax=plt.subplots(1,1,figsize=(10,10))
#     ax.scatter(stars_all[:,2],dvel_all,color='k',alpha=0.3)
#     ax.set_xlabel('[H]',fontsize=20)
#     ax.set_ylabel(r'$\sigma_{\vec {v}}$(arcsec/yr)',fontsize=20)
#     # ax.set_ylim(0,0.03)
#     # ax.axhline(l_min, color='r', linestyle='dashed', linewidth=3)
#     # ax.axhline(l_max, color='r', linestyle='dashed', linewidth=3)
#     if unc_v <1:
#         ax.axhline(unc_v, color='g', linestyle='dashed', linewidth=3,zorder=3)
#         ax.text(min(stars[:,2]),unc_v+unc_v/40,'#stars= %s'%(len(stars[:,2])), color='g',fontsize=20,weight='bold',zorder=3)
#     ax.tick_params(axis='x', labelsize=20)
#     ax.tick_params(axis='y', labelsize=20)
    
    
    
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





























