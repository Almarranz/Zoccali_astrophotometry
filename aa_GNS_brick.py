#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 10 13:07:44 2021

@author: amartinez
"""

import numpy as np
import matplotlib.pyplot as plt
import astroalign as aa
from astropy.io.fits import getheader
from astropy.io import fits
from scipy.spatial import distance
import sys
folder='im_jitter_NOgains/'
# chip=4
zone='Z2'

#%%
#%%
in_brick=0
if in_brick==1:
    GNS='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/field12/'
    lst=1#1 is the biggest list and 3 the smallest one
    if lst==1:
        x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK=np.loadtxt(GNS+'field12_on_brick.txt',unpack=True)
    elif lst==2:
            x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK=np.loadtxt(GNS+'field12_on_brick_accu.txt',unpack=True)
    elif lst==3:
            x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK=np.loadtxt(GNS+'field12_on_brick_reduced.txt',unpack=True)

elif in_brick==0:
    field=3# field 16(lst 2 ,3) field3(lst 1,4)
    lst=4
    GNS='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/field%s/'%(field)
    # x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK=np.loadtxt(GNS+'cat_Ban_%s_%s.txt'%(field,lst),unpack=True)
    x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK=np.loadtxt(GNS+'%s_cat_Ban_%s_%s.txt'%(zone,field,lst),unpack=True)#'Z1_..' hace referencia a una lista sobrre una zona del mismo tamaño de Zone A sobre el brock
 
    print('som')
 #%%       

print(len(x_gns))
# gns=np.loadtxt(GNS+'field12_on_brick_accu.txt')
valid=np.where((mH<90) & (mK<90) )
x_gns=x_gns[valid]
y_gns=y_gns[valid]
mH=mH[valid]
mK=mK[valid]
H_Ks=np.where(mH-mK>1.3)
x_gns=x_gns[H_Ks]
y_gns=y_gns[H_Ks]

x_gns=x_gns*0.5
y_gns=y_gns*0.5
GNS_xy=np.array([x_gns,y_gns]).T
print('Elements in GNS: %s'%(len(GNS_xy)))

#%%
if in_brick==1:
    tmp='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_H/dit_10/'+folder+'tmp_bs/'
    a ,d , m, dm, f, df,x,y,dx,dy= np.loadtxt(tmp+'BRICK_stars_calibrated_H_chip'+str(chip)+'_sirius.txt',unpack=True)
    zoc= np.loadtxt(tmp+'BRICK_stars_calibrated_H_chip'+str(chip)+'_sirius.txt')
    zoca=np.array([x,y]).T
    print('Elements in Zoc: %s'%(len(x)))
elif in_brick==0:
    tmp='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/058_H/dit_10/'+folder+'tmp_bs/'
    # a ,d , m, dm, f, df,x,y,dx,dy= np.loadtxt(tmp+'stars_calibrated_H_on_field%s_%s.txt'%(field,lst),unpack=True)
    a ,d , m, dm, f, df,x,y,dx,dy= np.loadtxt(tmp+'%s_stars_calibrated_H_on_field%s_%s.txt'%(zone,field,lst),unpack=True)#'Z1_..' hace referencia a una lista sobrre una zona del mismo tamaño de Zone A sobre el Brick
 

    # zoc= np.loadtxt(tmp+'stars_calibrated_H_on_field%s_%s.txt'%(field,lst))
    zoc= np.loadtxt(tmp+'%s_stars_calibrated_H_on_field%s_%s.txt'%(zone,field,lst))
    zoca=np.array([x,y]).T
    print('Elements in Zoc: %s'%(len(x)))
#%%

 # -638.350, 215.718 manually stimated xoff yoff
check_x,check_y=2,2
while abs(check_x) >1 or abs(check_y)>1  :
    m,(_,_)= aa.find_transform(zoca, GNS_xy,max_control_points=250)
    print('For fiel%s, chip%s'%(field,lst)+'\n'+"Translation: (x, y) = (%.2f, %.2f)"%(m.translation[0],m.translation[1]))
    print("Rotation: %.3f degrees"%(m.rotation * 180.0 / np.pi))
    print("Scale factor: %.4f"%(m.scale))
    
    test_gns = aa.matrix_transform(zoca, m.params)
            
    print(20*'#'+'\n'+'CHECKING'+'\n'+20*'#')
    t,(_,_)= aa.find_transform(test_gns,GNS_xy,max_control_points=250)
    print("Translation: (x, y) = (%.2f, %.2f)"%(t.translation[0],t.translation[1]))
    print("Rotation: %.3f degrees"%(t.rotation * 180.0 / np.pi))
    print("Scale factor: %.4f"%(t.scale))
    check_x= t.translation[0]
    check_y= t.translation[1]
print('___NOW TRANSFORMING ZOC___')  
trans=aa.matrix_transform(zoca, m.params)
zoc[:,6]=trans[:,0]
zoc[:,7]=trans[:,1]
# =============================================================================
# if in_brick==1:
#     np.savetxt(tmp+'%s_aa_BRICK_stars_calibrated_H_chip'+str(chip)+'_sirius.txt'%(zone),zoc, fmt='%.5f') 
# elif in_brick==0:
#     np.savetxt(tmp + '%s_aa_stars_calibrated_H_on_field%s_%s.txt'%(zone,field,lst),zoc, fmt='%.5f')
# 
# 
# =============================================================================





