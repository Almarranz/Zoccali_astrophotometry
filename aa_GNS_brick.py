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
chip=3
GNS='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/field12/'
tmp='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_H/dit_10/'+folder+'tmp_bs/'
#%%
#%%
lst=1#1 is the biggest list and 3 the smallest one
if chip == 3:
    if lst==1:
        x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK=np.loadtxt(GNS+'field12_on_brick.txt',unpack=True)
    elif lst==2:
            x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK=np.loadtxt(GNS+'field12_on_brick_accu.txt',unpack=True)
    elif lst==3:
            x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK=np.loadtxt(GNS+'field12_on_brick_reduced.txt',unpack=True)
if chip == 2:
     x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK=np.loadtxt(GNS+'field12_on_brick.txt',unpack=True)
#%%       

print(len(x_gns))

valid=np.where((mH<90) & (mK<90) )
x_gns=x_gns[valid]
y_gns=y_gns[valid]
mH=mH[valid]
mK=mK[valid]

#selecting only background stars for the alignemnt? uncoment
###################
H_Ks=np.where((mH-mK)>1.3)
x_gns=x_gns[H_Ks]
y_gns=y_gns[H_Ks]
####################

x_gns=x_gns*0.5
y_gns=y_gns*0.5
GNS_xy=np.array([x_gns,y_gns]).T
print('Elements in GNS: %s'%(len(GNS_xy)))

#%%
a ,d , m, dm, f, df,x,y,dx,dy= np.loadtxt(tmp+'BRICK_stars_calibrated_H_chip'+str(chip)+'_sirius.txt',unpack=True)
zoc= np.loadtxt(tmp+'BRICK_stars_calibrated_H_chip'+str(chip)+'_sirius.txt')
zoca=np.array([x,y]).T
print('Elements in Zoc: %s'%(len(x)))
#%%

 # -638.350, 215.718 manually stimated xoff yoff
check_x,check_y=2,2
while abs(check_x) >1 or abs(check_y)>1  :
    m,(_,_)= aa.find_transform(zoca, GNS_xy,max_control_points=250)
    print('For chip%s'%(3)+'\n'+"Translation: (x, y) = (%.2f, %.2f)"%(m.translation[0],m.translation[1]))
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
print('___NOW TRANSFORMING GNS___')  
trans=aa.matrix_transform(zoca, m.params)
zoc[:,6]=trans[:,0]
zoc[:,7]=trans[:,1]
np.savetxt(tmp+'aa_BRICK_stars_calibrated_H_chip'+str(chip)+'_sirius.txt',zoc, fmt='%.5f') 







