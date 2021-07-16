#!/usr/bin/env python
# coding: utf-8

# In[1]:



import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import distance
from astropy.stats import sigma_clip
from astropy.stats import sigma_clipped_stats
from scipy.stats import gaussian_kde
import sys

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
tmp='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/054_'+band+'/dit_'+str(exptime)+'/'+folder+'tmp/'
sirius='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/SIRIUS/'
GNS='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/field12/'
scripts='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/scripts/'

stream = open(scripts+'polywarp.py')
read_file = stream.read()
exec(read_file)

distancia=1
chip=3# in this case only chip 2 and chip 3 have common elements with the GNS on the brick
unc=0.01

# In[3]:




field12=np.loadtxt(GNS+'field12_on_brick.txt')

h_ks=field12[:,10]-field12[:,12]
field12=np.c_[field12,h_ks]
in_place=np.where((h_ks>1.3))
field12=field12[in_place]
field12[:,0]*=0.5
field12[:,2]*=0.5
np.savetxt(GNS+'field12_no_foreground.txt',field12,fmt='%.6f',header='x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK, H-Ks,')

# x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK,H_Ks=np.loadtxt(GNS+'field12_no_foreground.txt',unpack=True)
#sys.exit("STOP")
# Here we are to select GNS stars by their uncertainty.
dxy_gns=np.sqrt((field12[:,1]*0.106)**2+(field12[:,3]*0.106)**2)
low_g=np.where(dxy_gns<unc)
field12=field12[low_g]
x_gns=field12[:,0]
y_gns=field12[:,2]

gns_xy=np.array([x_gns,y_gns]).T

for chip in range(chip,chip+1):
    #a ,d , m, dm, f, df,x,y,dx,dy=np.loadtxt(tmp+'stars_calibrated_'+band+'_chip'+strn(chip)+'_sirius.txt',unpack=True)
    brick=np.loadtxt(tmp+'stars_calibrated_'+band+'_chip'+str(chip)+'_sirius.txt')
    # Here we are to select ZOC stars by their uncertainty.
    dxy_zoc=np.sqrt(brick[:,8]**2+brick[:,9]**2)
    low_z=np.where(dxy_zoc<unc)
    brick=brick[low_z]
    #We have to add the coordinates offset between the two lists
    #I this case I have choose this bu coomparing in Aladin
    if band=='H' and chip==3:
        ##########################################################
   
        xm_ref,ym_ref=264, 637# xm_ref is GNS
        xm,    ym    =901,425
 
        xoff = xm_ref - xm
        yoff = ym_ref - ym
        ##########################################################
    elif band=='H' and chip==2:
        ##########################################################
       
         


      
     
        xm_ref,ym_ref=332.838000 ,244.081000 # xm_ref is GNS
        xm,ym        = 967.8858643, 2221.4008789

    
        
    
        xoff = xm_ref - xm
        yoff = ym_ref - ym
        ##########################################################
    
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
   
    ciclo=20
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
        ciclo+=10
    ##########################################################
    gns_txt=[diff[i][0][0:] for i in range(len(diff))]
    zoc_txt=[diff[i][1][0:] for i in range(len(diff))]
    x_dis=[diff[i][1][6]-diff[i][0][0] for i in range(len(diff))]
    y_dis=[diff[i][1][7]-diff[i][0][2] for i in range(len(diff))]
    displa=np.array([x_dis,y_dis]).T
    zoc_txt=np.c_[zoc_txt,displa]
    np.savetxt(tmp+'GNS_commons_w_Zoc_c%s.txt'%(chip),gns_txt,header='x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK, H-Ks')
    np.savetxt(tmp+'Zoc_c%s_commons_w_GNS.txt'%(chip),zoc_txt,header='a ,d , m, dm, f, df,x,y,dx,dy,x_dis,y_dis. X and Y are the correspondig coorinates wit GNS, They are not the original ones!!!!')
    np.savetxt(tmp+'dis_xy_chip%s.txt'%(chip),displa,header='Displacement in pixels.')

    # In[5]:
    
    
    
    diff=[]
    for i in range(len(gns_xy)): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
        dist=distance.cdist(gns_xy[i:i+1,0:2],brick[:,6:8], 'euclidean')
        d=np.where(dist<distancia)
        if len(d[1])>0:
            diff.append((gns_xy[i],brick[d[1][np.argmin(dist[d])]]))
    diff=np.array(diff)
    l_mag=1
    
    x_shift=[(diff[t][0][0]-diff[t][1][6]) for t in range(len(diff))]# if abs(diff[t][0][brillo]-diff[t][1][2])<l_mag]
    y_shift=[(diff[t][0][1]-diff[t][1][7]) for t in range(len(diff))]# if abs(diff[t][0][brillo]-diff[t][1][2])<l_mag]
    print(len(x_shift))

    x_shift=np.array(x_shift)*0.106
    y_shift=np.array(y_shift)*0.106
    fig,ax= plt.subplots(1,2,figsize=(20,10))
    ls=[x_shift,y_shift]
    nam=['x diff (arcsec)','y diff (arcsec)']
    #ax[0]=plt.suptitle('Sigma Threshold at reconstruct = %s, CHIP  %s'%(th,ch),fontsize=20)
    #plt.clf()
    for h in range(len(ls)):
        sig_h=sigma_clipped_stats(ls[h],sigma=3,maxiters=20,cenfunc='mean')
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
    
    # In[10]:
    
    
    
   


    
    
    
    # In[6]:
    
    '''   
    from scipy.interpolate import fitpack
    from scipy import interpolate
    #x=[(diff[t][0][-2]) for t in range(len(diff))]
    #y=[(diff[t][0][-1]) for t in range(len(diff))]
    
    z=np.sqrt(np.array(x_shift)**2+np.array(y_shift)**2)
    x=[(diff[t][0][0]) for t in range(len(diff))]
    y=[(diff[t][0][1]) for t in range(len(diff))]
    xline = np.linspace(min(x), max(x), 2000)
    yline = np.linspace(min(y), max(y), 2000)
    #xline = np.linspace(0, 2548, 2549)
    #yline = np.linspace(0, 2548, 2549)
    xgrid,ygrid = np.meshgrid(xline, yline)
    #tck=fitpack.bisplrep(xs,ys,zs,s=len(ys)-2*np.sqrt(len(ys)))
    tck=fitpack.bisplrep(x,y,z,kx=3,ky=3)
    #znew = interpolate.bisplev(xnew, ynew, tck)
    znew = interpolate.bisplev(xline, yline, tck)
    
    plt.figure()
    #plt.pcolormesh(xnew, ynew, znew, shading='flat',vmin=min(zs), vmax=max(zs))
    #plt.pcolormesh(xnew, ynew, znew)
    plt.pcolormesh(xgrid, ygrid, znew)
    plt.colorbar()
    plt.title("Difference on position for the alignment")
    plt.xlabel("X diff (pixel)")
    plt.ylabel("Y diff (pixel)")
    plt.show()
    '''   
    
    # In[7]:
    
    
    #x_gns, dx_gns, y_gns, dy_gns, raH, draH, decH, ddecH, mJ, dmJ, mH, dmH, mK, dmK=np.loadtxt(GNS+'field12_on_brick.txt',unpack=True)
    #field12=np.loadtxt(GNS+'field12_on_brick.txt')
    #a ,d , m, dm, f, df,x,y,dx,dy=np.loadtxt(tmp+'stars_calibrated_'+band+'_chip'+strn(chip)+'_sirius.txt',unpack=True)
    #brick=np.loadtxt(tmp+'stars_calibrated_'+band+'_chip'+str(chip)+'_sirius.txt')
    sig=2

    diff=[]
    for i in range(len(gns_xy)): #compara las distancia entre los puntos y guarda las menores que a, si hay mas de dos puntos con distancias menores que a, guarda la más perqueña
        dist=distance.cdist(gns_xy[i:i+1,0:2],brick[:,6:8], 'euclidean')
        d=np.where(dist<distancia)
        if len(d[1])>0:
            diff.append((field12[i],brick[d[1][np.argmin(dist[d])]]))
    
    resta=[]   
    mags=[]
    for i in range(len(diff)):
        #print(diff[i][0][2]-diff[i][1][2])
        resta.append(diff[i][0][10]-diff[i][1][2])
        mags.append((diff[i][0][10],diff[i][1][2]))
        
    diff=np.array(diff)
    resta=np.array(resta)
    average=np.mean(np.array(resta))
    mags=np.array(mags)
    
    s=sigma_clipped_stats(resta,sigma=sig,maxiters=10)
    mask_sig=sigma_clip(resta,sigma=sig,maxiters=10)
    nope=np.where(mask_sig.mask==True)
    
    nbins=(np.ceil(max(mags[:,1]))-np.floor(min(mags[:,1])))
    mags_bins=np.floor(min(mags[1]))+np.arange(nbins)
    
    fig,ax=plt.subplots(1,1,figsize=(20,10))
    ax.scatter(mags[:,1],resta,color='k',alpha=0.2)
    ma=mags[:,1]
    ax.scatter(mags[:,1][nope],resta[nope],color='red',marker='x',s=100,alpha=0.7)
    ax.grid()
    ax.axhline(s[0],color='g',ls='--',lw=2)
    ax.set_ylim(-2,2)
    ax.tick_params(axis='x', labelsize=20)
    ax.tick_params(axis='y', labelsize=20)
    ax.legend(['Chip %s,Common stars #%s'%(chip,len(mags[:,1]))],fontsize=20,shadow=True,loc=1, handlelength=0,handletextpad=0)
    ax.set_xlabel('$[%s]_{Zoc}$'%(band),fontsize=20)
    ax.set_ylabel('$[%s]_{Zoc}-[%s]_{SIR}$'%(band,band),fontsize=20)
    #ax[i].set_xlim(x0,x1)
    mmag =np.zeros(shape=(int(nbins)))
    sig_mag =np.zeros(shape=(int(nbins)))
    for j in range(int(nbins)):
        thisbin=np.where((mags[:,1]>mags_bins[j])&(mags[:,1]<=mags_bins[j]+1))
        vals = resta[thisbin]
        nope_thisbin=np.where((mags[:,1][nope]>mags_bins[j])&(mags[:,1][nope]<mags_bins[j]+1))
        if len(vals)>1:
            bin_rej=sigma_clip(vals, sigma=2, maxiters=5,masked=True)
            sig_bin=sigma_clipped_stats(vals,sigma=2.0,maxiters=5)
            mmag[j]=sig_bin[0]
            sig_mag[j]=sig_bin[2]
        ax.errorbar(mags_bins[j]+0.5,mmag[j],sig_mag[j],color='red', elinewidth=3,capsize=10,capthick=2,barsabove=True,zorder=3)
        ax.text(min(mags[:,1])+1,1,'mean offset (%s$\sigma$) =%.3f, std= %.3f'%(sig,s[0],s[2]),color='green',fontsize=14) 
        ax.text(mags_bins[j]+0.4,-1.5,'%.3f'%(sig_mag[j]),color='red',fontsize=14)  
        ax.text(mags_bins[j]+0.4,-1.7,'%s'%(len(vals)),color='blue',fontsize=14)
        ax.text(mags_bins[j]+0.4,-1.9,'%s'%(len(nope_thisbin[0])),color='orange',fontsize=14)
    
  
#%%

    dx_cgns=[diff[c][0][5] for c in range(len(diff))]
    dy_cgns=[diff[c][0][7] for c in range(len(diff))]
    
    dx_czoc=[diff[c][1][8] for c in range(len(diff))]
    dy_czoc=[diff[c][1][9] for c in range(len(diff))]    
#%%

    # In[ ]:
    
    
    map_c='inferno'
    uncx=[dx_cgns,dx_czoc]
    uncy=[dy_cgns,dy_czoc]
    le=['GNS','Zoc.']
    
    fig, ax = plt.subplots(2,2,figsize=(15,10))
    
    for g in range(2):
        xy = np.vstack([mags[:,g],uncx[g]])
        xy1 = np.vstack([mags[:,g],uncy[g]])
        z = gaussian_kde(xy)(xy)
        z1= gaussian_kde(xy1)(xy1)
        
        
        #ax.invert_yaxis()
        #ax.set_ylim(16,10)
        ax[g,0].scatter(mags[:,g], uncx[g], c=z, s=5,alpha=1,cmap=map_c)
        ax[g,0].set_xlabel('[%s]'%(band),fontsize=20)
        ax[g,0].grid()
        ax[g,0].set_ylim(0,0.02)
        ax[g,1].set_ylim(0,0.02)
        ax[g,0].set_ylabel('dX (arcsec)',fontsize=20)
        ax[g,1].scatter(mags[:,1], uncy[g], c=z1, s=5,alpha=1,cmap=map_c)
        ax[g,1].set_xlabel('[%s]'%(band),fontsize=20)
        ax[g,1].set_ylabel('dY (arcsec)',fontsize=20)
        ax[g,1].grid()
        ax[g,0].set_ylim(0,0.0125)
        ax[g,1].set_ylim(0,0.0125)
        ax[g,1].legend([le[g]],fontsize=20,shadow=True,loc=1,markerscale=0, handlelength=0,handletextpad=0)
        plt.suptitle('%s common stars to Zocalli (Chip %s) and GNS field 12'%(len(uncx[g]),chip),fontsize=20)
    
