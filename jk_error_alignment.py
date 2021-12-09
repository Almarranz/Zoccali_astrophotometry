import numpy as np
import scipy as sp
import pylab
import matplotlib.pyplot as plt
from scipy.stats import norm
import matplotlib.mlab as mlab
import glob
import os

#new: deg3 alignment, errors of hst from starfinder, cutting also the hst values based on the position errors, and from filter 20

jackknive='/Users/amartinez/Desktop/PhD/HAWK/The_Brick/photometry/jack_knive/'

degree=1
line=[]



it=len(glob.glob(jackknive+'Z1_aa_NPL058_jackknife_degree_%s*'%(degree)))
#%%
for i in range(1,it+1):
		f=open(jackknive+'Z1_aa_NPL058_jackknife_degree_%s_%s.txt' %(degree,i)) ####instead of aaaa name of your file ; dmax=2, degree=1
		line.append(f.readlines())
        
# for i in range(1):	
#     v=np.loadtxt(jackknive+'Z1_aa_NPL058_jackknife_degree_%s_%s.txt' %(degree,i+1), usecols=[0])
#     w=np.loadtxt(jackknive+'Z1_aa_NPL058_jackknife_degree_%s_%s.txt' %(degree,i+1), usecols=[1])
#     plt.plot(v, w , "bo")	

tot=[]
for j in range(len(line[0])):
    for k in range(it): 
        tot.append(float(line[k][j].split()[0])) ### for x  split()[0] but for y split()[1]
#%%


mufitx=[]
sigmafitx=[] 
for u in range(len(line[0])):
# for u in range(1):
    x=  ((tot[u*it:(u+1)*it]))
    # plt.hist(x, bins='auto', density=True)
    (mu, sigma) = norm.fit(x)
	# ~ print "sigma=", sigma
	# ~ print "np.std(x)=", np.std(x)
    mufitx.append(mu)
    sigmafitx.append(sigma)
	# ~ print mu
	# ~ print sigma
    # y = norm.pdf(np.linspace(min(x),max(x),100), mu, sigma)
    # plt.plot(np.linspace(min(x),max(x),100), y, 'r--', linewidth=2)
np.savetxt(jackknive+'mufitx_degree%s.txt'%(degree),mufitx)
np.savetxt(jackknive+'sigmitx_degree%s.txt'%(degree),sigmafitx)


#%%
tot=[]
for j in range(len(line[0])):
    for k in range(it): 
        tot.append(float(line[k][j].split()[1])) ### for x  split()[0] but for y split()[1]
        
mufity=[]
sigmafity=[] 
for u in range(len(line[0])):
# for u in range(1):
    y=  ((tot[u*it:(u+1)*it]))
    # plt.hist(y, bins='auto', density=True)
    (mu, sigma) = norm.fit(y)
	# ~ print "sigma=", sigma
	# ~ print "np.std(x)=", np.std(x)
    mufity.append(mu)
    sigmafity.append(sigma)
	# ~ print mu
	# ~ print sigma
    # yprim = norm.pdf(np.linspace(min(y),max(y),100), mu, sigma)
    # plt.plot(np.linspace(min(x),max(x),100), yprim, 'r--', linewidth=2)
np.savetxt(jackknive+'mufity_degree%s.txt'%(degree),mufity)
np.savetxt(jackknive+'sigmity_degree%s.txt'%(degree),sigmafity)



   
#%%
sigmafitx = np.array(sigmafitx)
sigmafitx_as=sigmafitx*0.106
meanx=np.mean(sigmafitx_as)
print(np.mean(sigmafitx_as))

sigmafity = np.array(sigmafity)
sigmafity_as=sigmafity*0.106
meany=np.mean(sigmafity_as)
print(np.mean(sigmafity_as))

print('quadratic : %s'%(np.sqrt(meanx**2+meany**2)))
#%%
print(np.mean(sigmafitx_as)*1000)
print(np.mean(sigmafity_as)*1000)