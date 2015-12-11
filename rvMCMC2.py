# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 21:24:31 2015

@author: jessicaluna
"""

from pylab import*
import numpy as np
import matplotlib.pyplot as plt
import time

G=6.67*10**(-8)
Mjup=1.899e30
Msun= 1.989*10**(33)

starttime=time.time()


##change file names to appropriate system###
#systemname='HD16760'
#systemname='HD80606'
#systemname='16cyg'
systemname='HD209458'

rvfile='rvdata'+systemname+'.txt'
startvals=systemname+'initialvals.txt'
mcmcfilename=systemname+'MCMC.txt'


JDdate,RV_obs, RV_obs_err=loadtxt(rvfile,unpack=True,usecols=[0,1,2])
#RV_obs=RV_obs+183.0 #for HD 80606 the data is offset
Msini,P,To,e,w,Mstep,Pstep,Tstep,estep,wstep,Mstar=loadtxt(startvals,unpack=True,usecols=[0,1,2,3,4,5,6,7,8,9,10])


##starting values

RV_obs=RV_obs*1e2 #cm/s
RV_obs_err=RV_obs_err*1e2 #cm/s


Mstar=Mstar*Msun  #HD 16760
Msini=Msini*Mjup
Mstep=Mstep*Mjup
 
def calcrv(Msini,P,To,e,w):
    w=radians(w)
    Mm=(2.0*pi)*(JDdate-To)/P
    #setting initial values for trancendental equation
    Eo=Mm+e*sin(Mm)
    g=Eo-e*sin(Eo)-Mm
    gprime=1-e*cos(Eo)
    #solving trancendental equation
    for j in arange(1,30):
        E=Eo-g/gprime
        Eo=E
        g=Eo-e*sin(Eo)-Mm
        gprime=1-e*cos(Eo)    
    E=Eo
    f=2.0*np.arctan(sqrt((1.0+e)/(1.0-e))*tan(E/2.0))   
    K=(Msini/sqrt(1.0-e**2))*(2*pi*G/(P*86400.0*Mstar**2))**(1.0/3.0)
    rv=K*(cos(f+w)+e*cos(w))
    X2 = sum(((RV_obs-rv)/RV_obs_err)**2)
    return X2
    
niters=300000   
nburners=5000
accept=[]
Msamples= []
Psamples= []
Tsamples= []
esamples= []
wsamples= []
chi2vals=[]



naccept=0
nreject=0
chicurrent=calcrv(Msini,P,To,e,w)


for jj in range(niters+nburners):
#for jj in range(niters):
    #calculate trial model
    num=np.random.random_integers(1,high=5)    
    
    if num==1:
        Mtrial=np.random.normal(Msini,Mstep)
        Ptrial=P
        Ttrial=To
        etrial=e
        wtrial=w
    elif num==2:
        Mtrial=Msini
        Ptrial=np.random.normal(P,Pstep)
        Ttrial=To
        etrial=e
        wtrial=w
    elif num==3:
        Mtrial=Msini
        Ptrial=P
        Ttrial=np.random.normal(To,Tstep)
        etrial=e
        wtrial=w
    elif num==4:
        Mtrial=Msini
        Ptrial=P
        Ttrial=To 
        etrial=np.random.normal(e,estep)
        wtrial=w
    elif num==5:
        Mtrial=Msini
        Ptrial=P
        Ttrial=To
        etrial=e
        wtrial=np.random.normal(w,wstep)   
    
    chitrial=calcrv(Mtrial,Ptrial,Ttrial,etrial,wtrial)
    
    ratio=exp(-0.5*(chitrial-chicurrent))
    
    alpha=min(1,ratio)
    u=np.random.uniform(0.0,1.0)
    if u <=alpha:
        Msini=Mtrial
        P=Ptrial
        To=Ttrial
        e=etrial
        w=wtrial
        chicurrent=chitrial
        if nburners<jj:
            accept.append(1)
            naccept +=1
        else:
            accept.append(0)
            nreject +=0

    if nburners<jj:
        chi2vals.append(chicurrent)
        Msamples.append(Msini)
        Psamples.append(P)
        Tsamples.append(To)
        esamples.append(e)
        wsamples.append(w)
    #print('Chi trial',chitrial)    
        


print(float(naccept)/niters, naccept,nreject)
print('Chisq',chicurrent)
print(Msini/Mjup,P,To,e,w)

savetxt(mcmcfilename, zip(Msamples,Psamples,Tsamples,esamples,wsamples,chi2vals), fmt='%.18e', delimiter=' ', newline='\n', header='Msini,P,Tzero,e,w.chi2vals')




print(time.time()-starttime)





