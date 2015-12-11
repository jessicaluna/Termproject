# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 00:17:37 2015

@author: jessicaluna
"""

from pylab import*
import numpy as np
import matplotlib.pyplot as plt
import time
import corner
import math
from matplotlib import gridspec

Msun= 1.989*10**(33)
G=6.67*10**(-8)
Mjup=1.899e30


#systemname='HD16760'
#systemname='HD80606'
#systemname='16cyg'
systemname='HD209458'
rvfile='rvdata'+systemname+'.txt'
startvals=systemname+ 'initialvals.txt'
litvals=systemname+'litvals.txt'
MCMCsamples=systemname+'MCMC.txt'
medvalsfile=systemname+'Medvals.txt'
Mhist=systemname+'Mshist.pdf'
Phist=systemname+'Pshist.pdf'
Tohist=systemname+'Tohist.pdf'
ehist=systemname+'ehist.pdf'
whist=systemname+'whist.pdf'
histograms=systemname+'histogram.pdf'
rvphasedplot=systemname+'rvphased.pdf'
rvplot=systemname+'rv.pdf'

JDdate,RV_obs, RV_obs_err=loadtxt(rvfile,unpack=True,usecols=[0,1,2])
#RV_obs=RV_obs+183.0 #for HD 80606 the data is offset
Ms,Ps,Tos,es,ws,chi2=loadtxt(MCMCsamples,unpack=True,usecols=[0,1,2,3,4,5])
Msinia,Pa,Toa,ea,wa,Mstar=loadtxt(litvals,unpack=True,usecols=[0,1,2,3,4,5])

Mstar=Mstar*Msun  #HD 16760
Msinia=Msinia*Mjup

Msini,P,To,e,w=[np.median(x) for x in Ms,Ps,Tos,es,ws]
#print(Msini,P,To,e,w)
savetxt(medvalsfile, (Msini/Mjup,P,To,e,w), fmt='%.18e', delimiter=' ', newline='\n', header='medvals of Msini(mjup),P,Tzero,e,w')


plt.hist(Ms/Mjup,label='Msini',bins=15)
plt.title('Companion Mass')
plt.xlabel('Msini (Mjup)')
plt.savefig(Mhist)
plt.clf()

plt.hist(Ps,label='P',bins=15)
plt.title('Period')
plt.xlabel('P (days)')
plt.savefig(Phist)
plt.clf()

plt.hist(Tos,label='To',bins=15)
plt.title('Time of Periastron')
plt.xlabel('To (JD days)')
plt.savefig(Tohist)
plt.clf()

plt.hist(es,label='e',bins=15)
plt.title('eccentricity')
plt.xlabel('e')
plt.savefig(ehist)
plt.clf()

plt.hist(ws,label='w',bins=15)
plt.title('Argument of Periapsis')
plt.xlabel('w (deg)')
plt.savefig(whist)
plt.clf()


def calcrv2(Msini,P,To,e,w,dates):
    w=radians(w)
    Mm=(2.0*pi)*(dates-To)/P
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
    return rv  

def calcrvph(Msini,P,To,e,w,dates):
    w=radians(w)
    Mm=(2.0*pi)*(dates)/P
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
    return rv
 


   
t=linspace(0,Pa,2000)
times=linspace(min(JDdate)-100,max(JDdate)+100,5000)


rvactualph= calcrvph(Msinia,Pa,Toa,ea,wa,t)
rvactual2= calcrv2(Msinia,Pa,Toa,ea,wa,times)
rvmodelph= calcrvph(Msini,P,To,e,w,t)
rvmodel2=calcrv2(Msini,P,To,e,w,times)

phase =(JDdate - Toa)/ Pa %1.0


##phased rv curve
plt.errorbar(phase, RV_obs, yerr=RV_obs_err,marker='.',ls='none',label='data')
plt.plot(t/Pa,rvmodelph*1.0e-2,c='r',label='Model')
plt.plot(t/Pa,rvactualph*1.0e-2,c='g',label='literature')
plt.ylabel('RV m/s')
plt.xlabel('Phase')
legend(loc='best')
plt.savefig(rvphasedplot)
plt.clf()

#nonphased plot
plt.errorbar(JDdate, RV_obs, yerr=RV_obs_err,marker='.',ls='none',label='data')
plt.plot(times,rvmodel2*1.0e-2,c='r',label='Model')
plt.plot(times,rvactual2*1.0e-2,c='g',label='literature')
plt.ylabel('RV m/s')
plt.xlabel('JD dates')
legend(loc='best')
plt.savefig(rvplot)
plt.clf()
    
